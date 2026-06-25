"""Slide-tags per-cell mapping: DBSCAN (raw + norm weightings) + σ topographic
prominence + spatial profiles, in a single pass per cell. Merges the former
classify.py (data collection) and cell_profiles.py (spatial profiles) — one
DBSCAN and one KDE surface per cell, no intermediate file.

Collects raw per-cell measurements only — it assigns NO status (no
singlet/doublet/unmapped call); thresholding is left to downstream consumers.

Per-cell algorithm
  1. DBSCAN at fixed eps, run under TWO weightings (both emitted, prefixed
     raw_/norm_):
       raw  — sample_weight = raw per-bead UMI; descending min_samples
              escalation, first tier yielding ≥1 cluster.
       norm — sample_weight = u/∑u · norm_to; integer min_samples scanned
              high→low (50→1 by default), highest tier yielding ≥2 clusters.
     In each, the largest N_DBSCAN_CLUSTERS clusters (c1..c8) by ∑weight are
     recorded individually, a c9plus pool holds any beyond that, c0 = noise.
     cx1/cy1/cx2/cy2 alias the raw-mode c1/c2 centroids.
  2. Fixed-σ Gaussian-sum KDE at σ µm on a G×G grid spanning the cell's bbox
     (+ 3σ padding).  Surface units: UMI / µm².
  3. Topographic prominence via descending-flood union-find. Records the top
     N_TOP_PEAKS peak heights/prominences, the full top-`nprom` prominence and
     height rankings, and background stats. Each mode's c1 centroid maps to its
     nearest KDE peak.
  4. Leftmost-line obs-vs-projected log10 fold changes for peaks 1 & 2:
       peak{1,2}fclr    — fit log10(prom) vs rank from rank 3, prom > 0.1 excluded.
       peak{1,2}fclr_h  — same on raw height (peaks re-ranked by height), prom
                          > 0.1 still excluded.
     fclr/fclr_h are each fit in linear-rank and log10-rank x-space (`_logr`),
     recording R², n points fit, window start/stop rank, and which fit better.
  5. Cumulative UMI vs radius from the raw c1 centroid (r = 1..400 µm).

Inputs
  - Pair whitelist with columns ``cell_bc_10x``, ``bead_bc``, ``nUMI`` — the
    pair step's cell-barcode_coords.csv, which also carries ``x_um``/``y_um``
    (puck merged in pair). Coords are read straight from it.
  - Puck CSV (no header): ``bead_bc, x_um, y_um`` — only needed as a fallback
    for a legacy coord-free whitelist (df_whitelist.txt).
  (Near-)homopolymer beads (a single base in ≥12/14 positions) and beads with no
  coords are dropped automatically before clustering/KDE.

Output (--out-dir)
  - cell_coords.csv : one row per cell with all of the above (comma-separated).

CLI
  python map_cells.py --sample S --whitelist cell-barcode_coords.csv \\
      [--puck puck.csv] --out-dir OUT
  Tunables: --eps, --min-samples, --norm-to, --norm-ms-scan, --sigma, --grid,
            --nprom, --radii, --ncores.
"""
import argparse, os, time
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from joblib import Parallel, delayed

# ── tunables ────────────────────────────────────────────────────────────────
N_TOP_PEAKS = 8           # number of peakK_* columns emitted
N_DBSCAN_CLUSTERS = 8     # number of per-cluster column sets emitted (c1..c8)
MODES = ("raw", "norm")

DEFAULT_EPS = [50]
DEFAULT_MIN_SAMPLES = [10, 9, 8]      # raw-mode escalation (descending)
DEFAULT_NORM_TO = 1000               # norm-mode weights: w = u/∑u * norm_to
DEFAULT_NORM_MS_HI = 50               # norm-mode min_samples scan: high → low
DEFAULT_NORM_MS_LO = 1
DEFAULT_SIGMA = 40.0
DEFAULT_GRID = 120
DEFAULT_NPROM = 1000
DEFAULT_RADII = list(range(1, 401, 1))
DEFAULT_NCORES = 16

LO_FIT, WIN_FIT, CUT_FIT = 3, 4, 0.80   # leftmost-line: peaks 1-2 excluded, 4-pt R²>0.8
FIT_PROM_MAX_LR = 0.1                     # fclr/fclr_h fit: peaks above this prom excluded

PARAMS = dict()   # filled at runtime from CLI

NB = np.array([[-1,-1],[-1,0],[-1,1],[0,-1],[0,1],[1,-1],[1,0],[1,1]])


# ── KDE + topographic prominence ────────────────────────────────────────────
def topo_prominence(Z):
    """Descending-flood union-find topographic prominence.
    Returns (prom_dict: {flat_peak_idx: prom}, floor_height)."""
    nr, nc = Z.shape
    flat = Z.ravel()
    order = np.argsort(-flat, kind="stable")
    parent = -np.ones(nr*nc, dtype=np.int64)
    peak = -np.ones(nr*nc, dtype=np.int64)
    prom = {}
    seen = np.zeros(nr*nc, dtype=bool)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for idx in order:
        seen[idx] = True
        r, c = divmod(idx, nc)
        nbr_comps = set()
        for dr, dc in NB:
            rr, cc = r+dr, c+dc
            if 0 <= rr < nr and 0 <= cc < nc:
                ni = rr*nc + cc
                if seen[ni]: nbr_comps.add(find(ni))
        if not nbr_comps:
            parent[idx] = idx; peak[idx] = idx
        elif len(nbr_comps) == 1:
            parent[idx] = next(iter(nbr_comps))
        else:
            comps = list(nbr_comps)
            heights = [flat[peak[c]] for c in comps]
            tallest = comps[int(np.argmax(heights))]
            cur_val = flat[idx]
            for c in comps:
                if c == tallest: continue
                p = peak[c]
                prom[p] = flat[p] - cur_val
                parent[c] = tallest; peak[c] = -1
            parent[idx] = tallest
    floor = float(flat.min())
    for c in range(nr*nc):
        if parent[c] == c and peak[c] != -1:
            p = peak[c]; prom[p] = flat[p] - floor
    return prom, floor


def gauss_sum(qx, qy, x, y, u, sig):
    """∑ u_i · N(q − p_i; σI) = (∑ u_i exp(-d²/(2σ²))) / (2πσ²) → UMI/µm²."""
    inv2s2 = 1.0 / (2*sig*sig)
    norm = 1.0 / (2*np.pi * sig*sig)
    qx_r = qx.ravel(); qy_r = qy.ravel()
    out = np.zeros_like(qx_r)
    for xi, yi, ui in zip(x, y, u):
        out += ui * np.exp(-((qx_r-xi)**2 + (qy_r-yi)**2) * inv2s2)
    return (out * norm).reshape(qx.shape)


def kde_peaks(x, y, u):
    """σ KDE on a G×G grid (bbox + 3σ pad) → interior topo peaks (px, py, heights,
    proms, descending by prominence) + background stats over interior pixels."""
    sigma = PARAMS["sigma"]; G = PARAMS["grid"]
    pad = 3 * sigma
    xlim = (x.min() - pad, x.max() + pad)
    ylim = (y.min() - pad, y.max() + pad)
    gx, gy = np.mgrid[xlim[0]:xlim[1]:G*1j, ylim[0]:ylim[1]:G*1j]
    Z = gauss_sum(gx, gy, x, y, u, sigma)
    EDGE = 2 * sigma
    interior = ((gx >= xlim[0] + EDGE) & (gx <= xlim[1] - EDGE)
              & (gy >= ylim[0] + EDGE) & (gy <= ylim[1] - EDGE))
    Z_int = Z[interior] if interior.any() else Z
    bg_median = float(np.median(Z_int))
    bg_p05 = float(np.quantile(Z_int, 0.05))
    bg_p95 = float(np.quantile(Z_int, 0.95))
    prom_all, _ = topo_prominence(Z)
    empty = dict(px=np.array([]), py=np.array([]),
                 heights=np.array([]), proms=np.array([]),
                 bg_median=bg_median, bg_p05=bg_p05, bg_p95=bg_p95)
    keys_all = np.array(list(prom_all.keys()))
    if len(keys_all) == 0:
        return empty
    keep = interior.ravel()[keys_all]
    keys = keys_all[keep]
    if len(keys) == 0:
        return empty
    proms = np.array([prom_all[k] for k in keys])
    heights = Z.ravel()[keys]
    order = np.argsort(-proms)
    keys = keys[order]; proms = proms[order]; heights = heights[order]
    px = gx.ravel()[keys]; py = gy.ravel()[keys]
    return dict(px=px, py=py, heights=heights, proms=proms,
                bg_median=bg_median, bg_p05=bg_p05, bg_p95=bg_p95)


def _peak_cols(heights, proms, bg_median, n_top=N_TOP_PEAKS):
    """Per-peak columns: peakK_h, peakK_prom, peakK_prom_normbg (k=1..n_top)."""
    out = {}
    nb = max(bg_median, 1e-12)
    for k in range(n_top):
        if k < len(heights):
            out[f"peak{k+1}_h"] = float(heights[k])
            out[f"peak{k+1}_prom"] = float(proms[k])
            out[f"peak{k+1}_prom_normbg"] = float(proms[k] / nb)
        else:
            out[f"peak{k+1}_h"] = np.nan
            out[f"peak{k+1}_prom"] = np.nan
            out[f"peak{k+1}_prom_normbg"] = np.nan
    return out


# ── DBSCAN (raw + norm weightings) ──────────────────────────────────────────
def _dbscan_precomp(D, w, eps, min_samples):
    return DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed"
                  ).fit(D, sample_weight=w).labels_


def _clusters_from_labels(labels, x, y, u, w):
    """The largest N_DBSCAN_CLUSTERS clusters by ∑w (each dict with n_beads, umi,
    cx, cy), a c9plus pool for any beyond, c0 (noise, label -1), and n_clusters.
    c1..c8 + c9plus + c0 account for every bead/UMI in the cell."""
    uniq = [l for l in np.unique(labels) if l != -1]
    sizes = {l: float(w[labels == l].sum()) for l in uniq}
    ordered = sorted(sizes, key=sizes.get, reverse=True)
    def mk(lab):
        m = (labels == lab)
        return dict(n_beads=int(m.sum()), umi=float(u[m].sum()),
                    cx=float(x[m].mean()), cy=float(y[m].mean()))
    clusters = [mk(lab) for lab in ordered[:N_DBSCAN_CLUSTERS]]
    rest = ordered[N_DBSCAN_CLUSTERS:]
    if rest:
        mr = np.isin(labels, rest)
        crest = dict(n_beads=int(mr.sum()), umi=float(u[mr].sum()), n_clusters=len(rest))
    else:
        crest = dict(n_beads=0, umi=0.0, n_clusters=0)
    m0 = (labels == -1)
    c0 = dict(n_beads=int(m0.sum()), umi=float(u[m0].sum()))
    return clusters, crest, c0, len(uniq)


def _empty_db(n_pts, u):
    return dict(clusters=[], crest=dict(n_beads=0, umi=0.0, n_clusters=0),
                c0=dict(n_beads=int(n_pts), umi=float(u.sum())),
                n_clusters=0, min_samples_used=np.nan)


def _dbscan_escalate(D, x, y, u, w, eps, ms_steps, n_pts):
    """Raw mode: descending min_samples escalation — first tier with ≥1 cluster wins."""
    if n_pts < min(ms_steps):
        return _empty_db(n_pts, u)
    for ms in ms_steps:
        lbl = _dbscan_precomp(D, w, eps, ms)
        if any(l != -1 for l in np.unique(lbl)):
            clusters, crest, c0, nclu = _clusters_from_labels(lbl, x, y, u, w)
            return dict(clusters=clusters, crest=crest, c0=c0, n_clusters=nclu,
                        min_samples_used=float(ms))
    return _empty_db(n_pts, u)


def _dbscan_scan_2clusters(D, x, y, u, w, eps, ms_hi, ms_lo, n_pts):
    """Norm mode: scan integer min_samples from ms_hi down to ms_lo; take the
    highest that yields ≥2 clusters. NaN min_samples_used if none does."""
    for ms in range(int(ms_hi), int(ms_lo) - 1, -1):
        if ms < 1:
            break
        lbl = _dbscan_precomp(D, w, eps, ms)
        uniq = [l for l in np.unique(lbl) if l != -1]
        if len(uniq) >= 2:
            clusters, crest, c0, nclu = _clusters_from_labels(lbl, x, y, u, w)
            return dict(clusters=clusters, crest=crest, c0=c0, n_clusters=nclu,
                        min_samples_used=float(ms))
    return _empty_db(n_pts, u)


def dbscan_modes(x, y, u):
    """DBSCAN at fixed eps under raw + norm weightings → {'raw':…, 'norm':…},
    each dict(clusters, crest, c0, n_clusters, min_samples_used)."""
    eps = PARAMS["eps"][0]
    ms_steps = PARAMS["min_samples"]
    n_pts = len(x)
    if u.sum() <= 0 or n_pts < 2:
        return dict(raw=_empty_db(n_pts, u), norm=_empty_db(n_pts, u))
    xy = np.c_[x, y]
    dx = xy[:, None, 0] - xy[None, :, 0]; dy = xy[:, None, 1] - xy[None, :, 1]
    D = np.sqrt(dx*dx + dy*dy)
    raw = _dbscan_escalate(D, x, y, u, u.copy(), eps, ms_steps, n_pts)
    w_norm = u / u.sum() * PARAMS["norm_to"]
    norm = _dbscan_scan_2clusters(D, x, y, u, w_norm, eps,
                                  PARAMS["norm_ms_hi"], PARAMS["norm_ms_lo"], n_pts)
    return dict(raw=raw, norm=norm)


# ── leftmost-line fold change ───────────────────────────────────────────────
def find_section(x, y, start, win, cut):
    n = len(x)
    for s in range(start, n - win + 1):
        xx = x[s:s+win]; yy = y[s:s+win]
        m, b = np.polyfit(xx, yy, 1); pred = m*xx + b
        ss = float(((yy-pred)**2).sum()); sst = float(((yy-yy.mean())**2).sum())
        r2 = 1 - ss/sst if sst > 0 else 1.0
        if r2 > cut:
            return s, s+win-1, float(m), float(b), float(r2)
    return None


_FC_NAN = dict(fc1=np.nan, fc2=np.nan, sec_lo=np.nan, sec_hi=np.nan, r2=np.nan, npts=0)


def _fit_leftmost(val, prom, prom_max, xspace, lo_fit=LO_FIT, win=WIN_FIT, cut=CUT_FIT):
    """One leftmost-line obs-minus-projected log10 fold change for peaks 1 & 2.

    `val`  = per-peak quantity to fit (prominence or raw height), descending in its
             own rank order (val[0]/val[1] are peaks 1 & 2 of that ranking).
    `prom` = prominence aligned to `val`'s order; points with prom > `prom_max` are
             excluded from the fit (peaks 1..lo_fit-1 excluded by rank too).
    `xspace` = 'lin' fits log10(val) vs rank; 'log' fits log10(val) vs log10(rank).
    Returns dict(fc1, fc2, sec_lo, sec_hi, r2, npts); all-NaN if no qualifying window."""
    if len(val) < lo_fit + win - 1:
        return dict(_FC_NAN)
    rk = np.arange(lo_fit, len(val)+1, dtype=float)
    yl = np.log10(val[lo_fit-1:])
    fit_ok = prom[lo_fit-1:] <= prom_max
    rk_fit = rk[fit_ok]; yl_fit = yl[fit_ok]
    if len(rk_fit) < win:
        return dict(_FC_NAN)
    xr = (lambda r: np.log10(r)) if xspace == "log" else (lambda r: r)
    s = find_section(xr(rk_fit), yl_fit, 0, win, cut)
    if s is None:
        return dict(_FC_NAN)
    i0, i1, m, b, r2 = s
    fc1 = float(np.log10(val[0]) - (m*xr(1.0) + b))
    fc2 = float(np.log10(val[1]) - (m*xr(2.0) + b)) if len(val) >= 2 else np.nan
    return dict(fc1=fc1, fc2=fc2, sec_lo=float(rk_fit[i0]), sec_hi=float(rk_fit[i1]),
                r2=float(r2), npts=int(win))


def leftmost_peak_fc_both(val, prom, prom_max, lo_fit=LO_FIT, win=WIN_FIT, cut=CUT_FIT):
    """Run the leftmost-line fit in both rank x-spaces. Returns (lin, log, better)
    where lin/log are _fit_leftmost dicts and `better` is 'lin'/'log'/'' (higher R²)."""
    lin = _fit_leftmost(val, prom, prom_max, "lin", lo_fit, win, cut)
    log = _fit_leftmost(val, prom, prom_max, "log", lo_fit, win, cut)
    r2l, r2g = lin["r2"], log["r2"]
    if np.isnan(r2l) and np.isnan(r2g):
        better = ""
    elif np.isnan(r2g) or (not np.isnan(r2l) and r2l >= r2g):
        better = "lin"
    else:
        better = "log"
    return lin, log, better


def _fc_cols(fam, lin, log, better):
    """Column dict for one fold-change family (prefix `fam`, e.g. 'fclr'/'fclr_h').
    Linear-rank fit uses bare names; log-rank fit carries a `_logr` suffix."""
    return {
        f"peak1{fam}": lin["fc1"], f"peak2{fam}": lin["fc2"],
        f"{fam}_r2": lin["r2"], f"{fam}_npts": lin["npts"],
        f"{fam}_sec_lo": lin["sec_lo"], f"{fam}_sec_hi": lin["sec_hi"],
        f"peak1{fam}_logr": log["fc1"], f"peak2{fam}_logr": log["fc2"],
        f"{fam}_logr_r2": log["r2"], f"{fam}_logr_npts": log["npts"],
        f"{fam}_logr_sec_lo": log["sec_lo"], f"{fam}_logr_sec_hi": log["sec_hi"],
        f"{fam}_better": better,
    }


# ── bead loading ────────────────────────────────────────────────────────────
HOMOPOLYMER_MAX_MISMATCH = 2

def near_homopolymer_mask(bcs):
    """Boolean Series: True where a barcode is a (near-)homopolymer of any base —
    its most-frequent of A/C/G/T covers ≥ len − HOMOPOLYMER_MAX_MISMATCH positions."""
    s = bcs.astype(str)
    max_base = pd.concat([s.str.count(b) for b in "ACGT"], axis=1).max(axis=1)
    return max_base >= (s.str.len() - HOMOPOLYMER_MAX_MISMATCH)


def load_cb_bead_table(whitelist, puck=None):
    """Read the pair whitelist → DataFrame[CB, bead_bc, u, x_um, y_um].

    If the whitelist already carries x_um/y_um (cell-barcode_coords.csv with puck
    merged in pair), coords are taken from it and `puck` is not needed. Otherwise
    (legacy coord-free df_whitelist.txt) the puck CSV is merged in. Delimiter is
    auto-detected. Beads with no coords (off-puck) and (near-)homopolymer beads
    are dropped, so neither reaches DBSCAN or KDE."""
    hdr = pd.read_csv(whitelist, sep=None, engine="python", nrows=0).columns
    if "x_um" in hdr and "y_um" in hdr:
        wl = pd.read_csv(whitelist, sep=None, engine="python",
                         usecols=["cell_bc_10x", "bead_bc", "nUMI", "x_um", "y_um"])
        wl = wl.dropna(subset=["x_um", "y_um"])       # drop off-puck beads (blank coords)
        wl = (wl.groupby(["cell_bc_10x", "bead_bc"], as_index=False)
                .agg(u=("nUMI", "sum"), x_um=("x_um", "first"), y_um=("y_um", "first"))
                .rename(columns={"cell_bc_10x": "CB"}))
        return wl[~near_homopolymer_mask(wl["bead_bc"])]
    if puck is None:
        raise ValueError(f"{whitelist} has no x_um/y_um columns and no puck was given")
    wl = pd.read_csv(whitelist, sep=None, engine="python",
                     usecols=["cell_bc_10x", "bead_bc", "nUMI"])
    wl = wl.groupby(["cell_bc_10x", "bead_bc"], as_index=False)["nUMI"].sum().rename(
        columns={"cell_bc_10x": "CB", "nUMI": "u"})
    pk = pd.read_csv(puck, header=None, names=["bead_bc", "x_um", "y_um"])
    pk = pk[~near_homopolymer_mask(pk["bead_bc"])]
    return wl.merge(pk, on="bead_bc")


# ── per-cell ────────────────────────────────────────────────────────────────
def one_cell(cb, x, y, u):
    out = dict(cb=cb, n_beads=int(len(x)), tot_umi=float(u.sum()),
               eps_used=float(PARAMS["eps"][0]))

    # --- DBSCAN: raw + norm weightings, one column set per mode ---
    rest_tag = f"c{N_DBSCAN_CLUSTERS + 1}plus"
    db = dbscan_modes(x, y, u)
    for mode in MODES:
        r = db[mode]; clusters = r["clusters"]; crest = r["crest"]; c0 = r["c0"]
        out[f"{mode}_min_samples_used"] = r["min_samples_used"]
        out[f"{mode}_dbscan_clusters"]  = int(r["n_clusters"])
        for k in range(N_DBSCAN_CLUSTERS):
            c = clusters[k] if k < len(clusters) else None
            out[f"{mode}_c{k+1}_x"]       = c["cx"] if c else np.nan
            out[f"{mode}_c{k+1}_y"]       = c["cy"] if c else np.nan
            out[f"{mode}_c{k+1}_n_beads"] = c["n_beads"] if c else 0
            out[f"{mode}_c{k+1}_umi"]     = c["umi"] if c else 0.0
        out[f"{mode}_{rest_tag}_n_beads"]   = int(crest["n_beads"])
        out[f"{mode}_{rest_tag}_umi"]       = float(crest["umi"])
        out[f"{mode}_{rest_tag}_n_clusters"] = int(crest["n_clusters"])
        out[f"{mode}_c0_n_beads"] = int(c0["n_beads"])
        out[f"{mode}_c0_umi"]     = float(c0["umi"])

    # raw-mode c1/c2 centroids aliased as cx1/cy1/cx2/cy2 (cumumi uses cx1/cy1)
    rawc = db["raw"]["clusters"]
    cx1 = out["raw_c1_x"]; cy1 = out["raw_c1_y"]
    out["cx1"] = cx1; out["cy1"] = cy1
    out["cx2"] = out["raw_c2_x"]; out["cy2"] = out["raw_c2_y"]

    # --- KDE / topo peaks (weighting-independent) ---
    nprom = PARAMS["nprom"]
    have_kde = len(x) >= 2 and u.sum() > 0
    kp = kde_peaks(x, y, u) if have_kde else None
    if kp is None:
        proms = np.array([]); heights = np.array([]); px = np.array([]); py = np.array([])
        out.update(n_topo_peaks=0, bg_median=np.nan, bg_p05=np.nan, bg_p95=np.nan)
    else:
        proms = kp["proms"]; heights = kp["heights"]; px = kp["px"]; py = kp["py"]
        out["n_topo_peaks"] = int(len(proms))
        out["bg_median"] = kp["bg_median"]; out["bg_p05"] = kp["bg_p05"]; out["bg_p95"] = kp["bg_p95"]
    out.update(_peak_cols(heights, proms, kp["bg_median"] if kp else 1e-12))
    out["prom1_over_prom2"] = (float(proms[0] / proms[1])
                               if len(proms) >= 2 and proms[1] > 0 else np.nan)
    out["top_prom"] = float(proms[0]) if len(proms) else np.nan
    out["top_h"]    = float(heights[0]) if len(heights) else np.nan
    out["top_px"]   = float(px[0]) if len(px) else np.nan
    out["top_py"]   = float(py[0]) if len(py) else np.nan
    # each mode's c1 centroid → nearest KDE topo peak
    for mode in MODES:
        c1x = out[f"{mode}_c1_x"]; c1y = out[f"{mode}_c1_y"]
        if len(proms) and np.isfinite(c1x):
            d2 = (px - c1x)**2 + (py - c1y)**2
            k = int(np.argmin(d2))
            out[f"{mode}_c1_prom_db"]      = float(proms[k])
            out[f"{mode}_c1_h_db"]         = float(heights[k])
            out[f"{mode}_c1_peak_dist_db"] = float(np.sqrt(d2[k]))
        else:
            out[f"{mode}_c1_prom_db"] = np.nan
            out[f"{mode}_c1_h_db"] = np.nan
            out[f"{mode}_c1_peak_dist_db"] = np.nan

    # full prominence + height rankings
    p = proms[np.isfinite(proms) & (proms > 0)] if len(proms) else np.array([])
    h = heights[:len(p)] if len(p) else np.array([])
    for i in range(nprom):
        out[f"prom{i+1:04d}"] = float(p[i]) if i < len(p) else np.nan
        out[f"height{i+1:04d}"] = float(h[i]) if i < len(h) else np.nan

    # --- leftmost-line fold changes ---
    out.update(_fc_cols("fclr", dict(_FC_NAN), dict(_FC_NAN), ""))
    out.update(_fc_cols("fclr_h", dict(_FC_NAN), dict(_FC_NAN), ""))
    if len(p) >= LO_FIT + WIN_FIT - 1:
        # fclr: leftmost-line fit on prominence, prom > 0.1 excluded; lin + log rank
        lin, log, better = leftmost_peak_fc_both(p, p, FIT_PROM_MAX_LR)
        out.update(_fc_cols("fclr", lin, log, better))
        # fclr_h: same fit on raw height (re-ranked by height), prom > 0.1 still excluded
        hok = np.isfinite(h) & (h > 0)
        hh = h[hok]; ph = p[hok]
        order = np.argsort(hh)[::-1]
        hh = hh[order]; ph = ph[order]
        lin_h, log_h, better_h = leftmost_peak_fc_both(hh, ph, FIT_PROM_MAX_LR)
        out.update(_fc_cols("fclr_h", lin_h, log_h, better_h))

    # --- cumulative UMI vs radius from raw c1 centroid ---
    radii = PARAMS["radii"]
    if np.isfinite(cx1):
        rr = np.sqrt((x - cx1)**2 + (y - cy1)**2)
        for R in radii:
            out[f"cumumi_r{int(R):03d}"] = float(u[rr <= R].sum())
    else:
        for R in radii:
            out[f"cumumi_r{int(R):03d}"] = np.nan
    cumumi_vals = np.array([out[f"cumumi_r{int(R):03d}"] for R in radii if R <= 150], dtype=float)
    last = float(out.get("cumumi_r150", np.nan))
    out["cumumi_auc"] = float(np.nanmean(cumumi_vals / last)) if np.isfinite(last) and last > 0 else np.nan
    return out


def _safe(cb, xyu):
    if xyu is None or len(xyu[0]) < 1:
        return {"cb": cb, "n_beads": 0}
    try:
        return one_cell(cb, xyu[0], xyu[1], xyu[2])
    except Exception as e:
        return {"cb": cb, "n_beads": -1, "err": str(e)[:60]}


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("--sample", required=True)
    ap.add_argument("--whitelist", required=True,
                    help="pair output cell-barcode_coords.csv (legacy df_whitelist.txt ok)")
    ap.add_argument("--puck", default=None,
                    help="puck CSV (bead_bc, x_um, y_um); only needed for a legacy "
                         "coord-free whitelist")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--eps", nargs="+", type=float, default=DEFAULT_EPS,
                    help="DBSCAN eps (µm); only the first value is used")
    ap.add_argument("--min-samples", nargs="+", type=int, default=DEFAULT_MIN_SAMPLES,
                    help="raw-mode DBSCAN min_samples escalation (descending)")
    ap.add_argument("--norm-to", type=float, default=DEFAULT_NORM_TO,
                    help="norm-mode DBSCAN weights: w = u/∑u * norm_to")
    ap.add_argument("--norm-ms-scan", nargs=2, type=int,
                    default=[DEFAULT_NORM_MS_HI, DEFAULT_NORM_MS_LO], metavar=("HI", "LO"),
                    help="norm-mode integer min_samples scan (high→low); highest tier "
                         "yielding ≥2 clusters wins")
    ap.add_argument("--sigma", type=float, default=DEFAULT_SIGMA, help="KDE σ (µm)")
    ap.add_argument("--grid", type=int, default=DEFAULT_GRID, help="KDE grid (G×G)")
    ap.add_argument("--nprom", type=int, default=DEFAULT_NPROM,
                    help="number of top KDE prominence peaks to record")
    ap.add_argument("--radii", nargs="+", type=float, default=DEFAULT_RADII,
                    help="cumulative-UMI radii (µm) from the raw c1 centroid")
    ap.add_argument("--ncores", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", str(DEFAULT_NCORES))))
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    PARAMS.update(eps=list(args.eps), min_samples=list(args.min_samples),
                  norm_to=args.norm_to,
                  norm_ms_hi=args.norm_ms_scan[0], norm_ms_lo=args.norm_ms_scan[1],
                  sigma=args.sigma, grid=args.grid,
                  nprom=args.nprom, radii=list(args.radii), ncores=args.ncores)

    m = load_cb_bead_table(args.whitelist, args.puck)
    arr = {cb: (g["x_um"].to_numpy(float), g["y_um"].to_numpy(float),
                g["u"].to_numpy(float))
           for cb, g in m.groupby("CB")}
    cells = list(arr.keys())
    print(f"[{args.sample}] {len(cells)} cells, ncores={args.ncores}, "
          f"eps={PARAMS['eps']} min_samples={PARAMS['min_samples']} "
          f"norm_to={PARAMS['norm_to']} nprom={PARAMS['nprom']}")
    t1 = time.time()
    res = Parallel(n_jobs=args.ncores, verbose=5)(
        delayed(_safe)(cb, arr.get(cb)) for cb in cells)
    df = pd.DataFrame(res).set_index("cb")
    df.insert(0, "sample", args.sample)
    out = os.path.join(args.out_dir, "cell_coords.csv")
    df.to_csv(out)
    ok = int(df["raw_c1_umi"].notna().sum()) if "raw_c1_umi" in df else 0
    print(f"[{args.sample}] done {time.time()-t1:.0f}s  wrote {out}  "
          f"({len(df)} cells, {ok} with raw c1, {df.shape[1]} cols)")


if __name__ == "__main__":
    main()
