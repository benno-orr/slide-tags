#!/usr/bin/env python3
"""
Sequencing saturation curve for CB-SB matching output.

Downsamples analytically (no simulation) using Poisson approximation.
Fits Michaelis-Menten curve and projects library complexity.

Usage:
    python saturation.py <df_whitelist_fp> <odir> [spl] [plt_id] [plots_dir]

plots_dir defaults to <odir> (per-sample), so the pipeline writes no global plots dir.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

chunk_id = "saturation-curve"
version  = 1

# ── args ──────────────────────────────────────────────────────────────────────

df_fp   = sys.argv[1]
odir    = sys.argv[2]
spl     = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(odir)
plt_id  = sys.argv[4] if len(sys.argv) > 4 else "pltBOxx"
plt_dir = sys.argv[5] if len(sys.argv) > 5 else odir
os.makedirs(odir, exist_ok=True)
os.makedirs(plt_dir, exist_ok=True)

# ── load ──────────────────────────────────────────────────────────────────────

print("Loading cell-barcode_coords...", flush=True)
df = pd.read_csv(df_fp, sep=None, engine="python", usecols=["nUMI", "nReads"])
print(f"  {len(df):,} CB-SB pairs", flush=True)

r = df["nReads"].values.astype(np.float64)
u = df["nUMI"].values.astype(np.float64)
N_total = r.sum()
print(f"  Total reads: {N_total:,.0f}", flush=True)

# ── analytical downsampling ───────────────────────────────────────────────────
# fractions 0..1 for fitting; beyond 1 for extrapolation

fracs = np.concatenate([
    np.linspace(0.02, 1.0, 50),
    np.array([1.5, 2.0, 3.0, 5.0, 10.0]),
])

# Pre-UMI collapse: E[unique CB-SB pairs] = sum_i P(>=1 read from pair i at fraction f)
#   P(>=1) = 1 - (1-f)^r_i  [binomial]
# Post-UMI collapse: E[unique UMIs] = sum_i u_i * (1 - exp(-f * r_i/u_i))
#   (Poisson approx: reads distributed across UMIs within each pair)

pre_umi  = np.zeros(len(fracs))
post_umi = np.zeros(len(fracs))

print("Computing downsampling curves...", flush=True)
for i, f in enumerate(fracs):
    # use log-space to avoid overflow for large r: (1-f)^r = exp(r*log(1-f))
    log1mf = np.log1p(-f) if f < 1.0 else -np.inf
    pre_umi[i]  = np.sum(1 - np.exp(r * log1mf))
    rate = np.where(u > 0, f * r / u, f * r)
    post_umi[i] = np.nansum(u * (1 - np.exp(-rate)))

reads_axis = fracs * N_total  # absolute read counts

# ── Michaelis-Menten fit (on observed range f<=1) ─────────────────────────────

def michaelis_menten(x, vmax, km):
    return vmax * x / (km + x)

obs_mask = fracs <= 1.0
x_obs = reads_axis[obs_mask]

p0 = [pre_umi[obs_mask].max() * 2, N_total]

popt_pre,  _ = curve_fit(michaelis_menten, x_obs, pre_umi[obs_mask],
                          p0=p0, maxfev=10000)
popt_post, _ = curve_fit(michaelis_menten, x_obs, post_umi[obs_mask],
                          p0=p0, maxfev=10000)

vmax_pre,  km_pre  = popt_pre
vmax_post, km_post = popt_post

# extrapolation axis: 0 → 5× current depth
x_fit = np.linspace(0, reads_axis.max(), 500)
y_fit_pre  = michaelis_menten(x_fit, *popt_pre)
y_fit_post = michaelis_menten(x_fit, *popt_post)

print(f"  Pre-UMI  Vmax = {vmax_pre:,.0f}  Km = {km_pre:,.0f}", flush=True)
print(f"  Post-UMI Vmax = {vmax_post:,.0f}  Km = {km_post:,.0f}", flush=True)

# ── summary stats ─────────────────────────────────────────────────────────────

# reads to reach 75% of Vmax: solve Vmax*x/(Km+x) = 0.75*Vmax → x = 3*Km
reads_at_75 = 3 * km_post
additional_reads = max(0.0, reads_at_75 - N_total)
pairs_at_current = float(michaelis_menten(N_total, *popt_post))

stats = {
    "observed_reads":          int(N_total),
    "vmax_pairs":              int(vmax_post),
    "km":                      int(km_post),
    "pairs_at_current_depth":  int(pairs_at_current),
    "pct_of_vmax_current":     round(100 * pairs_at_current / vmax_post, 1),
    "reads_to_reach_75pct":    int(reads_at_75),
    "additional_reads_for_75pct": int(additional_reads),
}
for k, v in stats.items():
    print(f"  {k}: {v:,}" if isinstance(v, int) else f"  {k}: {v}", flush=True)

# saturation curve table
sat_df = pd.DataFrame({
    "fraction":       fracs,
    "total_reads":    reads_axis,
    "post_umi_pairs": post_umi,
})
sat_df.to_csv(os.path.join(odir, "saturation_curve.csv"), index=False)

pd.DataFrame([stats]).to_csv(os.path.join(odir, "saturation_stats.csv"), index=False)

# ── plot ──────────────────────────────────────────────────────────────────────

# ── chunk: saturation-curve ──────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(4.5, 3.5))

x_m     = reads_axis / 1e6
x_fit_m = x_fit / 1e6
c       = "#EE854A"
grey    = "#888888"

ax.plot(x_m[obs_mask], post_umi[obs_mask] / 1e6,
        "o", ms=3, color=c, label="observed", zorder=3)
ax.plot(x_fit_m, y_fit_post / 1e6, "-", color=c, lw=1.5, label="M-M fit")

# Vmax asymptote
ax.axhline(vmax_post / 1e6, ls="--", lw=0.8, color=c, alpha=0.6)
ax.text(x_fit_m[-1] * 0.98, vmax_post / 1e6 * 1.02,
        f"Vmax = {vmax_post/1e6:.1f}M", ha="right", va="bottom",
        fontsize=8, color=c)

# 75% line
y75 = 0.75 * vmax_post
ax.axhline(y75 / 1e6, ls="--", lw=0.8, color=grey, alpha=0.7)
ax.axvline(reads_at_75 / 1e6, ls="--", lw=0.8, color=grey, alpha=0.7)
ax.text(reads_at_75 / 1e6 * 1.01, ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else 0,
        f"75%: {reads_at_75/1e6:.1f}M reads", ha="left", va="bottom",
        fontsize=7, color=grey)

# current depth
ax.axvline(N_total / 1e6, ls=":", lw=0.8, color=grey)
ax.text(N_total / 1e6 * 1.01, pairs_at_current / 1e6,
        f"now: {N_total/1e6:.1f}M reads\n+{additional_reads/1e6:.1f}M needed",
        ha="left", va="center", fontsize=7, color=grey)

ax.set_xlabel("Total reads (M)", fontsize=11)
ax.set_ylabel("Unique CB-SB pairs (M)", fontsize=11)
ax.tick_params(labelsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(fontsize=8, frameon=False)

fig.tight_layout()

for ext in ("pdf", "png"):
    fig.savefig(
        os.path.join(plt_dir, f"{plt_id}_{chunk_id}_v{version}_{spl}.{ext}"),
        dpi=300, bbox_inches="tight",
    )

print(f"Saved plots to {plt_dir}", flush=True)
# ── end chunk ────────────────────────────────────────────────────────────────
