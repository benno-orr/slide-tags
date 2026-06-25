"""Lasso-select cells from the peak2 scatter, then render bead plots for a random subset.

Needs a graphical display (X11 forwarding or local display).

Usage:
  python lasso_cells.py \\
    --profiles-tsv  STAMP_SAMPLE_cell_profiles.tsv \\
    --whitelist     df_whitelist.txt \\
    --puck          puck.csv \\
    --sample        SAMPLE \\
    --out-dir       /path/out \\
    [--n 16] [--seed 42] [--cutoff p0005]

Draw a lasso on the scatter, then press Enter to render.
Saves:
  STAMP_SAMPLE_lasso_cbs.txt      — one CB per line (full selection)
  STAMP_SAMPLE_lasso_beads.png    — bead scatter grid for N random cells
"""
import argparse, os, sys, datetime as dt
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from map_cells import load_cb_bead_table

CLS_MAP = {
    "singlet":     "steelblue",
    "doublet":     "tomato",
    "unmapped_c1": "grey",
    "dbscan_fail": "black",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles-tsv", required=True)
    ap.add_argument("--whitelist",    required=True)
    ap.add_argument("--puck",         required=True)
    ap.add_argument("--sample",       required=True)
    ap.add_argument("--out-dir",      required=True)
    ap.add_argument("--n",    type=int, default=16, help="cells to render after lasso")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--cutoff", default="p0005",
                    help="class column suffix used for colouring (class_top_<cutoff>)")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    if args.seed is not None:
        np.random.seed(args.seed)

    d = pd.read_csv(args.profiles_tsv, sep="\t", index_col=0)
    lfc2_col = "peak2_obs-prj_lfc"
    prom_col = "prom002"
    cls_col  = f"class_top_{args.cutoff}"

    if lfc2_col not in d.columns or prom_col not in d.columns:
        sys.exit(f"ERROR: profiles TSV missing {lfc2_col!r} or {prom_col!r}")

    d2 = d[d[lfc2_col].notna() & d[prom_col].notna() & (d[prom_col] > 0)].copy()
    d2["_lx"] = np.log10(d2[prom_col])
    d2["_ly"] = d2[lfc2_col]
    has_cls  = cls_col in d2.columns
    xy = np.c_[d2["_lx"].to_numpy(), d2["_ly"].to_numpy()]

    # ── scatter window ────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 7))
    if has_cls:
        for cls, col in CLS_MAP.items():
            m = d2[cls_col] == cls
            if m.any():
                ax.scatter(d2.loc[m, "_lx"], d2.loc[m, "_ly"],
                           s=4, alpha=0.35, color=col, linewidth=0,
                           label=f"{cls} ({m.sum():,})")
        ax.legend(fontsize=8, markerscale=3)
    else:
        ax.scatter(d2["_lx"], d2["_ly"], s=4, alpha=0.35,
                   color="steelblue", linewidth=0)

    ax.axhline(0, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("log10(prom002)  [2nd KDE peak prominence]")
    ax.set_ylabel("peak2_lfc  [obs − projected]")
    ax.set_title(f"{args.sample}  —  draw lasso, then press Enter", fontweight="bold")
    ax.grid(alpha=0.3)
    fig.tight_layout()

    selected_idx = []

    def onselect(verts):
        mask = Path(verts).contains_points(xy)
        selected_idx.clear()
        selected_idx.extend(np.where(mask)[0].tolist())
        n_sel = len(selected_idx)
        n_rnd = min(args.n, n_sel)
        ax.set_title(f"{args.sample}  —  {n_sel:,} selected. "
                     f"Press Enter to render {n_rnd}.", fontweight="bold")
        fig.canvas.draw_idle()

    lasso = LassoSelector(ax, onselect, useblit=True)  # noqa: F841 (keep ref)

    def on_key(event):
        if event.key == "enter" and selected_idx:
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()

    if not selected_idx:
        print("No cells selected — exiting.")
        return

    sel_cbs  = d2.index[selected_idx].tolist()
    n_render = min(args.n, len(sel_cbs))
    chosen   = np.random.choice(sel_cbs, size=n_render, replace=False)
    stamp    = dt.datetime.now().strftime("%Y-%m-%d_%H-%M")

    sel_path = os.path.join(args.out_dir, f"{stamp}_{args.sample}_lasso_cbs.txt")
    pd.Series(sel_cbs).to_csv(sel_path, index=False, header=False)
    print(f"Selected {len(sel_cbs):,} cells  →  {sel_path}")

    # ── bead scatter grid ─────────────────────────────────────────────────────
    print("Loading bead table …")
    m   = load_cb_bead_table(args.whitelist, args.puck)
    arr = {cb: (g["x_um"].to_numpy(float), g["y_um"].to_numpy(float),
                g["u"].to_numpy(float))
           for cb, g in m.groupby("CB")}

    ncol = min(4, n_render)
    nrow = int(np.ceil(n_render / ncol))
    figR, axes = plt.subplots(nrow, ncol, figsize=(ncol * 3.5, nrow * 3.5))
    figR.suptitle(f"{args.sample}  —  {n_render} random lasso cells (bead scatter)",
                  fontsize=11, fontweight="bold")
    axes_flat = np.atleast_1d(axes).flat

    for i, cb in enumerate(chosen):
        axb = axes_flat[i]
        xyu = arr.get(cb)
        if xyu is None:
            axb.set_title(f"{cb[:16]}\n(no beads)", fontsize=7)
            axb.axis("off")
            continue
        x, y, u = xyu
        sc = axb.scatter(x, y, c=u, s=5, cmap="viridis", linewidth=0)
        axb.set_aspect("equal", adjustable="box")
        axb.set_xticks([]); axb.set_yticks([])
        lbl = d2.loc[cb, cls_col] if (has_cls and cb in d2.index) else ""
        axb.set_title(f"{cb[:16]}  {lbl}", fontsize=7)
        figR.colorbar(sc, ax=axb, fraction=0.04, pad=0.02)

    for i in range(n_render, nrow * ncol):
        axes_flat[i].set_visible(False)

    figR.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = os.path.join(args.out_dir, f"{stamp}_{args.sample}_lasso_beads.png")
    figR.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Wrote  {out_path}")


if __name__ == "__main__":
    main()
