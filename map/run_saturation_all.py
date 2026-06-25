#!/usr/bin/env python3
"""
run_saturation_all.py — tags pipeline: SB library-complexity over every matched sample.

Runs saturation.py (analytic downsampling + Michaelis-Menten projection) for every
<restored-dir>/<sample>/df_whitelist.txt produced by the cb-sb matcher.

Usage:
    python run_saturation_all.py --restored-dir /path/to/outputs/restored \
                                 [--plots-dir DIR] [--plt-id pltBOxx] \
                                 [--samples xBO203 ...] [--force]

Per sample writes saturation_curve.csv + saturation_stats.csv into its restored dir, and
a saturation curve plot into --plots-dir (default: each sample's own restored dir).

NOTE: --plt-id defaults to the placeholder "pltBOxx". Plot filenames carry the sample name;
registering real pltBO ids in tables/plots.csv is a separate step.
"""

import argparse
import glob
import os
import subprocess
import sys

PIPE_DIR = os.path.dirname(os.path.abspath(__file__))
SATURATION = os.path.join(PIPE_DIR, "saturation.py")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--restored-dir", required=True,
                   help="dir holding per-sample <sample>/cell-barcode_coords.csv")
    p.add_argument("--plots-dir", default=None,
                   help="dir for saturation plots (default: each sample's restored dir)")
    p.add_argument("--plt-id", default="pltBOxx",
                   help="pltBO id passed to saturation.py (default: placeholder pltBOxx)")
    p.add_argument("--samples", nargs="*", default=None,
                   help="restrict to these sample names (default: all found)")
    p.add_argument("--force", action="store_true",
                   help="recompute even if saturation_stats.csv already exists")
    args = p.parse_args()

    restored = os.path.abspath(args.restored_dir)
    whitelists = sorted(glob.glob(os.path.join(restored, "*", "cell-barcode_coords.csv"))
                        + glob.glob(os.path.join(restored, "*", "df_whitelist.txt")))
    if not whitelists:
        print(f"No cell-barcode_coords.csv found under {restored}", file=sys.stderr)
        sys.exit(1)

    n_run = n_skip = n_fail = 0
    for df_fp in whitelists:
        odir = os.path.dirname(df_fp)
        spl = os.path.basename(odir)
        if args.samples and spl not in args.samples:
            continue
        if not args.force and os.path.exists(os.path.join(odir, "saturation_stats.csv")):
            print(f"SKIP {spl}: saturation_stats.csv exists (use --force)")
            n_skip += 1
            continue
        plots_dir = os.path.abspath(args.plots_dir) if args.plots_dir else odir
        print(f"=== saturation: {spl} ===", flush=True)
        try:
            # saturation.py positional args: df_fp odir spl plt_id plots_dir
            subprocess.run(
                [sys.executable, SATURATION, df_fp, odir, spl, args.plt_id, plots_dir],
                check=True,
            )
            n_run += 1
        except subprocess.CalledProcessError as e:
            print(f"FAILED {spl}: {e}", file=sys.stderr)
            n_fail += 1

    print(f"\nDone. ran={n_run} skipped={n_skip} failed={n_fail}")


if __name__ == "__main__":
    main()
