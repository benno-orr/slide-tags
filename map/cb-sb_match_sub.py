#!/usr/bin/env python3
"""
cb-sb_match_sub.py — tags pipeline launcher.

For each row in a file_paths.csv where run == TRUE, submit a cb-sb_match sbatch job.
Output and log locations are caller-supplied so the pipeline is reusable across analyses.

Usage:
    python cb-sb_match_sub.py [--odir /path/to/outputs] \
                              [--csv config/file_paths.csv] \
                              [--log-dir /path/to/logs] [--dry-run]

Output base is per-row from the CSV 'odir' column, falling back to --odir for blank
rows; one of the two must supply it. Each invocation creates a timestamped run dir
<base>/<YYMMDD_HHMMSS>/, and matcher output lands in <run-dir>/<spl>/df_whitelist.txt.
SLURM logs (<log-dir>/%j_cb-sb_match.{out,err}) override the sbatch wrapper defaults.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
from datetime import datetime

PIPE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(PIPE_DIR, "config", "file_paths.csv")
SBATCH_SCRIPT = os.path.join(PIPE_DIR, "cb-sb_match_sbatch.sh")


def pipeline_version():
    """Git commit + dirty flag of the pipeline checkout, for output provenance."""
    def git(*a):
        return subprocess.check_output(
            ["git", "-C", PIPE_DIR, *a], text=True, stderr=subprocess.DEVNULL
        ).strip()
    try:
        commit = git("rev-parse", "HEAD")
        dirty = bool(git("status", "--porcelain"))
    except Exception:
        return {"commit": "unknown", "dirty": None}
    return {"commit": commit, "dirty": dirty}


def running_job_names():
    """Set of job names currently in this user's SLURM queue."""
    try:
        out = subprocess.check_output(
            ["squeue", "-u", os.environ.get("USER", "beo703"),
             "--format=%j", "--noheader"],
            text=True,
        )
        return set(out.split())
    except Exception:
        return set()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-odir", "--odir", "--odir-base", dest="odir_base", default=None,
                   help="fallback output base dir, used for rows whose CSV 'odir' "
                        "column is blank. Output lands in <odir>/<YYMMDD_HHMMSS>/<spl>/.")
    p.add_argument("--csv", default=DEFAULT_CSV,
                   help=f"sample table (default: {DEFAULT_CSV})")
    p.add_argument("--log-dir", default=None,
                   help="dir for SLURM logs (default: <run-dir>/logs)")
    p.add_argument("--dry-run", action="store_true",
                   help="print sbatch commands without submitting")
    args = p.parse_args()

    # One timestamp per invocation. Each output base (CSV 'odir' column, else --odir)
    # gets its own timestamped run dir <base>/<YYMMDD_HHMMSS>/, so re-runs never clobber.
    stamp = datetime.now().strftime("%y%m%d_%H%M%S")
    cli_base = os.path.abspath(args.odir_base) if args.odir_base else None
    ver = pipeline_version()
    print(f"pipeline commit: {ver['commit'][:12]}"
          f"{' (dirty)' if ver['dirty'] else ''}")

    run_dirs = {}  # output base -> its timestamped run dir, provisioned on first use

    def run_dir_for(base):
        """Resolve and (once) provision the timestamped run dir for an output base."""
        base = os.path.abspath(base)
        if base not in run_dirs:
            rdir = os.path.join(base, stamp)
            print(f"run dir: {rdir}")
            if not args.dry_run:
                os.makedirs(rdir, exist_ok=True)
                # Provenance: pin each run dir to the pipeline version + inputs.
                with open(os.path.join(rdir, "pipeline_version.json"), "w") as fh:
                    json.dump({
                        "timestamp": stamp,
                        "pipeline_dir": PIPE_DIR,
                        "git_commit": ver["commit"],
                        "git_dirty": ver["dirty"],
                        "csv": os.path.abspath(args.csv),
                    }, fh, indent=2)
                shutil.copy(args.csv, os.path.join(rdir, "file_paths.used.csv"))
            run_dirs[base] = rdir
        return run_dirs[base]

    queued = running_job_names()

    with open(args.csv) as fh:
        reader = csv.DictReader(fh, skipinitialspace=True)
        for row in reader:
            if row.get("run", "").strip() != "TRUE":
                continue
            spl      = row["spl"].strip()
            r1       = row["r1_fq"].strip().strip('"')
            r2       = row["r2_fq"].strip().strip('"')
            wl       = row["wl"].strip().strip('"')
            cb_fmats = row.get("cb_fmats", "").strip().strip('"')
            puck     = row.get("puck_csv", "").strip().strip('"')
            base     = row.get("odir", "").strip().strip('"') or cli_base
            if not wl:
                print(f"  SKIP {spl}: no whitelist")
                continue
            if not base:
                print(f"  SKIP {spl}: no odir (set the CSV 'odir' column or pass --odir)")
                continue
            if f"cbsb_{spl}" in queued:
                print(f"  SKIP {spl}: already queued")
                continue
            rdir = run_dir_for(base)
            log_dir = os.path.abspath(args.log_dir) if args.log_dir \
                else os.path.join(rdir, "logs")
            if not args.dry_run:
                os.makedirs(log_dir, exist_ok=True)
            odir = os.path.join(rdir, spl)
            cmd = ["sbatch", f"--job-name=cbsb_{spl}",
                   f"--export=ALL,CBSB_PIPE_DIR={PIPE_DIR}",
                   "-o", os.path.join(log_dir, "%j_cb-sb_match.out"),
                   "-e", os.path.join(log_dir, "%j_cb-sb_match.err"),
                   SBATCH_SCRIPT, r1, r2, wl, odir]
            if cb_fmats:
                cmd += ["--cb-fmats", cb_fmats]
            if puck:
                cmd += ["--puck", puck]
            print(" ".join(cmd))
            if not args.dry_run:
                subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
