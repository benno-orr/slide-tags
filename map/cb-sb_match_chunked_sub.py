#!/usr/bin/env python3
"""
cb-sb_match_chunked_sub.py — chunked launcher for CB-SB matching.

Splits each sample's FASTQs into chunks, submits one sbatch job per chunk,
then submits a merge job that depends on all chunk jobs completing successfully.

Usage:
    python cb-sb_match_chunked_sub.py [--odir <base>] [--csv config/file_paths.csv]
                                       [--chunk-size 20000000] [--scratch <dir>]
                                       [--log-dir <dir>] [--dry-run]

Run from a compute node — the split step decompresses both FASTQs with pigz
(~5 min per sample) and should not be run on a login node.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
from datetime import datetime

PIPE_DIR      = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV   = os.path.join(PIPE_DIR, "config", "file_paths.csv")
CHUNK_SBATCH  = os.path.join(PIPE_DIR, "cb-sb_match_chunk_sbatch.sh")
MERGE_SBATCH  = os.path.join(PIPE_DIR, "cb-sb_match_merge_sbatch.sh")
DEFAULT_SCRATCH = "/n/scratch/users/b/beo703"


def pipeline_version():
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


def split_fastqs(r1_fp, r2_fp, r1_dir, r2_dir, lines_per_chunk):
    """Decompress and split R1 and R2 in parallel using pigz. Returns (r1_chunks, r2_chunks)."""
    os.makedirs(r1_dir, exist_ok=True)
    os.makedirs(r2_dir, exist_ok=True)

    r1_decomp = subprocess.Popen(["pigz", "-dc", r1_fp], stdout=subprocess.PIPE)
    r1_split  = subprocess.Popen(
        ["split", "-l", str(lines_per_chunk), "-d", "--suffix-length=3",
         "-", os.path.join(r1_dir, "r1_")],
        stdin=r1_decomp.stdout,
    )
    r1_decomp.stdout.close()

    r2_decomp = subprocess.Popen(["pigz", "-dc", r2_fp], stdout=subprocess.PIPE)
    r2_split  = subprocess.Popen(
        ["split", "-l", str(lines_per_chunk), "-d", "--suffix-length=3",
         "-", os.path.join(r2_dir, "r2_")],
        stdin=r2_decomp.stdout,
    )
    r2_decomp.stdout.close()

    r1_decomp.wait(); r1_split.wait()
    r2_decomp.wait(); r2_split.wait()

    r1_chunks = sorted(os.path.join(r1_dir, f) for f in os.listdir(r1_dir) if f.startswith("r1_"))
    r2_chunks = sorted(os.path.join(r2_dir, f) for f in os.listdir(r2_dir) if f.startswith("r2_"))
    return r1_chunks, r2_chunks


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-odir", "--odir", "--odir-base", dest="odir_base", default=None,
                   help="fallback output base dir for rows with blank CSV 'odir' column")
    p.add_argument("--csv", default=DEFAULT_CSV,
                   help=f"sample table (default: {DEFAULT_CSV})")
    p.add_argument("--chunk-size", type=int, default=20_000_000,
                   help="reads per chunk (default: 20,000,000)")
    p.add_argument("--scratch", default=DEFAULT_SCRATCH,
                   help=f"scratch dir for FASTQ chunks (default: {DEFAULT_SCRATCH})")
    p.add_argument("--log-dir", default=None,
                   help="dir for SLURM logs (default: <run-dir>/logs)")
    p.add_argument("--dry-run", action="store_true",
                   help="print actions without splitting or submitting")
    args = p.parse_args()

    stamp    = datetime.now().strftime("%y%m%d_%H%M%S")
    cli_base = os.path.abspath(args.odir_base) if args.odir_base else None
    ver      = pipeline_version()
    print(f"pipeline commit: {ver['commit'][:12]}{' (dirty)' if ver['dirty'] else ''}")

    run_dirs = {}

    def run_dir_for(base):
        base = os.path.abspath(base)
        if base not in run_dirs:
            rdir = os.path.join(base, stamp)
            print(f"run dir: {rdir}")
            if not args.dry_run:
                os.makedirs(rdir, exist_ok=True)
                with open(os.path.join(rdir, "pipeline_version.json"), "w") as fh:
                    json.dump({
                        "timestamp": stamp, "pipeline_dir": PIPE_DIR,
                        "git_commit": ver["commit"], "git_dirty": ver["dirty"],
                        "csv": os.path.abspath(args.csv),
                    }, fh, indent=2)
                shutil.copy(args.csv, os.path.join(rdir, "file_paths.used.csv"))
            run_dirs[base] = rdir
        return run_dirs[base]

    lines_per_chunk = args.chunk_size * 4

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
                print(f"  SKIP {spl}: no whitelist"); continue
            if not base:
                print(f"  SKIP {spl}: no odir"); continue

            rdir    = run_dir_for(base)
            log_dir = os.path.abspath(args.log_dir) if args.log_dir \
                else os.path.join(rdir, "logs")
            odir    = os.path.join(rdir, spl)
            chunk_scratch = os.path.join(args.scratch, f"cbsb_chunks_{stamp}_{spl}")

            if not args.dry_run:
                os.makedirs(log_dir, exist_ok=True)
                os.makedirs(odir, exist_ok=True)

            # Split FASTQs
            print(f"{spl}: splitting into {args.chunk_size/1e6:.0f}M-read chunks "
                  f"(scratch: {chunk_scratch})...", flush=True)
            if args.dry_run:
                print(f"  [dry-run] would split {r1} and {r2}")
                continue

            r1_chunks, r2_chunks = split_fastqs(
                r1, r2,
                os.path.join(chunk_scratch, "r1"),
                os.path.join(chunk_scratch, "r2"),
                lines_per_chunk,
            )
            n_chunks = len(r1_chunks)
            if len(r1_chunks) != len(r2_chunks):
                raise RuntimeError(
                    f"R1 ({len(r1_chunks)}) and R2 ({len(r2_chunks)}) chunk counts differ"
                )
            print(f"  {n_chunks} chunks", flush=True)

            # Submit chunk jobs
            chunk_job_ids = []
            for i, (r1c, r2c) in enumerate(zip(r1_chunks, r2_chunks)):
                chunk_odir = os.path.join(chunk_scratch, f"chunk_{i:03d}")
                cmd = [
                    "sbatch",
                    f"--job-name=cbsb_{spl}_c{i:03d}",
                    f"--export=ALL,CBSB_PIPE_DIR={PIPE_DIR}",
                    "-o", os.path.join(log_dir, f"%j_chunk_{i:03d}.out"),
                    "-e", os.path.join(log_dir, f"%j_chunk_{i:03d}.err"),
                    CHUNK_SBATCH, r1c, r2c, wl, chunk_odir,
                ]
                if cb_fmats:
                    cmd += ["--cb-fmats", cb_fmats]
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                job_id = result.stdout.strip().split()[-1]
                chunk_job_ids.append(job_id)
                print(f"  chunk {i:03d}: job {job_id}")

            # Write manifest for merge job
            manifest_path = os.path.join(chunk_scratch, "partial_files.txt")
            with open(manifest_path, "w") as fh:
                for i in range(n_chunks):
                    fh.write(os.path.join(chunk_scratch, f"chunk_{i:03d}",
                                          "partial_triplets.tsv") + "\n")

            # Submit merge job depending on all chunk jobs
            dep = "afterok:" + ":".join(chunk_job_ids)
            merge_cmd = [
                "sbatch",
                f"--job-name=cbsb_{spl}_merge",
                f"--dependency={dep}",
                f"--export=ALL,CBSB_PIPE_DIR={PIPE_DIR}",
                "-o", os.path.join(log_dir, "%j_merge.out"),
                "-e", os.path.join(log_dir, "%j_merge.err"),
                MERGE_SBATCH, manifest_path, odir, spl, chunk_scratch, puck,
            ]
            result = subprocess.run(merge_cmd, check=True, capture_output=True, text=True)
            merge_job_id = result.stdout.strip().split()[-1]
            print(f"  merge job: {merge_job_id} (depends on {n_chunks} chunks)")


if __name__ == "__main__":
    main()
