#!/usr/bin/env python3
"""
cb-sb_match.py

Match 10x cell barcodes (CBs) from R1 FASTQ reads with spatial barcodes (SBs)
from R2 FASTQ reads, filtered against a cell-barcode whitelist.

Usage:
    python cb-sb_match.py <r1_fp> <r2_fp> <wl_fp> <odir> [--workers N] [--batch-size N]

Arguments:
    r1_fp         R1 FASTQ (.fastq.gz) — contains cell barcode (pos 1-16) + UMI (pos 17-28)
    r2_fp         R2 FASTQ (.fastq.gz) — contains spatial barcode (SBa pos 1-8, SBb pos 27-32)
    wl_fp         Cell barcode whitelist (one barcode per line, with or without -1 suffix)
    odir          Output directory

Optional:
    --workers     Number of parallel worker processes (default: all CPUs)
    --batch-size  Reads per batch dispatched to workers (default: 100000)
    --tmp-dir     Directory for sort temp files (default: odir)
"""

import argparse
import csv
import gzip
import io
import json
import os
import subprocess
from collections import defaultdict
from itertools import islice
from multiprocessing import Pool


# ── reference data ────────────────────────────────────────────────────────────

def load_cb_fmats(cb_fmats_fp: str) -> dict:
    """Returns atac_rc_to_gex dict from arc_cb_fmats.csv."""
    import csv as _csv
    mapping = {}
    with open(cb_fmats_fp) as fh:
        for row in _csv.DictReader(fh):
            mapping[row["cb_atac_rc"].strip('"')] = row["cb_gex"].strip('"')
    print(f"CB format map: {len(mapping):,} entries", flush=True)
    return mapping


def load_whitelist(wl_fp: str, cb_fmats_fp: str = None):
    """Returns (valid_cb_set, n_whitelist).

    If cb_fmats_fp is given, treats whitelist barcodes as cb_atac_rc and
    converts them to cb_gex via arc_cb_fmats.csv before building the set.
    Otherwise, strips -1 suffixes and uses barcodes as-is (GEX format).
    """
    mapping = load_cb_fmats(cb_fmats_fp) if cb_fmats_fp else None
    valid_cb_set = set()
    open_fn = gzip.open if wl_fp.endswith(".gz") else open
    skipped = 0
    with open_fn(wl_fp, "rt") as fh:
        for line in fh:
            bc = line.strip()
            if bc.endswith("-1"):
                bc = bc[:-2]
            if mapping is not None:
                gex = mapping.get(bc)
                if gex is None:
                    skipped += 1
                    continue
                bc = gex
            valid_cb_set.add(bc)
    n = len(valid_cb_set)
    if skipped:
        print(f"Whitelist: {n:,} barcodes ({skipped:,} not found in cb_fmats)", flush=True)
    else:
        print(f"Whitelist: {n:,} barcodes", flush=True)
    return valid_cb_set, n


# ── worker ────────────────────────────────────────────────────────────────────

_worker_valid_cb_set = None


def _init_worker(valid_cb_set):
    global _worker_valid_cb_set
    _worker_valid_cb_set = valid_cb_set


def process_batch(batch):
    """
    batch: list of (r1_seq, r2_seq) strings
    Returns list of (cb, sb, umi) tuples for reads passing all filters.
    """
    results = []
    valid_cb_set = _worker_valid_cb_set
    for r1_seq, r2_seq in batch:
        cb = r1_seq[0:16]
        if cb not in valid_cb_set:
            continue
        umi = r1_seq[16:28]
        sb  = r2_seq[0:8] + r2_seq[26:32]  # SBa (8bp) + SBb (6bp) = 14bp
        results.append((cb, sb, umi))
    return results


# ── FASTQ streaming ───────────────────────────────────────────────────────────

def fastq_pair_reader(r1_fp: str, r2_fp: str, batch_size: int):
    """Generator yielding batches of (r1_seq, r2_seq) tuples.

    Uses pigz for parallel FASTQ decompression; faster than single-threaded gzip.open.
    """
    r1_proc = subprocess.Popen(["pigz", "-dc", r1_fp], stdout=subprocess.PIPE)
    r2_proc = subprocess.Popen(["pigz", "-dc", r2_fp], stdout=subprocess.PIPE)
    try:
        r1_fh = io.TextIOWrapper(r1_proc.stdout)
        r2_fh = io.TextIOWrapper(r2_proc.stdout)
        while True:
            batch = []
            for _ in range(batch_size):
                r1_lines = list(islice(r1_fh, 4))
                r2_lines = list(islice(r2_fh, 4))
                if not r1_lines or not r2_lines:
                    break
                batch.append((r1_lines[1].strip(), r2_lines[1].strip()))
            if not batch:
                break
            yield batch
    finally:
        r1_proc.stdout.close()
        r2_proc.stdout.close()
        r1_proc.wait()
        r2_proc.wait()


# ── aggregation ───────────────────────────────────────────────────────────────

def _load_puck(puck_fp):
    """Read headerless puck CSV (bead_bc, x_um, y_um) → {bead_bc: (x_um, y_um)}.
    Returns None when no usable puck is given, so the writer falls back to the
    5-column (coord-free) format."""
    if not puck_fp or not os.path.exists(puck_fp):
        return None
    coords = {}
    with open(puck_fp, newline="") as fh:
        for row in csv.reader(fh):
            if len(row) >= 3:
                coords[row[0]] = (row[1], row[2])
    return coords


def aggregate_and_write(tmp_path: str, odir: str, tmp_dir: str,
                        total_reads: int, n_whitelist: int, puck=None):
    sorted_path = tmp_path + ".sorted"
    print("Sorting triplets...", flush=True)
    subprocess.run(
        [
            "sort",
            "--buffer-size=2G",
            "--parallel=4",
            f"-T{tmp_dir}",
            "-t\t",
            "-k1,1", "-k2,2", "-k3,3",
            "-o", sorted_path,
            tmp_path,
        ],
        check=True,
    )
    os.remove(tmp_path)

    print("Aggregating...", flush=True)
    cb_sb_nreads = defaultdict(int)
    cb_sb_numi   = defaultdict(int)
    sb_nreads    = defaultdict(int)
    matched_cbs  = set()
    prev_triplet = None

    with open(sorted_path) as fh:
        for line in fh:
            cb, sb, umi = line.rstrip("\n").split("\t")
            cb_sb   = (cb, sb)
            triplet = (cb, sb, umi)
            cb_sb_nreads[cb_sb] += 1
            sb_nreads[sb]       += 1
            matched_cbs.add(cb)
            if triplet != prev_triplet:
                cb_sb_numi[cb_sb] += 1
            prev_triplet = triplet

    os.remove(sorted_path)

    total_matched_reads = sum(cb_sb_nreads.values())
    cb_matched  = len(matched_cbs)
    cb_sb_unique = len(cb_sb_nreads)

    _write_df_whitelist(cb_sb_nreads, cb_sb_numi, odir, puck)
    _write_reads_per_sb(sb_nreads, odir)
    _write_summary(
        odir,
        UP_site_matches=total_reads,
        cell_bender_cb=n_whitelist,
        cb_matched=cb_matched,
        cb_matched_reads=total_matched_reads,
        CB_SB_unique_pairings=cb_sb_unique,
    )

    print(
        f"Done.\n"
        f"  Total reads:          {total_reads:,}\n"
        f"  CB-matched reads:     {total_matched_reads:,} "
        f"({100*total_matched_reads/max(total_reads,1):.1f}%)\n"
        f"  Unique CBs matched:   {cb_matched:,} / {n_whitelist:,}\n"
        f"  Unique CB-SB pairs:   {cb_sb_unique:,}",
        flush=True,
    )


# ── output writers ────────────────────────────────────────────────────────────

def _write_df_whitelist(cb_sb_nreads, cb_sb_numi, odir, puck=None):
    """Write cell-barcode_coords.csv. With a puck dict, left-join x_um,y_um by
    bead_bc (blank if the bead is not on the puck); without one, write the
    5-column coord-free format."""
    path = os.path.join(odir, "cell-barcode_coords.csv")
    header = ["CB_SB", "nUMI", "nReads", "cell_bc_10x", "bead_bc"]
    if puck is not None:
        header += ["x_um", "y_um"]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for (cb, sb), nreads in cb_sb_nreads.items():
            cb_sb = cb + sb
            bead_bc = cb_sb[16:30]
            row = [cb_sb, cb_sb_numi[(cb, sb)], nreads, cb_sb[0:16], bead_bc]
            if puck is not None:
                x, y = puck.get(bead_bc, ("", ""))
                row += [x, y]
            writer.writerow(row)


def _write_reads_per_sb(sb_nreads, odir):
    # Transposed matrix format: 1 header row + 1 data row, columns = SBs
    path = os.path.join(odir, "reads_per_SB.csv")
    sorted_sbs = sorted(sb_nreads.keys())
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, quoting=csv.QUOTE_ALL)
        writer.writerow([""] + sorted_sbs)
        writer.writerow(["1"] + [sb_nreads[sb] for sb in sorted_sbs])


def _write_summary(odir, **metrics):
    path = os.path.join(odir, "matcher_summary.txt")
    all_metrics = {"hamming_dist": 0, **metrics}
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(list(all_metrics.keys()))
        writer.writerow(list(all_metrics.values()))


# ── main ──────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("r1_fp")
    p.add_argument("r2_fp")
    p.add_argument("wl_fp")
    p.add_argument("odir")
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--batch-size", type=int, default=100_000)
    p.add_argument("--tmp-dir", default=None)
    p.add_argument("--cb-fmats", default=None,
                   help="arc_cb_fmats.csv path; if given, whitelist barcodes are "
                        "treated as cb_atac_rc and converted to cb_gex")
    p.add_argument("--puck", default=None,
                   help="headerless puck CSV (bead_bc, x_um, y_um); if given, "
                        "x_um/y_um are left-joined into cell-barcode_coords.csv")
    p.add_argument("--partial", action="store_true",
                   help="chunked mode: output sorted partial_triplets.tsv + "
                        "partial_stats.json; skip aggregation")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.odir, exist_ok=True)
    tmp_dir = args.tmp_dir or args.odir

    valid_cb_set, n_whitelist = load_whitelist(args.wl_fp, args.cb_fmats)

    tmp_path = os.path.join(args.odir, "_tmp_triplets.tsv")
    total_reads = 0

    print("Processing reads...", flush=True)
    with open(tmp_path, "w", buffering=1 << 20) as tmp_fh, \
         Pool(args.workers,
              initializer=_init_worker,
              initargs=(valid_cb_set,)) as pool:
        for batch_results in pool.imap(
            process_batch,
            fastq_pair_reader(args.r1_fp, args.r2_fp, args.batch_size),
            chunksize=1,
        ):
            total_reads += args.batch_size
            for cb, sb, umi in batch_results:
                tmp_fh.write(f"{cb}\t{sb}\t{umi}\n")
            if total_reads % 5_000_000 < args.batch_size:
                print(f"  {total_reads:,} reads processed...", flush=True)

    if args.partial:
        sorted_path = os.path.join(args.odir, "partial_triplets.tsv")
        print("Sorting triplets...", flush=True)
        subprocess.run(
            ["sort", "--buffer-size=2G", "--parallel=4",
             f"-T{tmp_dir}", "-t\t", "-k1,1", "-k2,2", "-k3,3",
             "-o", sorted_path, tmp_path],
            check=True,
        )
        os.remove(tmp_path)
        with open(os.path.join(args.odir, "partial_stats.json"), "w") as fh:
            json.dump({"total_reads": total_reads, "n_whitelist": n_whitelist}, fh)
        print(f"Partial done: {total_reads:,} reads → {sorted_path}", flush=True)
    else:
        puck = _load_puck(args.puck)
        aggregate_and_write(tmp_path, args.odir, tmp_dir, total_reads, n_whitelist, puck)


if __name__ == "__main__":
    main()
