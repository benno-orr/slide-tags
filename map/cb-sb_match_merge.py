#!/usr/bin/env python3
"""
Merge sorted partial triplets from chunk jobs → df_whitelist.txt + saturation.

Usage:
    python cb-sb_match_merge.py <manifest> <odir> <spl> [chunk_scratch] [puck]

    manifest      — file listing one partial_triplets.tsv path per line
    odir          — output dir (where cell-barcode_coords.csv is written)
    spl           — sample name
    chunk_scratch — scratch dir to delete after completion (optional; "" to skip)
    puck          — headerless puck CSV (bead_bc, x_um, y_um); if given, x_um/y_um
                    are left-joined into cell-barcode_coords.csv (optional)
"""

import csv
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict

PIPE_DIR = os.path.dirname(os.path.abspath(__file__))

manifest     = sys.argv[1]
odir         = sys.argv[2]
spl          = sys.argv[3]
chunk_scratch = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] else None
puck_fp      = sys.argv[5] if len(sys.argv) > 5 else None

os.makedirs(odir, exist_ok=True)


def _load_puck(fp):
    """Headerless puck CSV (bead_bc, x_um, y_um) → {bead_bc: (x_um, y_um)} or None."""
    if not fp or not os.path.exists(fp):
        return None
    coords = {}
    with open(fp, newline="") as fh:
        for row in csv.reader(fh):
            if len(row) >= 3:
                coords[row[0]] = (row[1], row[2])
    return coords


puck = _load_puck(puck_fp)

# Load per-chunk stats
with open(manifest) as fh:
    partial_files = [l.strip() for l in fh if l.strip()]

total_reads = 0
n_whitelist = None
for pf in partial_files:
    stats_path = os.path.join(os.path.dirname(pf), "partial_stats.json")
    if os.path.exists(stats_path):
        with open(stats_path) as fh:
            s = json.load(fh)
            total_reads += s.get("total_reads", 0)
            if n_whitelist is None:
                n_whitelist = s.get("n_whitelist", 0)

print(f"{len(partial_files)} chunks, {total_reads:,} total reads", flush=True)

# Merge-sort all sorted partial files
merged_path = os.path.join(odir, "_merged_triplets.tsv")
print("Merge-sorting partial triplets...", flush=True)
subprocess.run(
    ["sort", "--merge", "--buffer-size=4G", "--parallel=8",
     f"-T{odir}", "-t\t", "-k1,1", "-k2,2", "-k3,3",
     "-o", merged_path] + partial_files,
    check=True,
)

# Aggregate
print("Aggregating...", flush=True)
cb_sb_nreads = defaultdict(int)
cb_sb_numi   = defaultdict(int)
sb_nreads    = defaultdict(int)
matched_cbs  = set()
prev_triplet = None

with open(merged_path) as fh:
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

os.remove(merged_path)

total_matched_reads = sum(cb_sb_nreads.values())
cb_matched          = len(matched_cbs)
cb_sb_unique        = len(cb_sb_nreads)

# Write cell-barcode_coords.csv (left-join puck coords when available)
wl_path = os.path.join(odir, "cell-barcode_coords.csv")
header = ["CB_SB", "nUMI", "nReads", "cell_bc_10x", "bead_bc"]
if puck is not None:
    header += ["x_um", "y_um"]
with open(wl_path, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(header)
    for (cb, sb), nreads in cb_sb_nreads.items():
        cb_sb_str = cb + sb
        bead_bc = cb_sb_str[16:30]
        row = [cb_sb_str, cb_sb_numi[(cb, sb)], nreads, cb_sb_str[0:16], bead_bc]
        if puck is not None:
            x, y = puck.get(bead_bc, ("", ""))
            row += [x, y]
        writer.writerow(row)

# Write reads_per_SB.csv
sorted_sbs = sorted(sb_nreads.keys())
with open(os.path.join(odir, "reads_per_SB.csv"), "w", newline="") as fh:
    writer = csv.writer(fh, quoting=csv.QUOTE_ALL)
    writer.writerow([""] + sorted_sbs)
    writer.writerow(["1"] + [sb_nreads[sb] for sb in sorted_sbs])

# Write matcher_summary.txt
with open(os.path.join(odir, "matcher_summary.txt"), "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    metrics = {
        "hamming_dist": 0,
        "UP_site_matches": total_reads,
        "cell_bender_cb": n_whitelist,
        "cb_matched": cb_matched,
        "cb_matched_reads": total_matched_reads,
        "CB_SB_unique_pairings": cb_sb_unique,
    }
    writer.writerow(list(metrics.keys()))
    writer.writerow(list(metrics.values()))

print(
    f"Done.\n"
    f"  Total reads:          {total_reads:,}\n"
    f"  CB-matched reads:     {total_matched_reads:,} "
    f"({100*total_matched_reads/max(total_reads,1):.1f}%)\n"
    f"  Unique CBs matched:   {cb_matched:,} / {n_whitelist:,}\n"
    f"  Unique CB-SB pairs:   {cb_sb_unique:,}",
    flush=True,
)

# Saturation
subprocess.run(
    [sys.executable, os.path.join(PIPE_DIR, "saturation.py"),
     wl_path, odir, spl],
    check=True,
)

# Cleanup chunk scratch
if chunk_scratch and os.path.exists(chunk_scratch):
    print(f"Cleaning up {chunk_scratch}", flush=True)
    shutil.rmtree(chunk_scratch)
