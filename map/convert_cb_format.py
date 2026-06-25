#!/usr/bin/env python3
"""
Convert a cell barcode list between 10x ARC formats.

Reads a whitelist (one barcode per line), looks each up in arc_cb_fmats.csv,
and writes the converted barcodes to stdout or a file.

Usage:
    python convert_cb_format.py <whitelist> --from <col> --to <col> [--out <file>]

Columns available in arc_cb_fmats.csv:
    cb_gex       GEX barcode (as reported by CellRanger GEX / scRNA-seq)
    cb_atac      ATAC barcode (raw)
    cb_atac_rc   Reverse complement of ATAC barcode
    cb_arc       GEX barcode with -1 suffix

Example:
    # Convert atac-rc whitelist to gex barcodes
    python convert_cb_format.py xBO203_whitelist.txt --from cb_atac_rc --to cb_gex --out xBO203_wl_gex.txt
"""

import argparse
import csv
import gzip
import sys

CB_FMATS = "/n/data1/hms/scrb/chen/lab/bco/misc/arc_cb_fmats.csv"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("whitelist", help="Input whitelist (one barcode per line, plain or .gz)")
    p.add_argument("--from", dest="from_col", required=True,
                   choices=["cb_gex", "cb_atac", "cb_atac_rc", "cb_arc"],
                   help="Column in arc_cb_fmats.csv that the input barcodes match")
    p.add_argument("--to", dest="to_col", required=True,
                   choices=["cb_gex", "cb_atac", "cb_atac_rc", "cb_arc"],
                   help="Column in arc_cb_fmats.csv to convert to")
    p.add_argument("--cb-fmats", default=CB_FMATS,
                   help=f"Path to arc_cb_fmats.csv (default: {CB_FMATS})")
    p.add_argument("--out", default=None,
                   help="Output file (default: stdout)")
    args = p.parse_args()

    print(f"Loading {args.cb_fmats}...", file=sys.stderr, flush=True)
    mapping = {}
    with open(args.cb_fmats) as fh:
        for row in csv.DictReader(fh):
            key = row[args.from_col].strip('"')
            val = row[args.to_col].strip('"')
            mapping[key] = val
    print(f"  {len(mapping):,} entries", file=sys.stderr, flush=True)

    open_fn = gzip.open if args.whitelist.endswith(".gz") else open
    out_fh = open(args.out, "w") if args.out else sys.stdout

    found = skipped = 0
    with open_fn(args.whitelist, "rt") as fh:
        for line in fh:
            bc = line.strip()
            if bc.endswith("-1"):
                bc = bc[:-2]
            converted = mapping.get(bc)
            if converted is None:
                skipped += 1
                continue
            out_fh.write(converted + "\n")
            found += 1

    if args.out:
        out_fh.close()

    print(f"Converted: {found:,}  Not found: {skipped:,}", file=sys.stderr, flush=True)
    if args.out:
        print(f"Written to {args.out}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
