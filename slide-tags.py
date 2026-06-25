#!/usr/bin/env python3
"""
slide-tags.py — slide-tags pipeline launcher.

Pulls the sample sheet fresh from Google Sheets on every run (no local CSV).
Creates one timestamped run directory per experiment odir, then submits jobs:

  pair — CB↔SB matching  → {run_dir}/pair/{spl}/df_whitelist.txt  [pair == TRUE]
  map  — classify + cell_profiles + lfc_sweep → {run_dir}/map/    [map  == TRUE]

When both pair and map are TRUE for a sample, map is submitted with
--dependency=afterok on the pair job so it waits for the whitelist.

Usage:
    python slide-tags.py [--dry-run] [--pair-only | --map-only]
"""

import argparse, glob, os, subprocess, sys
from datetime import datetime

PIPE_DIR    = os.path.dirname(os.path.abspath(__file__))
PAIR_DIR    = os.path.join(PIPE_DIR, "pair")
MAP_DIR     = os.path.join(PIPE_DIR, "map")
PAIR_SBATCH = os.path.join(PAIR_DIR, "cb-sb_match_sbatch.sh")
MAP_SBATCH  = os.path.join(MAP_DIR,  "tags_mapping.sbatch")

SHEET_ID  = "1ctXseMOjzDodmE581xLaFzc_dJ5I4svlMGTofDjNRUU"
WORKSHEET = "Sheet1"
GS_KEY    = "/home/beo703/.config/gspread/bco-gs-scripts-3146376045cb.json"


def pull_sheet():
    import gspread
    from google.oauth2.service_account import Credentials
    creds = Credentials.from_service_account_file(
        GS_KEY, scopes=["https://www.googleapis.com/auth/spreadsheets"])
    gc = gspread.authorize(creds)
    ws = gc.open_by_key(SHEET_ID).worksheet(WORKSHEET)
    return ws.get_all_records()


def spl(row):
    x = str(row.get("x_spl", "")).strip()
    return f"{row['xid']}-{x}" if x else str(row["xid"]).strip()


def is_true(val):
    return str(val).strip().upper() == "TRUE"


def sbatch_submit(cmd, dry_run=False):
    """Run sbatch; return job-id string or None."""
    if dry_run:
        print(f"  [dry] {' '.join(cmd)}")
        return None
    res = subprocess.run(cmd, capture_output=True, text=True)
    msg = (res.stdout + res.stderr).strip()
    print(f"  {msg}")
    for tok in reversed(msg.split()):
        if tok.isdigit():
            return tok
    return None


def find_whitelist(odir, name, pair_txt=""):
    """Resolve df_whitelist.txt for a map-only run.

    pair_txt is a path to a df_whitelist.txt file (or the directory holding one).
    Blank falls back to the most recent prior pair run under odir.
    """
    if pair_txt:
        # explicit path to a whitelist file (or its directory)
        if pair_txt.endswith(".txt"):
            return pair_txt
        return os.path.join(pair_txt, "df_whitelist.txt")
    hits = sorted(glob.glob(os.path.join(odir, "*", "pair", name, "df_whitelist.txt")))
    if hits:
        return hits[-1]
    # legacy: single df_whitelist.txt placed directly in odir
    root = os.path.join(odir, "df_whitelist.txt")
    return root if os.path.exists(root) else None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dry-run",   action="store_true", help="print commands, do not submit")
    p.add_argument("--pair-only", action="store_true")
    p.add_argument("--map-only",  action="store_true")
    args = p.parse_args()

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M')}] pulling sample sheet…")
    rows = pull_sheet()
    print(f"  {len(rows)} rows\n")

    stamp    = datetime.now().strftime("%y%m%d_%H%M%S")
    run_dirs = {}   # name → timestamped run dir (one per sample)

    for r in rows:
        name    = spl(r)
        odir    = str(r["odir"]).strip()
        do_pair = is_true(r.get("pair")) and not args.map_only
        do_map  = is_true(r.get("map"))  and not args.pair_only

        if not (do_pair or do_map):
            continue

        # Provision one timestamped run dir per sample
        if name not in run_dirs:
            rdir = os.path.join(odir, f"{stamp}_{name}")
            run_dirs[name] = rdir
            if not args.dry_run:
                os.makedirs(os.path.join(rdir, "pair"), exist_ok=True)
                os.makedirs(os.path.join(rdir, "map"),  exist_ok=True)
                os.makedirs(os.path.join(rdir, "logs"), exist_ok=True)
                print(f"run dir: {rdir}")
        rdir = run_dirs[name]

        pair_job_id = None

        # ── pair ─────────────────────────────────────────────────────────────
        if do_pair:
            pair_out = os.path.join(rdir, "pair", name)
            if not args.dry_run:
                os.makedirs(pair_out, exist_ok=True)
            print(f"pair  {name}")
            cmd = [
                "sbatch", "-J", f"pair_{name}",
                f"--export=ALL,CBSB_PIPE_DIR={PAIR_DIR}",
                "-o", os.path.join(rdir, "logs", f"%j_pair_{name}.out"),
                "-e", os.path.join(rdir, "logs", f"%j_pair_{name}.err"),
                PAIR_SBATCH,
                str(r["r1_fastq"]), str(r["r2_fastq"]),
                str(r["whitelist_tsv"]), pair_out,
            ]
            pair_job_id = sbatch_submit(cmd, dry_run=args.dry_run)

        # ── map ──────────────────────────────────────────────────────────────
        if do_map:
            map_out = os.path.join(rdir, "map")
            puck    = str(r["puck_csv"]).strip()

            if do_pair:
                # whitelist will be written by the pair job above
                wl = os.path.join(rdir, "pair", name, "df_whitelist.txt")
            else:
                pair_txt = str(r.get("pair_txt", "")).strip()
                wl = find_whitelist(odir, name, pair_txt)
                if wl is None:
                    print(f"map   {name}  SKIP: no df_whitelist.txt under {odir} — run pair first")
                    continue

            print(f"map   {name}  wl={wl}")
            cmd = [
                "sbatch", "-J", f"map_{name}",
                f"--export=ALL,sample={name},whitelist={wl},puck={puck},"
                f"out_dir={map_out},SLIDE_TAGS_MAP_DIR={MAP_DIR}",
                "-o", os.path.join(rdir, "logs", f"%j_map_{name}.out"),
                "-e", os.path.join(rdir, "logs", f"%j_map_{name}.err"),
            ]
            if pair_job_id:
                cmd += ["--dependency", f"afterok:{pair_job_id}"]
            cmd.append(MAP_SBATCH)
            sbatch_submit(cmd, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
