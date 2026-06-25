# CLAUDE.md — `tags` pipeline

Reusable pipeline for Slide-tags **SB↔CB matching** and **SB library-complexity** projection.
Distilled from analysis aBO089 (canonical Python path). Runs on the **O2 SLURM** cluster.

Two concerns, nothing else (no spatial positioning, no RNA saturation):

1. **SB↔CB matching** — match 10x cell barcodes (CB, R1) to spatial/bead barcodes (SB, R2),
   filtered to a CB whitelist → per-sample `cell-barcode_coords.csv`.
2. **Library complexity** — analytic sequencing-saturation curve + Michaelis-Menten projection
   of the SB library, from the matching output.

Output, log, and plot locations are **caller-supplied** — the pipeline writes nothing inside
its own directory. Point it at any analysis's output tree.

## Layout

```
cb-sb_match.py          workhorse matcher (R1 CB+UMI × R2 SB, whitelist-filtered)
cb-sb_match_sub.py      launcher: reads a file_paths.csv, sbatch one job per run==TRUE row
cb-sb_match_sbatch.sh   O2 sbatch wrapper (env dropme, short, 20c/4GB, 8h)
convert_cb_format.py    ATAC→GEX CB conversion (misc/arc_cb_fmats.csv)
saturation.py           per-sample SB saturation + M-M complexity projection
run_saturation_all.py   runs saturation.py over every matched sample
config/file_paths.csv    sample-table template (worked examples, all run=FALSE)
config/reuploads.csv     append-only log of every GitHub push (written by push.sh)
push.sh                  commit + push to benno-orr/tags, logging the reupload
```

## Running

```bash
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh && conda activate dropme

# SB↔CB matching + saturation curve — one sbatch job per row with run==TRUE
python cb-sb_match_sub.py [--odir <ANALYSIS>/outputs/restored] \
                          [--csv config/file_paths.csv] \
                          [--log-dir <RUN>/logs] [--dry-run]
```

- Each sbatch job runs matching then immediately runs `saturation.py` on the output —
  `saturation_curve.csv`, `saturation_stats.csv`, and plots land in the same per-sample dir.
  `run_saturation_all.py` is still available for re-running saturation standalone (e.g. with
  a custom `--plt-id` or `--plots-dir`).
- Output base is per-row from the CSV `odir` column, falling back to `--odir` (aliases
  `--odir-base`, `-odir`) for blank rows — one of the two must supply it. Each invocation
  creates a fresh timestamped run dir `<base>/<YYMMDD_HHMMSS>/`, and per-sample output lands
  in `<run>/<spl>/` — so re-runs never clobber prior output. `--log-dir` defaults to
  `<run>/logs`. The launcher passes `-o/-e` to sbatch (overriding the wrapper defaults).
- Each run dir is stamped with provenance: `pipeline_version.json` (git commit + dirty flag
  of the pipeline checkout, timestamp, csv path) and `file_paths.used.csv` (a snapshot of the
  sample table), so any output traces back to the exact pipeline version that produced it.

## `config/file_paths.csv`

Columns: `spl, run, r1_fq, r2_fq, puck_csv, wl, cb_fmats, odir`

- `r1_fq` — R1 FASTQ (CB pos 1–16, UMI pos 17–28); `r2_fq` — R2 (SB = pos 1–8 + pos 27–32, 14 bp).
- `wl` — CB whitelist (CellRanger `barcodes.tsv.gz` for GEX; `-1` suffix stripped automatically).
- `cb_fmats` — set to `/n/data1/hms/scrb/chen/lab/bco/misc/arc_cb_fmats.csv` **only** for ATAC
  reverse-complement whitelists (triggers on-load CB conversion); blank for plain GEX.
- `odir` — per-row output base (`<odir>/<YYMMDD_HHMMSS>/<spl>/`); blank falls back to `--odir`.
- `puck_csv` is carried for reference; matching does not use it.

The shipped rows are aBO089 examples (all `run=FALSE`); some point at ephemeral `/n/scratch`
paths — they document the format. Add a sample = append a row, set `run=TRUE`.

## Output — `<odir>/<YYMMDD_HHMMSS>/<spl>/`

`cell-barcode_coords.csv` (comma-separated): `CB_SB` (30 bp = 16 bp CB + 14 bp SB), `nUMI` (unique CB,SB,UMI
triplets), `nReads` (raw reads per CB,SB pair), `cell_bc_10x` (16 bp CB), `bead_bc` (14 bp SB).
Plus `saturation_curve.csv`, `saturation_stats.csv`, and a `<plt_id>_saturation-curve_v1_<spl>.{pdf,png}` plot.
The run dir also holds `pipeline_version.json` + `file_paths.used.csv` (provenance, see Running).

## GitHub — `benno-orr/tags`

The pipeline is version-controlled and pushed to `https://github.com/benno-orr/tags.git`.
Reupload with `./push.sh "message"`: it commits, pushes, and appends a row
(`timestamp, commit, branch, message`) to `config/reuploads.csv`. The commit SHA logged there
is the same one the launcher writes into each run dir's `pipeline_version.json`.

## Environment

Conda env `dropme` (`/n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh`). Matching jobs:
`short`, 20 cores × 4 GB, 8 h; sort temp under `/n/scratch/users/b/beo703`.

## Notes

- All scripts resolve siblings via their own location, so the pipeline is relocatable.
- Analysis **aBO134** holds an independent copy of this code (its own `file_paths.csv` + outputs);
  keep the two in sync manually if you change one.
- pltBO ids default to placeholder `pltBOxx`; registering rows in `tables/plots.csv` is separate.
