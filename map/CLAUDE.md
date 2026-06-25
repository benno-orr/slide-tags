# CLAUDE.md — `tags` pipeline

This repo holds two co-dependent pipelines for Slide-tags spatial mapping:

1. **SB↔CB matching** (`cb-sb_match*`) — match 10x cell barcodes to spatial bead barcodes → `cell-barcode_coords.csv`
2. **Spatial mapping** (`map_cells.py`) — DBSCAN + KDE + topographic prominence → per-cell spatial-profile data (`cell_coords.csv`)

Both run on the **O2 SLURM** cluster.

---

# Part 1 — SB↔CB Matching

Reusable pipeline for Slide-tags **SB↔CB matching** and **SB library-complexity** projection.
Distilled from analysis aBO089 (canonical Python path).

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

## Environment

Conda env `dropme` (`/n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh`). Matching jobs:
`short`, 20 cores × 4 GB, 8 h; sort temp under `/n/scratch/users/b/beo703`.

## Notes

- All scripts resolve siblings via their own location, so the pipeline is relocatable.
- Analysis **aBO134** holds an independent copy of this code (its own `file_paths.csv` + outputs);
  keep the two in sync manually if you change one.
- pltBO ids default to placeholder `pltBOxx`; registering rows in `tables/plots.csv` is separate.

---

# Part 2 — Spatial Mapping

Stand-alone Slide-tags spatial-mapping pipeline. `map_cells.py` does all per-cell
work (DBSCAN + KDE + spatial profiles) in a single pass.
(Formerly two scripts, `classify.py` + `cell_profiles.py`, now merged.)

Pipeline order:
```
map_cells.py       → cell_coords.csv   (all per-cell data; comma-separated)
cell_fit_metrics.R → merges fit metrics into cell_coords.csv (in place)
```

## Running

**Single sample, direct execution** (on a compute node):
```bash
cd /n/data1/hms/scrb/chen/lab/bco/pipelines/slide-tags/map

python map_cells.py --sample S --whitelist /path/cell-barcode_coords.csv --out-dir /path/out
```
`--whitelist` auto-detects delimiter and reads x/y straight from cell-barcode_coords.csv;
pass `--puck` only for a legacy coord-free `df_whitelist.txt`.

**Batch via SLURM** (short partition, 20 cores, 128 GB, 12 h):
1. Edit `tags_mapping.csv` — set `run=TRUE` on rows to submit (columns: `run,sample,whitelist,puck,out_dir`). Set `out_dir` to the run directory; outputs are written directly into it.
2. Run:
```bash
bash /n/data1/hms/scrb/chen/lab/bco/pipelines/tags_mapping/tags_mapping.sh
```
Each `TRUE` row submits one job (`tagsmap_{sample}`).

## Key tunables

**map_cells.py** (collects per-cell DATA only — emits no singlet/doublet/unmapped call)
- `--eps 50` — fixed DBSCAN eps (µm); only the first value is used
- DBSCAN runs under **two weightings**, both emitted (columns prefixed `raw_`/`norm_`):
  - `--min-samples 10 9 8` — raw mode: descending escalation, first tier with ≥1 cluster
  - `--norm-to 1000` — norm mode weights `w = u/Σu × norm_to`
  - `--norm-ms-scan 50 1` — norm mode: integer min_samples scanned high→low, highest tier yielding **≥2 clusters** wins (NaN if none)
  - the largest 8 clusters per mode are recorded as c1..c8, with a c9plus pool + c0 noise
- `--sigma 40` — KDE bandwidth (µm); `--grid 120` — KDE grid resolution
- `--nprom 1000` — number of KDE prominence/height peaks recorded per cell
- `--radii 1 2 ... 400` — radii for cumulative-UMI columns (from the raw c1 centroid)
- `--ncores 16` — joblib parallelism

## Output naming

- `map_cells.py` writes `{out_dir}/cell_coords.csv` (comma-separated): one row per cell with the DBSCAN (raw/norm) cluster stats, KDE peak heights/prominences, leftmost-line fold changes, and cumulative-UMI profile. `cell_fit_metrics.R` merges fit metrics into it in place.

## Algorithm summary

`map_cells.py` per-cell (DATA only — no classification), single pass:
- **DBSCAN** at fixed eps under two weightings — **raw** (per-bead UMI, descending min_samples escalation, first tier ≥1 cluster) and **norm** (`u/Σu·norm_to`, integer min_samples scanned 50→1, highest tier yielding ≥2 clusters). Each emits the largest 8 clusters (c1..c8), a c9plus pool, and c0 noise. cx1/cy1/cx2/cy2 alias the raw c1/c2 centroids. No singlet/doublet/unmapped status is assigned — thresholding is left to downstream consumers.
- **KDE**: fixed-σ Gaussian-sum on a bbox grid → topographic prominence (descending-flood union-find) → top-8 peak heights/prominences + background stats, the full top-`nprom` (default 1000) prominence and height rankings, and each mode's c1→nearest-peak prom/height.
- **Leftmost-line fold changes**: `peak{1,2}fclr` — fit log10(prom) vs rank from rank 3 (peaks 1-2 excluded), prom > 0.1 excluded, 4-pt window R² ≥ 0.80; metric = log10(promK) − projected line. `peak{1,2}fclr_h` is the same fit on raw height (peaks re-ranked by height, prom > 0.1 still excluded). Each is fit in linear-rank and log-rank space (`_logr`), recording R²/n-points/window-span and which space fit better.
- **Cumulative-UMI** vs radius (1–400 µm) from the raw c1 centroid.

(Near-)homopolymer beads — any base in ≥12/14 positions (within 2 mismatches of a pure homopolymer of A/C/G/T) — and beads with no coords are dropped at load time, before DBSCAN and KDE.
