# tags_mapping

Stand-alone slide-tags spatial-mapping pipeline. Three Python scripts; every
input and output path is supplied on the command line — no project-specific
hard-coded paths.

```
df_whitelist.txt + puck.csv
         │
         ▼
    classify.py        ──► {sample}_classify.tsv         (per-cell DBSCAN + KDE topo peaks)
         │
         ▼
  cell_profiles.py     ──► {stamp}_{sample}_cell_profiles.tsv
                              (per-cell prom rankings, peak1/peak2 leftmost-line
                               log-fold-changes lfc1/lfc2, radial cumulative-UMI)
         │
         ▼
   lfc_sweep.py        ──► {stamp}_{sample}_lfc_{spatial_sweep,stats,hexhist}.png
                              (singlet/multiplet calling vs lfc threshold)
```

## Inputs

- **Whitelist TSV** — Slide-tags `df_whitelist.txt` with columns
  `cell_bc_10x, bead_bc, nUMI`.
- **Puck CSV** — no header, columns `bead_bc, x_um, y_um`. Homopolymer
  contaminant beads (`GGGGGGGGGGGGGG`, `CCCCCCCCCCCCCC`) are dropped
  automatically.

## Step 1 — classify

```bash
python classify.py \
    --sample S \
    --whitelist /path/df_whitelist.txt \
    --puck /path/puck.csv \
    --out-dir /path/out
```

Writes `{out}/S_classify.tsv`. Multi-sample mode: pass `--manifest m.tsv`
(columns `sample, whitelist, puck`) instead of `--sample/--whitelist/--puck` —
this also emits a `{stamp}_classify_summary.{tsv,png}`.

Key tunables: `--eps 50 75 100` (DBSCAN eps escalation), `--min-samples 10`,
`--sigma 40` (KDE bandwidth µm), `--grid 120`, `--cutoffs 0.0001 0.0002 0.0005 0.001`
(topo-prominence cutoffs for singlet/doublet calling), `--ncores 16`.

## Step 2 — cell_profiles

```bash
python cell_profiles.py \
    --sample S \
    --whitelist /path/df_whitelist.txt \
    --puck /path/puck.csv \
    --classify-tsv /path/out/S_classify.tsv \
    --out-dir /path/out
```

Writes `{out}/{stamp}_S_cell_profiles.tsv`. The classify TSV supplies the
cell-barcode list; cells are re-DBSCAN'd with **raw-UMI-weighted** sample
weights (no /1000 normalization) at `--eps 50 --min-samples 8`. For each cell
records:
- `c1_umi, c2_umi, c0_umi` and centroids `cx1, cy1, cx2, cy2`
- top `--nprom 100` KDE prominences (`prom001..prom100`)
- `peak1_obs-prj_lfc`, `peak2_obs-prj_lfc` — log10 obs-vs-projected fold change
  for peaks 1 & 2, where the projection is a 4-pt leftmost-line fit of
  `log10(prom)` vs rank starting at rank 3 (peaks 1-2 excluded), first window
  with R² ≥ 0.80
- cumulative-UMI vs radius from the c1 centroid (`cumumi_r020..cumumi_r400`)

## Step 3 — lfc_sweep

```bash
python lfc_sweep.py \
    --sample S \
    --profiles-tsv /path/out/{stamp}_S_cell_profiles.tsv \
    --out-dir /path/out
```

Writes three PNGs (timestamped): `lfc_spatial_sweep`, `lfc_stats`,
`lfc_hexhist`. Sweeps thresholds 0.10..0.85 (step 0.05) on `lfc1 = peak1_lfc`;
classifies cells per threshold as
- **mapped** if `lfc1 ≥ t`
- **singlet** if `lfc1 ≥ t AND lfc2 < t`
- **multiplet** if `lfc1 ≥ t AND lfc2 ≥ t`

The dashed green line in figure B is the expected multiplet rate at
`--expected-mult-per-1k 0.8` (%/1k cells × #cells in the input).

## Batch submission (table-driven)

Same pattern as `bco/scripts/10x/cellranger-count`. Edit `tags_mapping.csv`
(columns `run,sample,whitelist,puck,out_dir`), set `run=TRUE` on the rows you
want, then:

```bash
bash /n/data1/hms/scrb/chen/lab/bco/pipelines/tags_mapping/tags_mapping.sh
```

Each `TRUE` row submits one sbatch job (`tagsmap_{sample}`) on the `short`
partition that runs all three steps (`classify → cell_profiles → lfc_sweep`)
for that sample. **`out_dir` should be the spatial run directory** (e.g.
`xBO368/spatial/260530_141031`); the sbatch appends `/spl` automatically so
all outputs land in `out_dir/spl/`. Rows where `run` is not literally `TRUE`
are skipped.

## Notes

- `cell_profiles.py` re-uses `classify.PARAMS / load_cb_bead_table /
  kde_peaks_full` directly — keep the three scripts together.
- All randomness is seeded only where it matters (none in this pipeline);
  re-running on identical inputs is deterministic up to joblib ordering.
- For multi-sample fan-out wrap each step in your own shell/sbatch loop —
  no project manifest paths are baked in.
