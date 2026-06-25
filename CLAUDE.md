# slide-tags pipeline

Top-level launcher for CB↔SB matching (pair) and spatial mapping (map).

## Sample sheet

Stored in Google Sheets — no local CSV. Sheet ID: `1ctXseMOjzDodmE581xLaFzc_dJ5I4svlMGTofDjNRUU`

Schema: `pair, map, xid, x_spl, s_dst, r_dst, puck_id, pair_txt, r1_fastq, r2_fastq, whitelist_tsv, puck_csv, odir`

- `pair` / `map` — `TRUE` = submit this step; `FALSE` = skip.
- Sample name: `{xid}-{x_spl}` when x_spl is set (e.g. xBO368-1), else `{xid}`.
- `whitelist_tsv` — CellRanger `barcodes.tsv.gz` (input to pair step).
- `s_dst` / `r_dst` — spatial / RNA dataset ids; used to build the output dataset folder `{xid}[-{x_spl}]_S-{s_dst}_T-{r_dst}`.
- The output base is built from `xid`: `<experiments>/{xid}/slide-tags/{dataset}/`. The `odir` column is carried in the sheet but ignored by the launcher.
- `pair_txt` — (optional) whitelist source when `pair=FALSE`. A path to a `cell-barcode_coords.csv` file (or its directory; legacy `df_whitelist.txt` also accepted). Blank = auto-select most recent prior pair run under the dataset base.

To sync manually: `python3 sheet_sync.py pull` or `python3 sheet_sync.py push`

## Running

```bash
source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh && conda activate dropme
python slide-tags.py              # pair (pair=FALSE) then map (map=FALSE)
python slide-tags.py --pair-only  # pair only
python slide-tags.py --map-only   # map only
python slide-tags.py --dry-run    # print sbatch commands, do not submit
```

## Output structure

```
<experiments>/{xid}/slide-tags/{xid}[-{x_spl}]_S-{s_dst}_T-{r_dst}/{YYYY-MM-DD_HH-MM-SS}/
  cell-barcode_coords.csv   (pair: CB↔SB match table + puck x/y; comma-separated)
  cell_coords.csv           (map: final per-cell table; comma-separated)
  {puck_basename}.csv       (copy of puck used)
  pair_whitelist.txt        (path to cell-barcode_coords.csv used)
  logs/%j_pair_{spl}.out / %j_map_{spl}.out
```

Both pair and map write directly into the fork dir (no `pair/` or `map/` subdirs).

One timestamped fork dir is created per sample (row). When both pair and map are TRUE, map is
submitted with `--dependency=afterok` on the pair job. When only map=TRUE, the
launcher uses `pair_txt` if set, otherwise picks the most recent prior pair run
under the dataset base.

## Subdirectories

- `pair/` — CB↔SB matcher (from `tags/`); runs `cb-sb_match_sub.py`
- `map/`  — `map_cells.py` + `cell_fit_metrics.R`; runs `tags_mapping.sbatch`
