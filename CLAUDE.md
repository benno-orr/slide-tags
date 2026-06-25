# slide-tags pipeline

Top-level launcher for CB↔SB matching (pair) and spatial mapping (map).

## Sample sheet

Stored in Google Sheets — no local CSV. Sheet ID: `1ctXseMOjzDodmE581xLaFzc_dJ5I4svlMGTofDjNRUU`

Schema (13 cols): `pair, map, xid, x_spl, s_seq, r_seq, puck_id, pair_txt, r1_fastq, r2_fastq, whitelist_tsv, puck_csv, odir`

- `pair` / `map` — `TRUE` = submit this step; `FALSE` = skip.
- Sample name: `{xid}-{x_spl}` when x_spl is set (e.g. xBO368-1), else `{xid}`.
- `whitelist_tsv` — CellRanger `barcodes.tsv.gz` (input to pair step).
- `odir` — experiment's slide-tags output root (e.g. `.../xBO180/slide-tags`).
- `pair_txt` — (optional) whitelist source when `pair=FALSE`. A path to a `df_whitelist.txt` file (or its directory). Blank = auto-select most recent prior pair run under `odir`.

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
{odir}/{YYMMDD_HHMMSS}_{spl}/
  pair/{spl}/df_whitelist.txt
  map/mapping.tsv
  map/{spl}_cb_sb_umi_xy.csv     (merged CB-SB-nUMI-x_um-y_um table fed to classify)
  map/{spl}_fc_*.png
  map/{puck_basename}.csv        (copy of puck used)
  map/pair_whitelist.txt         (path to df_whitelist.txt used)
  logs/%j_pair_{spl}.out / %j_map_{spl}.out
```

One fork dir is created per sample (row). When both pair and map are TRUE, map is
submitted with `--dependency=afterok` on the pair job. When only map=TRUE, the
launcher uses `pair_txt` if set, otherwise picks the most recent prior pair run
under `odir`.

## Subdirectories

- `pair/` — CB↔SB matcher (from `tags/`); runs `cb-sb_match_sub.py`
- `map/`  — classify + cell_profiles + lfc_sweep (from `tags_mapping/`); runs `tags_mapping.sbatch`
