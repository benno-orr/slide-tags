#!/bin/bash
#SBATCH -c 20
#SBATCH --mem=128G
#SBATCH -p short
#SBATCH -t 0-12:00
#SBATCH -o /home/beo703/reports/%N_%j.out
#SBATCH -e /home/beo703/reports/%N_%j.err
#SBATCH --mail-user='benno_orr@g.harvard.edu'
#SBATCH --mail-type=ALL

set -euo pipefail

# Required env vars exported by tags_mapping.sh: sample, whitelist, puck, out_dir
# Outputs (cell_coords.csv, plots, provenance) are written directly into out_dir.
: "${sample:?sample not set}"
: "${whitelist:?whitelist not set}"
: "${puck:?puck not set}"
: "${out_dir:?out_dir not set}"

pdir="${SLIDE_TAGS_MAP_DIR:-/n/data1/hms/scrb/chen/lab/bco/pipelines/slide-tags/map}"
mkdir -p "$out_dir"
cd "$pdir"

# record provenance (skip the copy when puck already lives in out_dir, e.g. an
# in-place rerun — cp onto itself would abort under set -e)
[ "$puck" -ef "$out_dir/$(basename "$puck")" ] || cp "$puck" "$out_dir/$(basename "$puck")"
echo "$whitelist" > "$out_dir/pair_whitelist.txt"

echo "=== map_cells $sample ==="
python3 map_cells.py \
    --sample "$sample" \
    --whitelist "$whitelist" \
    --puck "$puck" \
    --out-dir "$out_dir" \
    --ncores "${SLURM_CPUS_PER_TASK:-20}"

profiles_tsv="$out_dir/cell_coords.csv"

echo "=== cell_fit_metrics $sample ==="
# Per-cell gamma + intden fit metrics, merged in place into cell_coords.csv.
Rscript="${RSCRIPT:-/n/app/R/4.4.2-gcc-14.2.0/bin/Rscript}"
"$Rscript" cell_fit_metrics.R "$profiles_tsv"

echo "=== done $sample ==="
