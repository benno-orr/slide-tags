#!/bin/bash
#SBATCH -J cbsbmach
#SBATCH -p short
#SBATCH -t 0-08:00:00
#SBATCH -c 20
#SBATCH --mem-per-cpu=4G
# default log target (relative to submit cwd); cb-sb_match_sub.py overrides with -o/-e
#SBATCH -o slurm-%j_cb-sb_match.out
#SBATCH -e slurm-%j_cb-sb_match.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=borr@broadinstitute.org

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate dropme

r1_path=$1
r2_path=$2
wl_path=$3
o_dir=$4
shift 4

mkdir -p "$o_dir"

# cb-sb_match.py lives beside this wrapper, but under sbatch SLURM copies the batch
# script to a spool dir, so BASH_SOURCE is unreliable. The launcher passes the real
# pipeline dir via CBSB_PIPE_DIR; fall back to BASH_SOURCE for direct `bash` runs.
SCRIPT_DIR="${CBSB_PIPE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

python "$SCRIPT_DIR/cb-sb_match.py" \
    "$r1_path" "$r2_path" "$wl_path" "$o_dir" \
    --workers 20 \
    --batch-size 500000 \
    --tmp-dir /n/scratch/users/b/beo703 \
    "$@"

python "$SCRIPT_DIR/saturation.py" \
    "$o_dir/cell-barcode_coords.csv" \
    "$o_dir" \
    "$(basename "$o_dir")"

# Normally invoked by cb-sb_match_sub.py. Direct submission example:
# sbatch cb-sb_match_sbatch.sh <R1.fastq.gz> <R2.fastq.gz> <whitelist> <odir> [--cb-fmats <arc_cb_fmats.csv>]
