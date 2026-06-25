#!/bin/bash
#SBATCH -J cbsbmerge
#SBATCH -p short
#SBATCH -t 0-02:00:00
#SBATCH -c 10
#SBATCH --mem-per-cpu=8G
#SBATCH -o slurm-%j_cb-sb_merge.out
#SBATCH -e slurm-%j_cb-sb_merge.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=borr@broadinstitute.org

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate dropme

manifest=$1
odir=$2
spl=$3
chunk_scratch=$4
puck=${5:-}

SCRIPT_DIR="${CBSB_PIPE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

python "$SCRIPT_DIR/cb-sb_match_merge.py" \
    "$manifest" "$odir" "$spl" "$chunk_scratch" "$puck"
