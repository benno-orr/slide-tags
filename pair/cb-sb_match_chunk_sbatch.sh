#!/bin/bash
#SBATCH -J cbsbchunk
#SBATCH -p short
#SBATCH -t 0-00:30:00
#SBATCH -c 10
#SBATCH --mem-per-cpu=4G
#SBATCH -o slurm-%j_cb-sb_chunk.out
#SBATCH -e slurm-%j_cb-sb_chunk.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=borr@broadinstitute.org

source /n/app/conda/miniforge3/24.11.3-0/etc/profile.d/conda.sh
conda activate dropme

r1_chunk=$1
r2_chunk=$2
wl_path=$3
chunk_odir=$4
shift 4

SCRIPT_DIR="${CBSB_PIPE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

mkdir -p "$chunk_odir"

python "$SCRIPT_DIR/cb-sb_match.py" \
    "$r1_chunk" "$r2_chunk" "$wl_path" "$chunk_odir" \
    --workers 10 \
    --batch-size 500000 \
    --tmp-dir /n/scratch/users/b/beo703 \
    --partial \
    "$@"

# Free scratch space — chunk FASTQs no longer needed after matching
rm -f "$r1_chunk" "$r2_chunk"
