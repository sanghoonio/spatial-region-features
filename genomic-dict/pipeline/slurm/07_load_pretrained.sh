#!/bin/bash
# SLURM wrapper for stage 07 (load pretrained R2V universe).
#
# Submit via yoke REPL with absolute path:
#
#     sbatch /home/sp5fd/spatial-region-features/genomic-dict/pipeline/slurm/07_load_pretrained.sh
#
# Downloads the databio/r2v-encode-hg38 model from HuggingFace, extracts
# vocabulary regions + embeddings, intersects with SCREEN metadata from
# stage 02. No training, no 85k-BED corpus needed.
# Runtime: ~5–15 minutes (HF download + overlap compute).

#SBATCH --job-name='load_pretrained'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/07_load_pretrained_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/07_load_pretrained_%j.err'
#SBATCH --mem='16GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=01:00:00

set -euo pipefail

WORKSPACE="/home/sp5fd/spatial-region-features"

echo "=== stage 07 load_pretrained ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo

source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"

cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/07_load_pretrained.py
