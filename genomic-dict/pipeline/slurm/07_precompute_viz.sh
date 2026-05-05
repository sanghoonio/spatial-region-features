#!/bin/bash
# SLURM wrapper for stage 07 (precompute viz artifacts).
#
# Submit via yoke REPL with absolute path:
#
#     sbatch /home/sp5fd/spatial-region-features/genomic-dict/pipeline/slurm/07_precompute_viz.sh
#
# UMAP + kNN on ~36k chr16 regions × 100-dim pretrained embeddings.
# Runtime: ~3–5 min.

#SBATCH --job-name='precompute_viz'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/07_precompute_viz_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/07_precompute_viz_%j.err'
#SBATCH --mem='16GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=00:30:00

set -euo pipefail

WORKSPACE="/home/sp5fd/spatial-region-features"

echo "=== stage 07 precompute_viz ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo

source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"

cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/07_precompute_viz.py
