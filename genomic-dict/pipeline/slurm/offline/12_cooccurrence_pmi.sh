#!/bin/bash
# Stage 12 — cooccurrence NPMI per stratum + per-token marginals. Chunked
# sparse matmul over ~17k files × ~36k chr16 tokens × 18 strata. ~30 min.
#SBATCH --job-name='cooc_pmi'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/12_cooccurrence_pmi_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/12_cooccurrence_pmi_%j.err'
#SBATCH --mem='32GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=02:00:00

set -euo pipefail
WORKSPACE="/home/sp5fd/spatial-region-features"
echo "=== stage 12 cooccurrence_pmi ==="
echo "host: $(hostname); date: $(date)"
source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"
cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/12_cooccurrence_pmi.py
