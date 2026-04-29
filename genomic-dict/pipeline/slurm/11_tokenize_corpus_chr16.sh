#!/bin/bash
# Stage 11 — tokenize 17k corpus, emit chr16 tokens + file embeddings. ~30-45 min expected.
#SBATCH --job-name='tokenize_chr16'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/11_tokenize_corpus_chr16_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/11_tokenize_corpus_chr16_%j.err'
#SBATCH --mem='16GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=06:00:00

set -euo pipefail
WORKSPACE="/home/sp5fd/spatial-region-features"
echo "=== stage 11 tokenize_corpus_chr16 ==="
echo "host: $(hostname); date: $(date)"
source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"
cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/11_tokenize_corpus_chr16.py
