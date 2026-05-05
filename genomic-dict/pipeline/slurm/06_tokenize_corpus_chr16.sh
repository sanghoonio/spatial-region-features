#!/bin/bash
# Stage 06 — tokenize 17k corpus, emit chr16 token lists, file embeddings,
# and the file UMAP (former stage 09 folded in). ~30-45 min expected.
#SBATCH --job-name='tokenize_chr16'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/06_tokenize_corpus_chr16_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/06_tokenize_corpus_chr16_%j.err'
#SBATCH --mem='16GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=06:00:00

set -euo pipefail
WORKSPACE="/home/sp5fd/spatial-region-features"
echo "=== stage 06 tokenize_corpus_chr16 ==="
echo "host: $(hostname); date: $(date)"
source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"
cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/06_tokenize_corpus_chr16.py
