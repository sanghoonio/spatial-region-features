#!/bin/bash
# Stage 05 — sequence-based intrinsic features. ~1 hr.
#SBATCH --job-name='extract_sequence'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/05_extract_sequence_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/05_extract_sequence_%j.err'
#SBATCH --mem='8GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=02:00:00

set -euo pipefail
WORKSPACE="/home/sp5fd/spatial-region-features"
echo "=== stage 05 extract_sequence ==="
echo "host: $(hostname); date: $(date)"
source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"
cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/05_extract_sequence.py
