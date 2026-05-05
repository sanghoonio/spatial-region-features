#!/bin/bash
# Stage 14 — featured-interval bigwig signal slices.
# Reads featured_files / featured_intervals + bigwig_manifest, samples N bins
# per (file, interval). Runtime <5 min when bigwigs are local.
#SBATCH --job-name='featured_signal'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/09_featured_signal_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/09_featured_signal_%j.err'
#SBATCH --mem='4GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=00:30:00

set -euo pipefail
WORKSPACE="/home/sp5fd/spatial-region-features"
echo "=== stage 09 featured_signal ==="
echo "host: $(hostname); date: $(date)"
source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"
cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/09_featured_signal.py
