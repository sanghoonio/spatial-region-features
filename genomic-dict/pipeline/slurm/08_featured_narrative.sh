#!/bin/bash
# SLURM wrapper for stage 08 (featured-narrative data prep).
#
# Submit via yoke REPL with absolute path:
#
#     sbatch /home/sp5fd/spatial-region-features/genomic-dict/pipeline/slurm/08_featured_narrative.sh
#
# Tokenizes ~22 BED files (18 featured + 4 mystery) via bbcache, slices by
# featured intervals. Runtime ~5 min. Uses BBCLIENT_CACHE.

#SBATCH --job-name='featured_narrative'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/08_featured_narrative_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/08_featured_narrative_%j.err'
#SBATCH --mem='8GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=00:30:00

set -euo pipefail

WORKSPACE="/home/sp5fd/spatial-region-features"

echo "=== stage 08 featured_narrative ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo

source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"

cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/08_featured_narrative.py
