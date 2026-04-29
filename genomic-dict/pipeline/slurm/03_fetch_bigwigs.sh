#!/bin/bash
# SLURM wrapper for stage 03 (fetch ENCODE bigwigs).
#
# Submit via the yoke REPL with absolute paths:
#
#     sbatch /home/sp5fd/spatial-region-features/genomic-dict/pipeline/slurm/03_fetch_bigwigs.sh
#
# All paths in this script are absolute because yoke forwards sbatch from its
# own working directory (typically brickyard), not from the project root. Using
# absolute paths sidesteps any $SLURM_SUBMIT_DIR confusion.
#
# One-time prep before first submission (sbatch can't create its own log dir):
#
#     mkdir -p /home/sp5fd/spatial-region-features/genomic-dict/logs
#
# Resources: network-bound, CPU-cheap. Two hours is generous; typical
# completion is 15–30 minutes on Rivanna.

#SBATCH --job-name='fetch_bigwigs'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/03_fetch_bigwigs_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/03_fetch_bigwigs_%j.err'
#SBATCH --mem='4GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=02:00:00

set -euo pipefail

WORKSPACE="/home/sp5fd/spatial-region-features"

echo "=== stage 03 fetch_bigwigs ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo "submit dir: $SLURM_SUBMIT_DIR"
echo "workspace:  $WORKSPACE"
echo

source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"

cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/03_fetch_bigwigs.py
