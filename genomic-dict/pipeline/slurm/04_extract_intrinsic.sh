#!/bin/bash
# SLURM wrapper for stage 04 (bigwig-based intrinsic features).
#
# Submit via yoke REPL with absolute path:
#
#     sbatch /home/sp5fd/spatial-region-features/genomic-dict/pipeline/slurm/04_extract_intrinsic.sh
#
# Reads the 18 bigwigs from stage 03 and produces
# data/annotations/intrinsic_bigwig.parquet.
#
# Compute: ~1.1M per-interval pyBigWig `stats` calls. Typical runtime
# 5–15 minutes on one CPU. One hour is a generous cap.

#SBATCH --job-name='extract_intrinsic'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/04_extract_intrinsic_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/04_extract_intrinsic_%j.err'
#SBATCH --mem='8GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=01:00:00

set -euo pipefail

WORKSPACE="/home/sp5fd/spatial-region-features"

echo "=== stage 04 extract_intrinsic ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo

source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"

cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/04_extract_intrinsic.py
