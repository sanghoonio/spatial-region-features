#!/bin/bash
# SLURM wrapper for stage 09 (file-level viz data prep).
#
# Submit via yoke REPL with absolute path:
#
#     sbatch /home/sp5fd/spatial-region-features/genomic-dict/pipeline/slurm/09_prepare_file_viz.sh
#
# Reads stage-00 HF parquet + stage-01 manifest; produces viz_files.parquet
# in ~10 seconds. Trivially small but needs to run on Rivanna so downstream
# stage 10 can find the parquet.

#SBATCH --job-name='prepare_file_viz'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/09_prepare_file_viz_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/09_prepare_file_viz_%j.err'
#SBATCH --mem='4GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=00:15:00

set -euo pipefail

WORKSPACE="/home/sp5fd/spatial-region-features"

echo "=== stage 09 prepare_file_viz ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo

source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"

cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/09_prepare_file_viz.py
