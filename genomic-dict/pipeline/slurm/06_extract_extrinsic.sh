#!/bin/bash
# Stage 06 — extrinsic annotations. ~1.5 hr (GTEx parsing dominates).
#SBATCH --job-name='extract_extrinsic'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/06_extract_extrinsic_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/06_extract_extrinsic_%j.err'
#SBATCH --mem='16GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=03:00:00

set -euo pipefail
WORKSPACE="/home/sp5fd/spatial-region-features"
echo "=== stage 06 extract_extrinsic ==="
echo "host: $(hostname); date: $(date)"
source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"
cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/06_extract_extrinsic.py
