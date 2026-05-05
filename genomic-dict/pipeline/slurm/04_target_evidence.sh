#!/bin/bash
# Stage 04 — ENCODE V4 cCRE-Gene Links → target evidence. Coord-based match
# from viz_chr16 to V4 EH38E accessions, joins 3D-Chromatin + eQTL evidence
# rows. ~5-10 min depending on filesystem; modest memory.
#SBATCH --job-name='target_evidence'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/04_target_evidence_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/04_target_evidence_%j.err'
#SBATCH --mem='8GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=00:30:00

set -euo pipefail
WORKSPACE="/home/sp5fd/spatial-region-features"
echo "=== stage 04 target_evidence ==="
echo "host: $(hostname); date: $(date)"
source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"
cd "$WORKSPACE"
uv run python genomic-dict/pipeline/scripts/04_target_evidence.py
