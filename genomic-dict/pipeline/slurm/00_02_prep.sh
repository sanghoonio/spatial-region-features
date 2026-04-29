#!/bin/bash
# SLURM wrapper — run stages 00, 01, 02 in sequence.
#
# Submit via yoke REPL with absolute path:
#
#     sbatch /home/sp5fd/spatial-region-features/genomic-dict/pipeline/slurm/00_02_prep.sh
#
# Stage 00 — fetch HF metadata parquets (tiny).
# Stage 01 — curate corpus manifest from stage 00 output.
# Stage 02 — download + filter SCREEN V4 cCRE universe.
#
# All three stages are fast and network-bound (no compute). Total runtime
# expected: 1–3 minutes. 30-minute time limit is generous.
#
# One-time prep before first submission:
#     mkdir -p /home/sp5fd/spatial-region-features/genomic-dict/logs

#SBATCH --job-name='prep_00_02'
#SBATCH --output='/home/sp5fd/spatial-region-features/genomic-dict/logs/00_02_prep_%j.out'
#SBATCH --error='/home/sp5fd/spatial-region-features/genomic-dict/logs/00_02_prep_%j.err'
#SBATCH --mem='4GB'
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --partition='standard'
#SBATCH --account='shefflab'
#SBATCH --time=00:30:00

set -euo pipefail

WORKSPACE="/home/sp5fd/spatial-region-features"

echo "=== prep 00_02 ==="
echo "host: $(hostname)"
echo "date: $(date)"
echo

source /project/shefflab/rivanna_config/env.sh
export PATH="$HOME/.local/bin:$PATH"

cd "$WORKSPACE"

echo "--- stage 00: inspect metadata ---"
uv run python genomic-dict/pipeline/scripts/00_inspect_metadata.py
echo

echo "--- stage 01: curate corpus ---"
uv run python genomic-dict/pipeline/scripts/01_curate_corpus.py
echo

echo "--- stage 02: prepare universe ---"
uv run python genomic-dict/pipeline/scripts/02_prepare_universe.py
echo

echo "=== prep complete ==="
