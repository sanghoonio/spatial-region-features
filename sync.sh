#!/usr/bin/env bash
# sync.sh — Deploy spatial-region-features workspace to Rivanna via rsync over SSH.
# Pulls results first so newer remote outputs don't get clobbered.
# Usage: ./sync.sh [--dry-run] [other rsync flags]
set -euo pipefail

REMOTE="sp5fd@login.hpc.virginia.edu:/home/sp5fd/spatial-region-features"
HERE="$(cd "$(dirname "$0")" && pwd)"

echo "=== Pulling results from rivanna ==="
"$HERE/pull-results.sh" "$@"

echo
echo "=== Pushing code to rivanna ==="
rsync -avzu "$@" \
  --exclude='.venv/' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  --exclude='.pytest_cache/' \
  --exclude='.ruff_cache/' \
  --exclude='.mypy_cache/' \
  --exclude='.ipynb_checkpoints/' \
  --exclude='.DS_Store' \
  --exclude='genomic-dict/data/' \
  --exclude='genomic-dict/logs/' \
  --exclude='genomic-dict/models/*.pt' \
  --exclude='genomic-dict/models/*.joblib' \
  --exclude='genomic-dict/models/*.bin' \
  --exclude='genomic-dict/models/*.safetensors' \
  --exclude='scratch/chr22_demo/data/' \
  --exclude='scratch/ccre_examples/data/' \
  "$HERE/" "$REMOTE/"
