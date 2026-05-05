#!/usr/bin/env bash
# pull-results.sh — Pull pipeline summaries and model checkpoints from rivanna.
# Usage: ./pull-results.sh [--dry-run]
set -euo pipefail

REMOTE="sp5fd@login.hpc.virginia.edu:/home/sp5fd/spatial-region-features"
HERE="$(cd "$(dirname "$0")" && pwd)"

# Pull stage summaries + any rendered HTML from genomic-dict/results/**
# || true handles first-run when remote dir doesn't exist yet.
mkdir -p "$HERE/genomic-dict/results"
rsync -avz "$@" \
  --include='*/' \
  --include='*.json' \
  --include='*.md' \
  --include='*.yaml' \
  --include='*.html' \
  --exclude='*' \
  "$REMOTE/genomic-dict/results/" "$HERE/genomic-dict/results/" || true

# Pull SLURM logs (small text files; useful for debugging).
mkdir -p "$HERE/genomic-dict/logs"
rsync -avz "$@" \
  --include='*.out' \
  --include='*.err' \
  --exclude='*' \
  "$REMOTE/genomic-dict/logs/" "$HERE/genomic-dict/logs/" || true

# Pull trained R2V model checkpoints (can be large — only after training runs).
mkdir -p "$HERE/genomic-dict/models"
rsync -avz "$@" \
  "$REMOTE/genomic-dict/models/" "$HERE/genomic-dict/models/" || true
