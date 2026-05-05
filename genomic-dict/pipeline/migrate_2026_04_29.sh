#!/usr/bin/env bash
# One-time migration after the 2026-04-29 pipeline cleanup.
# Idempotent — safe to re-run.
#
# Run with REMOTE=1 to apply on Rivanna via ssh; otherwise applies locally.
#
#   ./genomic-dict/pipeline/migrate_2026_04_29.sh           # local Mac
#   REMOTE=1 ./genomic-dict/pipeline/migrate_2026_04_29.sh  # Rivanna
set -euo pipefail

remote_cmd() {
  ssh sp5fd@login.hpc.virginia.edu "cd /home/sp5fd/spatial-region-features/genomic-dict && $*"
}

run() {
  if [[ "${REMOTE:-0}" == "1" ]]; then
    remote_cmd "$@"
  else
    cd "$(dirname "$0")/../.." # genomic-dict/
    bash -c "$*"
  fi
}

OLD_DIRS=(
  04_extract_intrinsic
  05_extract_sequence
  06_target_evidence
  07_load_pretrained
  08_precompute_viz
  09_prepare_file_viz
  10_featured_narrative
  11_tokenize_corpus_chr16
  14_featured_signal
  embedding_features
)

archive_dirs="${OLD_DIRS[*]}"

run "
  echo '== migrate_2026_04_29 =='
  if [ -f data/annotations/tokenized_corpus_chr16.parquet ]; then
    mv data/annotations/tokenized_corpus_chr16.parquet data/precomputed/ &&
      echo '  moved tokenized_corpus_chr16.parquet → precomputed/'
  fi
  if [ -f data/annotations/intrinsic_bigwig.parquet ]; then
    rm data/annotations/intrinsic_bigwig.parquet &&
      echo '  removed intrinsic_bigwig.parquet'
  fi
  mkdir -p results/_archived_2026_04_29
  moved=0
  for d in $archive_dirs; do
    if [ -d \"results/\$d\" ]; then
      mv \"results/\$d\" results/_archived_2026_04_29/
      moved=\$((moved+1))
    fi
  done
  echo \"  archived \$moved old results dirs\"
  echo '  remaining results/:'
  ls results/ | sed 's|^|    |'
"
