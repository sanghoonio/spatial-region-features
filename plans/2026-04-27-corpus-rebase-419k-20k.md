---
date: 2026-04-27
status: complete
description: Rebase the pipeline corpus onto the full 419k UMAP parquet, balanced-cap to ~17k, recompute file embeddings + UMAP from R2V tokens
---

> **Resolved during planning (2026-04-27):**
> - Strategy = **balanced + 2-pass reserve** (16,799 files). Proportional 20k misses 10/18 stage-10 cells; balanced misses only 1 (which is unavoidable — see below).
> - Embedding = `model.encode(rs, pooling="mean")`. Confirmed via geniml source that R2V is Word2Vec-derived (no Doc2Vec inference); `encode(pooling="mean")` is mathematically identical to mean-of-token-embeddings, just better wired since the model is already loaded for the tokenizer.
> - `(GM12878, H3K9me3)` has zero files in the 84,698 source pool. Drop that cell from stage 10's featured grid (asymmetric: GM12878 gets 5 marks, K562 and HepG2 get all 6).
> - Reserve = 5 files per (cell_line, target) cell, ranked by quality score (count of populated rich-meta fields).

# Corpus rebase: 419k → balanced 17k, recompute embeddings + UMAP

## Motivation

Current pipeline has two parallel corpora that drift apart:
- **84k Tight** (`manifest.parquet`): filter from `hg38_umap_all_3_13.parquet` — 5 assays + cell_line known. Used by stage 10 (featured narrative).
- **12.8k rich** (`rich_corpus_manifest.parquet`): filter from t1+t2 rich metadata — adds species + ChIP-target requirements. Used by stages 09/11.

This split causes (a) the file UMAP, the network, and the featured-file activation overlay to draw from different pools, and (b) stage 09 had to fall back to `hg38_geometry.parquet` (49k coords-only) for a denser layout, which is itself just a precomputed UMAP on a different curated subset — not a substantive metadata source.

**Goal:** one corpus, one set of file embeddings, one UMAP, all computed by us from the same R2V model that the chr16-region viz uses.

## Decisions

1. **Source = 419k `hg38_umap_all_3_13.parquet`.** Drop `hg38_geometry.parquet` entirely (it has only `id, x, y` — no metadata).
2. **Filter:** `assay ∈ {DNase-seq, ATAC-seq, ChIP-seq, TF ChIP-seq, Histone ChIP-seq}` AND `cell_line` non-empty AND `cell_line != "UNKNOWN"`. Yields 84,698 rows.
3. **Cap:** balanced per-assay — keep all 799 Histone ChIP-seq + 4,000 each of the other four → ~17k total. Random sample within each assay (seeded). This trades ~3k corpus size for real cross-assay representation; histone ChIP would otherwise be starved (190 files in proportional 20k sample).
4. **Rich-metadata left-join.** Keep the `t1` + `t2` left-join onto survivors for `target`, `antibody`, `treatment`, `global_experiment_id`. Stage 10's featured-file selector needs `target`; the 46,651 overlap with the rich subset will cover all (K562, GM12878, HepG2) × histone-mark cells.
5. **File embeddings = mean of per-token embeddings.** Looked up from `pretrained_universe.parquet` (which already has 100-dim embedding per token, via stage 07). No `model.encode()` inference needed — saves ~3× wall time on Rivanna and matches the embedding space of the chr16-region viz exactly.
6. **UMAP recomputed locally** from those embeddings using `umap-learn` with `n_neighbors=15, min_dist=0.1, seed=42` (same params as stage 08).

## Stage changes

### Stage 01 (rewrite, light)
- Already filters from `hg38_umap_all_3_13.parquet`. Add the balanced per-assay cap. Keep the rich-metadata join.
- Output: `data/corpus/manifest.parquet` (~17k rows, columns: `id, name, description, assay, cell_line, cell_type, tissue, target, antibody, treatment, global_experiment_id, global_sample_id`).

### Stage 01b (delete)
- Redundant — its filter is now subsumed by stage 01.
- Remove the script and the config block.

### Stage 09 (rewrite)
- Read `manifest.parquet` + `data/precomputed/file_embeddings.parquet` (new — output of stage 11).
- Compute UMAP locally on the 100-dim embeddings; emit `viz_files.parquet` with `id, name, description, umap_x, umap_y, assay, cell_line, cell_type, tissue, is_unlabeled=false`.
- Drop `hg38_geometry.parquet` from inputs.

### Stage 11 (extend)
- Tokenize the ~17k whole-genome (not just chr16-active filtering, but the full token list per file).
- For each file, compute mean embedding across all non-UNK tokens (lookup from `pretrained_universe.parquet`).
- Emit two parquets:
  - `tokenized_corpus_chr16.parquet` — id + chr16 token list (existing schema, for stages 12/13 region work)
  - `file_embeddings.parquet` — id + 100-dim mean embedding (new, consumed by stage 09)
- Streaming write, batch_size=5000.

### Stage 10 (no code change)
- Reads `manifest.parquet` which is now ~17k. Featured-file selection still works (cell_line × target hits are common). `viz_files.parquet` no longer has `is_unlabeled=true` rows, so the mystery-files branch yields 0 — already the case.

### Config edits
- Remove `01b_curate_rich_corpus` block.
- Stage 01: add `per_assay_cap` (4000) and `histone_chip_assay` marker (since histone has fewer than the cap).
- Stage 09: drop `source_parquet`, `manifest_relative` keys; add `umap_n_neighbors`, `umap_min_dist`, `umap_seed`.
- Stage 11: add `embeddings_output_relative` and drop the chr16-only restriction from the tokenizer pass.

## Execution order

1. Edit stage 01 + config; run locally → produces ~17k manifest.
2. Edit stage 11; run on Rivanna → produces both parquets. Estimated ~5 min for tokenization on 17k.
3. Edit stage 09; run locally → produces UMAP + viz_files.parquet.
4. Re-run stage 10 → produces featured-narrative parquets against new manifest.
5. Verify mattress / Observable viz still loads with the new viz_files.parquet.

## Risks / things to check

- **(K562, GM12878, HepG2) × histone marks coverage in the cap.** Stage 01's random-within-assay sampling could in theory miss specific (cell_line, target) cells we need for stage 10. Mitigation: stage 01 should *first* set aside any file matching a featured (cell_line, target) cell, then sample the rest within each assay. (Need to encode this constraint.)
- **Mean-of-token embeddings ≠ `model.encode()`.** For visualization this is fine — the gross structure of the file UMAP will be similar — but if we later use these embeddings as features for downstream models, we should note the choice in provenance.
- **Stage 08 cooccurrence dependency.** Stage 08 may have a corpus-dependent cooccurrence pass that becomes stale. Check whether anything still consumes its output before declaring stage 13 fully replaces it.

## Open questions for follow-up sessions

- Stage 13 stratification axes — `tissue` is mostly empty in the 419k; `cell_line` is the natural axis. ChIP target stratification (CTCF / H3K27me3 / etc.) would be richer but requires the t1+t2 join survives the cap.
- Whether to also recompute `featured_intervals.parquet` chr16-token slicing if stage 11's universe coverage changes.

## Implementation log (2026-04-27)

**Done locally:**
- `genomic-dict/pipeline/scripts/01_curate_corpus.py` — rewritten with quality score, two-pass reserve, balanced per-assay cap. Verified locally: 419,756 → 84,698 → 16,799, all 17 stage-10 cells covered (GM12878 H3K9me3 silently skipped via `skip_cells`).
- `genomic-dict/pipeline/scripts/01b_curate_rich_corpus.py` — deleted.
- `genomic-dict/pipeline/scripts/11_tokenize_corpus_chr16.py` — extended to also emit `file_embeddings.parquet` with per-file mean embedding (lookup from `pretrained_universe.parquet`'s embedding column).
- `genomic-dict/pipeline/scripts/09_prepare_file_viz.py` — rewritten to compute UMAP locally (`umap-learn`, cosine metric, n_neighbors=15, min_dist=0.1, seed=42) from stage 11's embeddings. No more dependency on `hg38_geometry.parquet`.
- `genomic-dict/pipeline/scripts/10_featured_narrative.py` — `select_featured_files()` now honors `skip_cells`.
- `genomic-dict/config.yaml` — stage 01 reshaped (added `per_assay_cap=4000`, `reserve_per_featured_cell=5`); stage 09 reshaped (UMAP params, drop geometry input); stage 11 added `embedding_pooling`; stage 10 added `skip_cells: [["GM12878", "H3K9me3"]]`; stage 01b block removed.
- `genomic-dict/pipeline/slurm/11_tokenize_corpus_chr16.sh` — comment updated (~30-45 min vs ~3 hr).

**Pending — needs Rivanna:**
- ~~Run stage 11 on Rivanna against the 17k manifest.~~ Done — job 12251912, 99 min, 16,794 emb + 16,755 chr16 lists, 0 failures.
- ~~Sync embeddings parquet locally → run stage 09 locally → run stage 10 locally.~~ Done — stage 09 produced 16,794-row viz_files.parquet; stage 10 (job 12262747, 27 s) produced 17 featured files + 40 track records.

**Stale output files (not deleted, left for user to git-clean):**
- `data/corpus/rich_corpus_manifest.parquet`
- `data/corpus/umap_corpus_manifest.parquet`
- `results/01b_curate_rich_corpus/`

**Sanity check at handoff:** stage 01 summary at `results/01_curate_corpus/summary.json` shows status=success, n_final=16,799, all 17 cells covered, mean quality_score=1.70, target coverage 28.1%.
