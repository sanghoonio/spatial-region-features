---
date: 2026-04-29
status: draft
description: Pivot the React demo to the BED-file-embedding route (file UMAP as lens designer + on-demand DuckDB NPMI) and strip all heavy precompute parquets except what Section 1 needs. Bigwigs become Section-1-only.
---

# BED-file route + minimum-precompute strip

## Why

Two separate but mutually reinforcing changes:

1. **The interesting question is "what semantics does the embedding learn"** — answered cleanly by exploring file embeddings + the corresponding region partnerships. The current demo treats the region UMAP as primary and named strata as the lens vocabulary. Pivoting to a file-UMAP-driven lens designer makes the question "which biology does THIS file cluster define, and which regions partner under that biology?" — directly addresses what the embedding learns vs what cooccurrence learns.
2. **Precompute is overweight relative to value**. `region_cooccurrence_pmi.parquet` is 97 MB and locks us into 18 named strata. Dropping it in favor of on-demand DuckDB NPMI gives unlimited lens flexibility (any file subset becomes a valid stratum) and slashes ~100 MB of static data.

These align: dropping the cooc parquet is what makes the file-brush lens designer feasible — without it, every brush would need a precomputed stratum, defeating the point.

## Goals

- **Demo surface**: file UMAP becomes the primary interactive surface; region UMAP stays as a complementary view.
- **Lens designer**: brushing the file UMAP defines a custom file pool; the picked region's NPMI partners are computed live against that pool.
- **Precompute strip**: ship only the substrate (region/file embeddings + tokenized corpus + target evidence + UMAPs + Section 1 data). Everything else recomputed on demand or at page load.
- **Bigwig isolation**: bigwig signal is consumed only by the Section 1 continuous-mode plot. Concept axes are recomputed from BED peak counts, not bigwig means.

## What we keep (the substrate)

| Parquet | Purpose | Kept? |
|---|---|---|
| `viz_chr16.parquet` | Region embeddings + UMAP + class. **Without** the per-(cell, mark) bigwig mean columns. | YES (slimmed) |
| `viz_files.parquet` | File embeddings + UMAP + metadata. | YES |
| `tokenized_corpus_chr16.parquet` | The cooccurrence substrate (per-file active token lists). | YES |
| `featured_intervals.parquet` | Section 1 interval metadata. | YES |
| `featured_files.parquet` | Section 1 file curation + per-file chr16_active_token_ids. | YES |
| `featured_tracks.parquet` | Section 1 peaks mode. | YES |
| `featured_signal.parquet` | Section 1 continuous mode (bigwig). | YES |
| `region_target_evidence.parquet` | ENCODE V4 Tier A — external annotation, can't be derived. | YES |

## What we drop

| Parquet | Replacement |
|---|---|
| `region_cooccurrence_pmi.parquet` (97 MB) | On-demand DuckDB NPMI, anchored to picked token + chosen file pool. |
| `region_stratum_marginals.parquet` | On-demand: `SELECT COUNT(*) FROM tokenized_corpus WHERE … AND token_id IN <active>`. |
| `region_concept_axes.parquet` | Recomputed at page load from embeddings + BED-derived case/control sets. |
| `region_class_prototypes.parquet` | Recomputed at page load via `GROUP BY class, AVG(embedding)`. |
| `region_stats.parquet` | On-demand `COUNT(*)` GROUP BY token from tokenized_corpus. |
| `region_target_evidence_summary.parquet` | On-demand aggregate over `region_target_evidence`. |
| `region_modules.parquet` + `module_summary.parquet` | **Defer the module-as-sentence feature**. Re-running Leiden in browser per arbitrary file pool isn't tractable; if we want named-stratum modules later, precompute only those. |

Net storage savings: ~106 MB.

## What changes pipeline-side

### `embedding_features.py` rewrite (BED-only concept axes)

```python
# Before — uses bigwig signal means:
case_mask = max(K562_H3K27ac_mean, K562_H3K4me1_mean, K562_H3K4me3_mean) >= q75
ctrl_mask = K562_H3K27me3_mean >= q75

# After — uses BED-derived peak-call counts:
# Each token's "active in K562 active marks" count = number of K562 active-mark
# files (H3K27ac + H3K4me1 + H3K4me3) where it's peak-called. Top 25% of tokens
# by that count = case set.
case_mask = (n_files_active in K562_active_pool) >= q75
ctrl_mask = (n_files_active in K562_polycomb_pool) >= q75
```

The case/ctrl pools are file lists derived from manifest filters — same patterns as stage 12's stratum filter rules. This makes concept axes derivable from the same substrate (tokenized corpus + file metadata) as everything else, no bigwig dependency.

### `viz_chr16.parquet` schema slim

Drop the per-(cell, mark) signal columns:

```
- K562__ATAC__mean, K562__ATAC__thumb
- K562__H3K4me3__mean, K562__H3K4me3__thumb
- K562__H3K27ac__mean, K562__H3K27ac__thumb
- ...same for GM12878 and HepG2
```

These were used for:
- Concept-axis case/ctrl definitions (replaced above with BED-derived).
- Reference page color toggles like `K562_H3K4me3_mean` (replaced with concept axis projections, which now read off the new BED-derived axes — cell-line specificity already worked this way, just now end-to-end consistent).
- Reference page signal thumbnails (Section 1 alternative tiny-thumb mode in the original Observable page). Not currently scoped in the React demo; can re-add later if needed by querying `featured_signal` directly.

After the slim: `viz_chr16.parquet` shrinks from ~28 MB to ~12 MB.

### Stages affected

- **Stage 08 (`08_precompute_viz.py`)**: drop the bigwig-mean and bigwig-thumb columns from the viz_chr16 emit. This cuts the parquet size and the pipeline runtime (no more bigwig sampling for viz_chr16).
- **`embedding_features.py`**: replace bigwig-signal case/ctrl masks with BED-derived counts. Re-emit `region_concept_axes.parquet`. (For the strip, this parquet itself can be optional — but emitting it keeps the React demo's load fast since axes don't recompute every page load.)
- **Stage 12 (cooc)**: not regenerated, but the resulting parquet **is no longer shipped** to the viz. Pipeline can still produce it for offline analysis. Same with stage 13 (modules), stage 06 summary.

### Files to copy to viz `public/data/dictionary/`

- `viz_chr16.parquet` (slimmed)
- `viz_files.parquet`
- `tokenized_corpus_chr16.parquet`
- `featured_intervals.parquet`
- `featured_files.parquet`
- `featured_tracks.parquet`
- `featured_signal.parquet`
- `region_target_evidence.parquet`
- `region_concept_axes.parquet` (small — keep precomputed for fast page load)

Drop from viz: `region_cooccurrence_pmi.parquet`, `region_stratum_marginals.parquet`, `region_stats.parquet`, `region_target_evidence_summary.parquet`, `region_class_prototypes.parquet`, `region_modules.parquet`, `module_summary.parquet`.

## What changes viz-side

### `MosaicCoordinatorProvider`

- Remove `dict_cooc`, `dict_modules`, `dict_module_summary`, `dict_class_proto`, `dict_region_stats`, `dict_target_summary`, `dict_stratum_marginals` from the table registry and `initializeData()`.
- Keep the `regionsClassed` view extension (the `cclass_category` column for embedding-atlas).
- Add a `filesClassed` view (or just `viz_files` as-is) for the file UMAP.

### `DictCard` rewrite

The kNN section stays as-is (reads from `viz_chr16.knn_token_ids`). The NPMI section is rewritten to compute on demand:

```sql
-- Anchored NPMI: for picked token T against every other token,
-- restricted to files in <pool>. Pool is either a named-stratum file
-- list (for backward compat with the LensPicker) or a brush-derived
-- list of file_ids.
WITH pool AS (SELECT id FROM dict_files WHERE <stratum_filter OR id IN brushed_ids>),
sub AS (
  SELECT id, UNNEST(chr16_active_token_ids) AS token_id
  FROM dict_tokenized_corpus
  WHERE id IN (SELECT id FROM pool)
),
n_total AS (SELECT COUNT(*)::DOUBLE AS N FROM pool),
n_anchor AS (SELECT COUNT(*)::DOUBLE AS n_a FROM sub WHERE token_id = :anchor),
joints AS (
  SELECT s.token_id, COUNT(*)::DOUBLE AS n_ab
  FROM sub s
  WHERE s.id IN (SELECT id FROM sub WHERE token_id = :anchor)
    AND s.token_id != :anchor
  GROUP BY s.token_id
),
marginals AS (SELECT token_id, COUNT(*)::DOUBLE AS n_b FROM sub GROUP BY token_id)
SELECT j.token_id,
       j.n_ab,
       LN(j.n_ab * (SELECT N FROM n_total) / ((SELECT n_a FROM n_anchor) * m.n_b)) AS pmi,
       GREATEST(0, LN(j.n_ab * (SELECT N FROM n_total) / ((SELECT n_a FROM n_anchor) * m.n_b))) /
         (-LN(j.n_ab / (SELECT N FROM n_total))) AS npmi
FROM joints j JOIN marginals m USING (token_id)
ORDER BY npmi DESC LIMIT 30;
```

Cost: anchor-restricted, scales with `O(pool_size × avg_tokens_per_file)`. For our 17K corpus that's ~85M intermediate rows max — DuckDB-WASM handles this in well under a second.

### `LensPicker` becomes a "named stratum" shortcut

Keep the LensPicker for the 18 named biology questions (most common use). It now produces a *file_id list* (translated from the manifest filters at page load) that feeds the on-demand NPMI query as the file pool.

### File UMAP (NEW)

New `<FileUMAP>` component using `EmbeddingViewMosaic` against `viz_files`:

- Categorical color: `assay` or `cell_line` (toggleable later)
- Brush: `vg.intervalXY` selection
- Selected file set drives a "custom lens" overriding the LensPicker

### DictCard "lens" source

Add a small indicator in the DictCard NPMI section showing which lens is active:

- "NPMI · tf_bound_pan_cell" (named lens from picker)
- "NPMI · custom (327 files brushed)" (file UMAP brush)

Brushing clears the named-lens picker (set to "Custom"); picking from the named-lens picker clears the brush. Mutual exclusion via shared state.

### Concept axes computation

At page load, run a small async function that:

1. Loads region embeddings (`viz_chr16` numerical columns or a separate embedding parquet).
2. Loads BED-derived case/ctrl masks via the marginals queries (lifted to JS once).
3. Computes axis vectors via `mean(case) - mean(ctrl)`.
4. Stores them in a context for the DictCard's concept-axis section.

Or — keep a small `region_concept_axes.parquet` if we want the axes precomputed (recommended for fast first paint). Pipeline-side rewrite still applies.

## Migration phases

### Phase A — pipeline cleanup (1 day)

1. Rewrite `embedding_features.py` to use BED-derived case/ctrl. Verify axes still produce sensible UMAP color toggles (eyeball K562_specificity).
2. Drop bigwig mean/thumb columns from stage 08's `viz_chr16` emit.
3. Run on Rivanna: `08_precompute_viz` → `embedding_features`.
4. Pull updated parquets local + sync to `genomic-regions/public/data/dictionary/`.

### Phase B — viz strip (½ day)

1. Remove `dict_cooc`, `dict_modules`, etc. from `lib/duckdb.ts` registry and provider's `initializeData`.
2. Delete the dropped parquets from `public/data/dictionary/`.
3. Verify build still succeeds and Section 1 / RegionUMAP / DictCard still work (the chips will momentarily break — fixed in phase C).

### Phase C — on-demand NPMI rewrite (1 day)

1. Replace `useTokenNpmiPartners` with a hook that builds the on-demand NPMI SQL from the picked token + file pool.
2. Translate named-stratum LensPicker selections into file_id lists at page load (cache).
3. Wire DictCard NPMI section to the new hook.
4. Smoke-test against canvas.md's known cases (FTO → IRX should still surface under tf_bound).

### Phase D — file UMAP + brush lens (1-2 days)

1. New `<FileUMAP>` component (EmbeddingViewMosaic on `viz_files`).
2. Brush selection state + bridge to a "custom file pool" Mosaic Selection.
3. Mutual-exclusion bridge: brush clears LensPicker, named lens clears brush.
4. DictCard lens-source badge.
5. Layout: file UMAP becomes a peer to the region UMAP (side-by-side or stacked).

### Phase E — chr16 strip + concept axes update (½ day)

1. Update `useChr16PartnerPositions` to use the on-demand NPMI hook.
2. Concept axes recomputed at page load OR keep precomputed parquet (decide based on size).

## Risks

1. **On-demand NPMI performance in DuckDB-WASM**. If the anchor-restricted query exceeds ~2s, UX degrades. Mitigation: pre-build a per-(file, token) inverted index in DuckDB at init time. Validate in phase C.
2. **Concept axis quality regression**. BED-derived case/ctrl might separate less cleanly than bigwig-signal-derived (peak calling is a thresholded version of signal). Eyeball comparison in phase A.
3. **File UMAP brush cardinality**. Brushing a huge file set (>10K files) might make the cooccurrence query slow. Mitigation: cap the pool size with a soft warning, or downsample.
4. **Loss of module catalogue / soft class profiles**. Currently not in the React demo; deferred. Mark as acceptable since the demo's primary goal shifts to the file-UMAP exploration.

## Effort estimate

**~4 working days** for phases A–E. Pipeline change (A) is the biggest bottleneck — one Rivanna run + sync. Viz changes are mostly mechanical once the data shape is right.

## Plan comprehension quiz topics

When implementation starts:

1. Why drop `region_cooccurrence_pmi.parquet` instead of keeping it as a fast cache for the named strata? *Expected answer*: the named strata can still hit fast on-demand queries (anchor-restricted), and keeping the parquet means we can't cleanly extend to brush-defined custom lenses without also recomputing on demand for those — so we may as well unify on the on-demand path and gain ~100 MB.
2. Why move concept axes from bigwig signal to BED peak counts? *Expected answer*: methodological cleanliness — the embedding was trained on BED cooccurrence, so axes derived from BED are interpretable in the embedding's own data substrate. Bigwig signal mixes a different (continuous) data source into a model that only saw thresholded peaks.
3. What's the riskiest part of this pivot? *Expected answer*: on-demand NPMI performance under brush lenses. The anchor-restricted query is fine for ~17K files; brushing a much smaller subset is faster, but a degenerate "brush almost everything" case could be slow. Validate in phase C with a couple of canvas.md test cases.
