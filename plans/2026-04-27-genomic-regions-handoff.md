---
date: 2026-04-27
status: complete
description: Handoff doc for the genomic-regions Observable Framework session — bundle inventory, schemas, decisions locked in, and open viz design choices
---

# Handoff to genomic-regions session

This document is a self-contained brief for a fresh Claude Code session working in the **genomic-regions** repo (https://github.com/sanghoonio/genomic-regions.git, Observable Framework migration target). The data bundle is finalized; this doc explains what's in it, what's been decided, and what choices remain for the viz layer.

Source pipeline: `genomic-dict/pipeline/scripts/` in `spatial-region-features`. Prior plan: `plans/2026-04-27-corpus-rebase-419k-20k.md` (status: complete).

---

## TL;DR

- 8 parquet files, ~117 MB total, sitting in `mattress/public/data/dictionary/` (copy these to genomic-regions' `src/data/`).
- Corpus is **16,799 BED files** filtered from the bedbase 419k UMAP parquet, balanced 4k/assay + 5-per-cell reserve for the featured 17.
- File-level UMAP is **genome-wide** (mean of all-chrom token embeddings → UMAP). Region-level UMAP is **chr16-only**. Don't conflate these.
- Cooccurrence network is **corpus-derived (Jaccard)**, not embedding-kNN. Stratified by `all` / `K562` / `GM12878` / `HepG2`. The contrast across contexts is the Atacformer-motivation evidence layer.
- The featured grid is **asymmetric**: K562 and HepG2 each get ATAC + 5 histone marks (6 cells); GM12878 gets only 5 (no H3K9me3 — zero files in the source pool, can't fix by sampling).

---

## Bundle inventory

All paths relative to the data dir.

### `viz_files.parquet` — 580 KB, 16,794 rows

One row per BED file in the corpus. Drives the file UMAP (Step 4 dot field).

| col | type | notes |
|---|---|---|
| `id` | String | md5 hash; primary key. Joins to `featured_files.file_id`, `tokenized_corpus_chr16.id`, manifest |
| `name` | String | original ENCODE / lab filename |
| `description` | String | `""` if missing (filled, never null) |
| `umap_x`, `umap_y` | Float32 | genome-wide R2V mean → UMAP (cosine, n_neighbors=15, min_dist=0.1, seed=42) |
| `assay` | String | one of: `DNase-seq`, `ATAC-seq`, `ChIP-seq`, `TF ChIP-seq`, `Histone ChIP-seq` |
| `cell_line` | String | non-null, never `"UNKNOWN"` (filtered out at stage 01) |
| `cell_type` | String | `""` if missing |
| `tissue` | String | `""` if missing |
| `is_unlabeled` | Boolean | always `false` in the new corpus (legacy column kept for schema stability) |

Common queries: `SELECT * FROM viz_files WHERE id = ?` (lookup); `SELECT * FROM viz_files` (full UMAP draw — 16k dots is fine for vgplot).

### `viz_chr16.parquet` — 28 MB, 35,934 rows

One row per chr16 universe token. Drives the region UMAP (Step 4 region panel) and per-region tooltip data.

Columns include `token_id`, `region` (e.g., `chr16:10001-10217`), `chrom`, `start`, `end`, `cclass` (SCREEN class: `PLS`, `pELS`, `dELS`, `CA-CTCF`, `CA-H3K4me3`, `unclassed`), `overlaps_screen`, three tissues × `{ATAC__mean, ATAC__thumb (list), H3K4me1__mean, H3K4me3__mean, H3K27ac__mean, H3K27me3__mean, H3K9me3__mean}`, `umap_x`, `umap_y`, `knn_token_ids` (list[i64]), `knn_distances` (list[f64]).

⚠️ **`knn_*` is kNN in 100-dim embedding space**, NOT corpus cooccurrence. Use it for "which tokens *look* similar to R2V" (Step 3 narrative). Use `region_cooccurrence.parquet` for "which tokens *actually appear together*" (Step 5 / network panel).

This parquet is the output of stage 08, which has not been rerun since the corpus rebase. The bigwig signal columns are unchanged (they only depend on the chr16 cCREs, not the corpus). The kNN edges are also unchanged.

### `region_stats.parquet` — 736 KB, 35,934 rows

One row per chr16 token, per-token activation stats across the 16,799-file corpus.

| col | type | meaning |
|---|---|---|
| `token_id` | Int64 | joins `viz_chr16.token_id` |
| `n_files_total` | UInt32 | files where this token is active |
| `freq_corpus` | Float32 | `n_files_total / 16799` |
| `n_cell_lines_active` | UInt32 | distinct cell_lines where the token activates ≥1 time (pleiotropy proxy) |
| `n_assays_active` | UInt32 | distinct assays where active |
| `n_K562`, `n_GM12878`, `n_HepG2` | Int64 | per-cell-line activation counts |
| `freq_K562`, `freq_GM12878`, `freq_HepG2` | Float32 | normalized by total files of that cell_line in the corpus |
| `n_DNase`, `n_ATAC`, `n_ChIP`, `n_TF_ChIP`, `n_Histone_ChIP` | Int64 | per-assay activation counts |

Use this for: per-region tooltip, sizing/coloring nodes in the network, the "this region activates in N% of files" callout for Step 5.

### `region_cooccurrence.parquet` — 26 MB, 143,736 rows

For each (token, context) pair, the top-30 partner tokens ranked by Jaccard similarity in that context.

| col | type | meaning |
|---|---|---|
| `token_id` | Int64 | the focal region |
| `context` | String | `"all"` / `"K562"` / `"GM12878"` / `"HepG2"` |
| `n_files_active` | Int64 | files (within this context) where `token_id` is active |
| `partner_token_ids` | List[Int64] | top-30 partners, sorted by Jaccard desc |
| `weights_jaccard` | List[Float64] | Jaccard similarity, parallel to `partner_token_ids` |
| `counts` | List[Int64] | raw cooccurrence count, parallel |

35,934 tokens × 4 contexts = 143,736 rows. Each list is up to 30 elements (shorter if the token has fewer than 30 nonzero partners in that context).

Jaccard formula: `J(a,b) = c(a,b) / (n_a + n_b − c(a,b))` where `c(a,b)` = files where both active, `n_a`, `n_b` = files where each is active. Symmetric: `(a → b)`'s Jaccard equals `(b → a)`'s.

Use this for: the cooccurrence network panel; the per-context comparison panel (Atacformer evidence).

### `tokenized_corpus_chr16.parquet` — 62 MB, 16,755 rows

Per-file chr16 token activation lists. Powers the click-on-any-file → highlight-on-region-UMAP path.

| col | type | meaning |
|---|---|---|
| `id` | String | joins `viz_files.id` |
| `chr16_active_token_ids` | List[Int64] | distinct chr16 tokens active in this file (sorted) |
| `n_chr16_active` | Int64 | length of the list |

16,755 < 16,794 — 39 files have a non-trivial mean embedding (so they're in `viz_files`) but no chr16-active tokens (so they're not in this parquet). The viz should treat missing rows as "no chr16 activation."

⚠️ The chr16 token list per file averages ~3,000 tokens; some files reach 13,000+. Don't hold all 16k lists in memory at once — query keyed on `id`.

### `featured_files.parquet` — 108 KB, 17 rows

The 17 hand-picked files for the narrative arc (3 cell_lines × 6 marks − GM12878/H3K9me3).

| col | type | meaning |
|---|---|---|
| `file_id` | String | joins `viz_files.id` |
| `name`, `assay`, `cell_line`, `target` | String | display metadata |
| `role` | String | `"featured"` (mystery role currently unused) |
| `n_chr16_active_tokens` | Int64 | duplicates `tokenized_corpus_chr16.n_chr16_active` |
| `chr16_active_token_ids` | List[Int64] | duplicates the same column from the tokenized corpus |

Strictly redundant with `tokenized_corpus_chr16` for the 17 featured files, but kept for fast small-bundle lookups (DuckDB-WASM SELECT on a 17-row parquet is instant; on a 65 MB parquet it takes a moment).

### `featured_intervals.parquet` — 8 KB, 4 rows

The 4 chr16 narrative hotspots:

| `interval_id` | label | bp | what it shows |
|---|---|---|---|
| `alpha_globin_genes` | α-globin gene cluster | chr16:218,000–238,000 | active-promoter signature |
| `ciita_promoters` | CIITA alternative promoters | chr16:10,960,000–10,980,000 | cell-type-specific promoter use |
| `cdh1_promoter` | CDH1 / E-cadherin promoter | chr16:68,725,000–68,745,000 | epithelial active promoter |
| `fto_enhancer` | FTO intron 1 enhancer | chr16:53,760,000–53,780,000 | the "nearest gene is wrong" case (acts on IRX3 ~1 Mb away) |

Plus `narrative_caption` (Markdown-friendly text), `universe_token_ids` (list of chr16 tokens overlapping this interval, ~14–23 each).

### `featured_tracks.parquet` — 8 KB, 40 rows

For each (`file_id`, `interval_id`) combo where there's signal: chr16 peak coords + active universe tokens within the interval. Drives Step 1 (per-file peaks) and Step 2 (per-file token activation overlay) in mattress's narrative.

| col | meaning |
|---|---|
| `file_id` | foreign key to `featured_files.file_id` |
| `interval_id` | foreign key to `featured_intervals.interval_id` |
| `n_active_tokens`, `active_token_ids` | tokens *within this interval* that are active in this file |
| `n_peaks`, `peak_starts`, `peak_ends` | original BED peak coords (chr16-only, within interval) |

40 = 17 files × 4 intervals minus pairs with no signal in that window.

---

## Mapping bundle → 5-step narrative

| Step | What it shows | Parquets needed |
|---|---|---|
| 1 — peaks per file | Featured BEDs as raw peaks at the 4 hotspots | `featured_intervals`, `featured_files`, `featured_tracks` |
| 2 — universe + tokenization | Universe band + per-file activation | `featured_intervals`, `featured_tracks` |
| 3 — embedding emerges | Region UMAP colored by SCREEN class; kNN connections in embedding space | `viz_chr16` (`umap_x/y`, `cclass`, `knn_token_ids`) |
| 4 — file UMAP cross-link | Click file dot → highlight its tokens on region UMAP | `viz_files`, `viz_chr16`, `tokenized_corpus_chr16` (or `featured_files` if click is restricted to featured 17) |
| 5 — hypothesis / network | Region's cooccurrence partners; per-context contrast | `region_stats`, `region_cooccurrence`, `viz_chr16` (for partner coords + bigwig overlays) |

---

## Decisions locked in (with rationale)

These shaped the bundle. Document them so the genomic-regions session can revisit if needed; don't silently change without re-running upstream stages.

1. **Corpus = balanced 4k/assay from 419k**, not the natural 84k Tight or the 12k rich subset. Reason: file UMAP needs cross-assay representation; histone ChIP only has 799 in the source so it's the binding constraint. Random within each assay, but reserve top-5-by-quality per (cell_line × target) cell first so the featured grid survives the cut.

2. **Featured grid asymmetric**: GM12878 has no H3K9me3 file in the entire 84,698-row pool. The viz must handle a 5-cell row for GM12878 vs 6-cell rows for K562 / HepG2. Don't render an empty cell — skip it. `featured_files.parquet` already drops it (17 rows = 3·6 − 1).

3. **File embedding = mean of all-chrom non-UNK token embeddings** (model.encode(rs, pooling="mean") and pretrained_universe parquet lookup are mathematically identical because R2V is Word2Vec-derived; no Doc2Vec inference). Genome-wide, not chr16. The chr16 restriction is *only* for the per-file token list output and the region UMAP, not the file UMAP.

4. **Region UMAP and region cooccurrence are different views**: kNN in `viz_chr16.knn_*` is "embedding-similar"; partners in `region_cooccurrence` are "actually co-active in the corpus." Both are meaningful; they're not interchangeable.

5. **Context stratification = `all` + 3 cell_lines**, not assay or tissue. Tissue field is mostly empty in the 419k. Assay-stratified cooccurrence is feasible but wasn't computed; the cell_line contrast is the cleanest Atacformer-motivation argument.

6. **Jaccard chosen over lift / PMI**. Bounded [0, 1], interpretable, less misleading at low counts than lift. Raw `counts` are also stored if a different metric is needed.

7. **Stage 08's `viz_chr16.parquet` was NOT rerun** with the new corpus. Bigwig signals and kNN edges are corpus-independent so this is fine; if the new session wants corpus-derived metrics on this parquet, that's an additional stage 08 run.

---

## Known gaps

- **Stage 05 outputs** (FIMO motifs, ChromHMM states, phastCons) — not in the bundle. Pipeline scripts are written but not all run. If Step 3 narrative needs "this region has a CTCF motif" or "this region is in K562's H3K27ac chromHMM state," wire stage 05 first.
- **Stage 06 outputs** (nearest GENCODE gene, GWAS Catalog overlap, GTEx eQTL) — not in the bundle. Same situation.
- **All-chrom region UMAP** — not computed. Stage 08 only ran for chr16. Doing genome-wide would 30× the bundle (~800 MB) so it needs deliberate decision before running.
- **viz_chr16's bigwig signals are still K562/GM12878/HepG2 + 5 histones + ATAC**, the same lineup as the corpus rebase. No change needed unless the narrative expands beyond these.
- **Per-assay cooccurrence** — only the 4 cell-line-stratified contexts exist. Per-assay would be a natural addition; rerun stage 13 with `contexts: [DNase-seq, ATAC-seq, ChIP-seq, ...]`.

---

## Open viz design choices (genomic-regions session decides)

These were not pre-resolved because they're rendering choices, not data choices.

1. **Cooccurrence network rendering**: force-directed graph (familiar but expensive at 35k nodes; would need to subsample to top-N by `freq_corpus` or `n_cell_lines_active`), adjacency matrix heatmap (scales better, less narrative), or "ego network" (center on clicked region, draw its top-30 partners as a radial layout). Recommended: ego network — pairs naturally with the click-a-region UX of Step 5 and keeps the visual budget small.

2. **Per-context contrast UI**: side-by-side small multiples (4 ego networks for the same focal region across `all` / `K562` / `GM12878` / `HepG2`), or a "diff" view that highlights edges that change rank between contexts. The latter is closer to the Atacformer argument but harder to read.

3. **Click-any-file vs click-featured-only**: the 62 MB `tokenized_corpus_chr16.parquet` enables any-file click-through; the 108 KB `featured_files.parquet` is faster but limited to 17. Mattress currently uses the featured path. Recommended for genomic-regions: any-file via DuckDB-WASM lazy query; first-load cost ~5s but then queries are sub-second and the UX is much richer.

4. **Featured grid rendering**: handle the GM12878 H3K9me3 gap explicitly. Options: render a `—` placeholder with a hover note ("not available in source corpus"), or just lay out the 17 cells in a 6+6+5 row pattern. Either way, don't drop GM12878 from the row entirely.

5. **Step 3 kNN visualization** is the easiest place to confuse "embedding similarity" with "corpus cooccurrence." Add a sentence in the panel narrative; don't let users assume `viz_chr16.knn_*` and `region_cooccurrence.partners` are the same thing.

---

## Loading patterns (DuckDB-WASM)

DuckDB-WASM is the assumed runtime. List columns require unnesting for some operations.

```sql
-- Step 4: lookup a clicked file's chr16 activation
SELECT chr16_active_token_ids FROM 'tokenized_corpus_chr16.parquet' WHERE id = ?

-- Step 5: ego network for a clicked region in K562 context
SELECT partner_token_ids, weights_jaccard, counts
FROM 'region_cooccurrence.parquet'
WHERE token_id = ? AND context = 'K562'

-- Region-level summary card (join cooccurrence row with stats + universe coords)
SELECT s.*, u.region, u.cclass, u.umap_x, u.umap_y
FROM 'region_stats.parquet' s
JOIN 'viz_chr16.parquet' u USING (token_id)
WHERE s.token_id = ?

-- Unnest cooccurrence partners for a join (Observable / Plot rendering)
SELECT token_id, context, p.partner_id, p.weight, p.count
FROM 'region_cooccurrence.parquet',
     UNNEST(partner_token_ids, weights_jaccard, counts) AS p(partner_id, weight, count)
WHERE token_id = ? AND context = 'all'
```

DuckDB-WASM cold-start on the 62 MB tokenized parquet is ~5s; warm queries with `WHERE id = ?` are <100 ms. Use `httpfs` to fetch parquets lazily; don't `COPY INTO` upfront.

---

## Pipeline state (reference)

| Stage | Status | Output (in bundle?) |
|---|---|---|
| 00 inspect | done | n/a |
| 01 corpus curate | done | `manifest.parquet` (not in bundle; not needed by viz) |
| 02 SCREEN universe | done | `ccre.bed` (not in bundle) |
| 03 ENCODE bigwigs | done | `bigwigs/...` (not in bundle; consumed by stage 04) |
| 04 intrinsic bigwig | done | columns inside `viz_chr16.parquet` (✓) |
| 05 sequence features | written, not all run | not in bundle (gap) |
| 06 extrinsic | written, not all run | not in bundle (gap) |
| 07 load pretrained R2V | done | `pretrained_universe.parquet` (not in bundle; 498 MB; per-token embeddings live here if the viz wants them) |
| 08 region UMAP + kNN | done (pre-rebase, but corpus-independent so still valid) | `viz_chr16.parquet` (✓) |
| 09 file UMAP | done | `viz_files.parquet` (✓) |
| 10 featured narrative | done | `featured_*.parquet` (✓) |
| 11 tokenize corpus | done | `tokenized_corpus_chr16.parquet` (✓), `file_embeddings.parquet` (not in bundle; intermediate to stage 09) |
| 12 corpus stats | done | `region_stats.parquet` (✓) |
| 13 cooccurrence | done | `region_cooccurrence.parquet` (✓) |

Pipeline source-of-truth: `genomic-dict/config.yaml` and `genomic-dict/pipeline/scripts/`. Read those if a schema or filter behavior is unclear.
