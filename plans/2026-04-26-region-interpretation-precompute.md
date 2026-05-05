---
date: 2026-04-26
status: draft
description: Data-side plan for richer region interpretation in the genomic-dict viz. Three pillars — Path A (external annotation lookup), Path B (corpus-derived BED-level statistics), and region co-occurrence networks (with cell-line context stratification for "contextual grammar" visualization). Assumes the Observable Framework port at https://github.com/sanghoonio/genomic-regions.git is already in place; this plan only adds new precomputed parquets to feed into existing or new panels.
---

# Region interpretation precompute — Path A + Path B + co-occurrence networks

## Goal

Add the data layer that backs three interpretation patterns the viz can render:

1. **Path A — external annotation lookup**: per region, what's been catalogued (TF motifs, ChromHMM state, conservation, nearest gene, GWAS, eQTL). For Step 5 hypothesis cards and richer Step 1/2 entry tooltips.
2. **Path B — corpus-derived statistics**: per region, what the 84k-file corpus says (cell-line activation frequency, assay activation frequency, tissue activation frequency, "pleiotropy" / specificity scores). For sub-class structure visible in the UMAP and quantitative claims in hypothesis text.
3. **Region co-occurrence networks**: per region, top-K most-co-activating partners. Both corpus-wide and **stratified by cell line** so the viz can show *the same region's neighborhood reorganizing across contexts* — visual evidence for context-dependent grammar.

This plan stays in the data layer. The Framework viz code is out of scope — it consumes the parquets we produce.

## Pillars and dependencies

```
┌── Path A ──────────────────┐
│  Stage 05 sequence features │  (FIMO motifs + ChromHMM + phastCons)
│  Stage 06 extrinsic         │  (nearest gene + GWAS + eQTL)
└─────────────────────────────┘
              ↓ (chr16 region annotations)
              join with viz_chr16.parquet downstream

┌── Path B + Networks ───────┐
│  Stage 11 tokenize 84k corpus  (chr16 only; deferred from earlier pivot)
│       ↓
│  Stage 12 corpus statistics    (cell-line / assay / tissue frequency tables)
│       ↓
│  Stage 13 region co-occurrence networks  (corpus-wide + stratified)
└────────────────────────────┘
```

Path A can run independently of Path B. Path B's tokenization is the biggest single compute item (~3 hours on Rivanna login node, network-bound). Stage 11 → 12 → 13 are sequential.

## Reference paths

- **Pipeline scripts**: `/Users/sam/Documents/Work/spatial-region-features/genomic-dict/pipeline/scripts/`
- **Existing parquets**: `/Users/sam/Documents/Work/spatial-region-features/genomic-dict/data/`
- **Observable Framework target repo**: https://github.com/sanghoonio/genomic-regions.git
- **Rivanna env**: `source /project/shefflab/rivanna_config/env.sh` (provides `BBCLIENT_CACHE`)
- **Existing config**: `genomic-dict/config.yaml` (extend with new stages)

## Path A — Stage 05: sequence-based intrinsic features

Per chr16 region, derive features from sequence + integrative annotations.

### Inputs

- `data/annotations/pretrained_universe.parquet` (chr16 subset; coords)
- New downloads: JASPAR 2024 vertebrates MEME (~25 MB), Roadmap ChromHMM 15-state for K562/GM12878/HepG2 (3 BEDs, ~50 MB total), UCSC phastCons100way bigwig (~10 GB)
- Reference genome FASTA for chr16 (need for FIMO sequence input). hg38 chr16 ~90 MB.

### Outputs

`data/annotations/intrinsic_sequence.parquet` — one row per chr16 region:

| column | type | source |
|---|---|---|
| `token_id` | int | join key with regions table |
| `top_motifs` | list[string] (top-10) | JASPAR + FIMO p<1e-4 |
| `top_motif_pvalues` | list[float] | matching pvalues |
| `phastcons_mean` | float | mean phastCons100way over the region |
| `chromhmm_K562` | string | majority-vote state |
| `chromhmm_GM12878` | string | majority-vote state |
| `chromhmm_HepG2` | string | majority-vote state |

### Compute

- FIMO scan: 35,934 regions × ~250 bp avg sequence × ~800 motifs. With pymemesuite, parallelizable. ~30 min on a single node.
- ChromHMM intersect: bedtools intersect, fast (<5 min).
- phastCons mean: `bigWigAverageOverBed` per region, fast (<5 min).

Total runtime: ~45 min.

### SLURM

`pipeline/slurm/05_extract_sequence.sh` — 4 GB mem, 1 cpu, 2 hr cap.

## Path A — Stage 06: extrinsic annotations

Per chr16 region, derive nearest-gene + disease associations.

### Inputs

- `data/annotations/pretrained_universe.parquet` (chr16 subset)
- New downloads:
  - GENCODE v46 GTF (~50 MB)
  - GWAS Catalog full TSV (~10 MB; from EBI)
  - GTEx v8 cis-eQTL TAR (~5 GB; only need significant pairs subset)

### Outputs

`data/annotations/extrinsic.parquet` — one row per chr16 region:

| column | type | source |
|---|---|---|
| `token_id` | int | join key |
| `nearest_gene_symbol` | string | bedtools closest, GENCODE |
| `nearest_gene_distance` | int | bedtools closest |
| `nearest_gene_top_go` | string | mygene.info lookup of top GO BP term |
| `gwas_traits` | list[string] | GWAS Catalog overlap |
| `gtex_eqtl_genes` | list[string] | GTEx significant eQTLs targeting which genes |

### Compute

- bedtools closest: <1 min
- mygene.info lookup: 35k API calls — batchable, ~30 min with rate limiting
- GWAS Catalog intersect: bedtools intersect, <1 min
- GTEx tar processing: parse + filter to chr16 + intersect, ~30 min

Total runtime: ~1 hour.

### SLURM

`pipeline/slurm/06_extract_extrinsic.sh` — 16 GB mem (GTEx tar parsing), 1 cpu, 2 hr cap.

## Path B — Stage 11: tokenize 84k corpus (chr16 only)

The deferred step. Tokenize all corpus BED files against the pretrained chr16 universe via bbcache.

### Inputs

- `data/corpus/manifest.parquet` (84,698 file IDs)
- Pretrained tokenizer (`databio/r2v-encode-hg38` via geniml; same as stage 10)
- chr16 universe filter (token IDs from `pretrained_universe.parquet` where chrom = "chr16")

### Outputs

`data/annotations/tokenized_corpus_chr16.parquet` — one row per BED file:

| column | type |
|---|---|
| `id` | string (BEDbase file ID) |
| `chr16_active_token_ids` | list[int] (sorted, deduplicated) |
| `n_chr16_active` | int |

Streaming-write in batches of 5000 to keep memory bounded (same fix we applied to the original stage 07 tokenization).

### Compute

- bbcache load_bed: ~0.1s per file × 84k = ~2.5 hr
- Tokenizer encode: ~10ms per file × 84k = ~15 min
- Total: ~3 hr on Rivanna login node, network-bound

### SLURM

`pipeline/slurm/11_tokenize_corpus_chr16.sh` — 16 GB mem, 1 cpu, 6 hr cap.

### Why chr16-only

The full-genome tokenization would be ~6× slower (the chr16 universe is ~17% of the full universe by token count, but the bottleneck is `bbc.load_bed`, which loads the whole BED regardless). For the demo we only need chr16 statistics; restricting the *output* to chr16 cuts memory and disk by 5×.

## Path B — Stage 12: corpus statistics

Aggregate the tokenized corpus into per-region frequency vectors stratified by metadata.

### Inputs

- `data/annotations/tokenized_corpus_chr16.parquet` (from stage 11)
- `data/precomputed/viz_files.parquet` (file metadata: cell_line, assay, tissue)

### Outputs

`data/annotations/region_corpus_stats.parquet` — one row per chr16 region:

| column | type | meaning |
|---|---|---|
| `token_id` | int | join key |
| `n_files_total` | int | total files activating this region (out of 84k) |
| `freq_total` | float | n_files_total / 84k |
| `freq_by_cell_line` | struct{K562, GM12878, HepG2, MCF-7, ...} | per-cell-line activation freq |
| `freq_by_assay` | struct{DNase-seq, ATAC-seq, ChIP-seq, ...} | per-assay activation freq |
| `freq_by_tissue` | struct{blood, brain, ...} | per-tissue activation freq |
| `dominant_cell_line` | string | argmax of freq_by_cell_line |
| `dominant_assay` | string | argmax of freq_by_assay |
| `tissue_pleiotropy` | int | number of distinct tissues with freq > threshold |
| `tissue_specificity` | float | normalized entropy of freq_by_tissue (0=specific, 1=uniform) |

### Compute

Inverted index: file → token list → per-token list of files. Then groupby per (region, cell_line) etc. Polars can do this in memory if we shape it right; otherwise sparse matrix via scipy.

Estimated runtime: ~20 min for the full aggregation given a sparse-matrix approach.

### SLURM

`pipeline/slurm/12_corpus_stats.sh` — 32 GB mem (sparse matrix + groupby), 1 cpu, 1 hr cap.

## Path B — Stage 13: region co-occurrence networks

For each chr16 region, compute its top-K most-co-activating partner regions, both corpus-wide and stratified by cell line.

### Inputs

- `data/annotations/tokenized_corpus_chr16.parquet`
- `data/precomputed/viz_files.parquet` (for cell-line stratification)

### Outputs

`data/annotations/region_cooccurrence.parquet` — one row per (region, context) pair:

| column | type | meaning |
|---|---|---|
| `token_id` | int | source region |
| `context` | string | "all" or cell_line name |
| `partner_token_ids` | list[int] (top-K=30) | partner regions, sorted by co-occurrence count |
| `partner_counts` | list[int] | matching co-activation counts |
| `partner_jaccard` | list[float] | normalized score: count / |files where either is active| |

So for 35k regions × 4 contexts (all + 3 main cell lines) = 140k rows. Compact since each row has 30 entries max.

### Compute

Sparse file × token matrix, then matrix multiplication with top-K selection per row:

- file × token matrix: 84k × 35k binary, sparse (~0.5% density)
- Multiply X.T @ X gives token × token co-occurrence (35k × 35k = 1.2B entries, but sparse)
- For top-K: row-wise sort + slice
- Stratified version: subset rows of X by cell_line first, then multiply

Approximate runtime:
- All-files matrix mult + top-K: ~30 min on a single node with scipy sparse
- Per cell line × 3: another ~30 min
- Total: ~1 hr

### SLURM

`pipeline/slurm/13_cooccurrence.sh` — 32 GB mem, 4 cpu (parallel cell-line strata), 2 hr cap.

## Summary table — all new stages

| stage | runtime | mem | output | size |
|---|---|---|---|---|
| 05 sequence | 45 min | 4 GB | intrinsic_sequence.parquet | ~3 MB |
| 06 extrinsic | 1 hr | 16 GB | extrinsic.parquet | ~5 MB |
| 11 tokenize 84k | 3 hr | 16 GB | tokenized_corpus_chr16.parquet | ~500 MB |
| 12 corpus stats | 20 min | 32 GB | region_corpus_stats.parquet | ~5 MB |
| 13 cooccurrence | 1 hr | 32 GB | region_cooccurrence.parquet | ~30 MB |

Total wall-time: ~6.5 hr if run sequentially. Path A and stages 11→12→13 can run in parallel since Path A doesn't need stage 11.

## Integration with the Framework viz

Three new parquets to drop into the genomic-regions repo's `src/data/` (or wherever it loads from):

- `intrinsic_sequence.parquet` + `extrinsic.parquet` → join into the regions table; populate Path 1 entry tooltips and Step 5 hypothesis text
- `region_corpus_stats.parquet` → drives a new "color UMAP by tissue bias / assay bias / pleiotropy" toggle on Step 3
- `region_cooccurrence.parquet` → drives a new "co-occurrence network around selected region" panel; allows toggling between corpus-wide / K562-only / GM12878-only / HepG2-only views

Schema note: all three keyed on `token_id` from `pretrained_universe.parquet`. The Framework code joins by token_id where needed.

## Order of operations

Suggested execution sequence in a fresh session:

1. **Path A in parallel** (cheaper, lower-risk): stage 05 + 06. Outputs in ~2 hr.
2. **Path B sequential**: stage 11 (tokenize, 3 hr) → stage 12 (stats, 20 min) → stage 13 (cooccurrence, 1 hr). Total ~4.5 hr.
3. Pull all 5 parquets back to local; drop into the genomic-regions repo.

If time-constrained: do stages 05 + 11 + 13 first (the most distinctive viz capabilities). Stages 06 + 12 add finer granularity but the demo can ship without them.

## Caveats

- **GTEx tar is ~5 GB**. Large initial download. Worth caching at `data/annotations/gtex_v8.tar` so reruns don't re-download.
- **mygene.info rate limiting** for stage 06's GO term lookups. Use the batch API (`/v3/gene` POST with multiple IDs) and respect ~1000 reqs/min.
- **Stage 11 OOM risk**. The original 85k-tokenization OOM'd at the unbounded-list stage. Re-use the streaming `pq.ParquetWriter` pattern with batch_size=5000. Per-file lists should never accumulate beyond one batch.
- **Cell-line subset sizes** for stage 13: K562 ~14k files, HepG2 ~7k, GM12878 ~4k. Ratios skewed but each subset is large enough to compute reliable top-K co-occurrence.
- **Stage 13 size on disk**: 30 MB compressed seems reasonable but list-of-int columns can balloon. If too large, drop `partner_counts`/`partner_jaccard` (compute on the fly in viz from token_id pairs).

## Status / next steps

**Draft.** Self-contained for fresh-session execution. Steps:

1. Read this plan
2. Read the existing pipeline conventions in `genomic-dict/pipeline/scripts/_common.py` (config-driven, summary.json per stage)
3. Add config sections for stages 05, 06, 11, 12, 13 in `genomic-dict/config.yaml` (model after existing stages)
4. Write the 5 new scripts + their SLURM wrappers
5. `./sync.sh` and submit jobs
6. After each lands: pull the parquet, drop it into the genomic-regions repo
