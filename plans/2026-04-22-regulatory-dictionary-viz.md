---
date: 2026-04-22
status: in-progress
description: Course data-viz project scaffolding an interactive "regulatory genomics dictionary" in Observable. Pivoted 2026-04-25 from a retrain-R2V-with-Shape-A/C plan to a pretrained-model demo with a 5-step narrative arc for an informed-layperson audience. Uses databio/r2v-encode-hg38 directly + SCREEN cCRE class tagging; visualized for chr16 with FTO / IRX3, α-globin, MC1R, and CDH1 as featured intervals. See 2026-04-26 update at bottom for current direction.
---

# Regulatory dictionary viz — course project (Shape A)

## Deliverable and deadline

- Target delivery: **2026-05-06** (2 weeks from 2026-04-22).
- Format: hosted Observable notebook + short accompanying writeup (exact format TBD by course requirements).
- Audience: data-viz-literate but not bioinformatics experts. The notebook must read narratively top-to-bottom while supporting interactive exploration.

## The premise (one sentence)

Treat learned region embeddings (R2V-style) as entries in a **dictionary of regulatory genomics**: each token is a place in the genome, with an intrinsic definition (the regulatory machinery present there) and a corpus-level context (which BED-file experiments it appears in alongside which other tokens).

**Motivating hook for the audience**: most GWAS variants sit in noncoding intronic/intergenic space. Interpreting them requires characterizing what regulatory machinery is at those positions — which is what this dictionary is built to surface.

## Goals

Shape A is the **verification step** — build enough of the dictionary metaphor to judge:

1. Does the intrinsic-vs-extrinsic annotation framing read cleanly when applied to real tokens?
2. Does R2V-embedding proximity match regulatory-mechanism similarity (intrinsic annotation overlap)?
3. Is the pipeline sound enough to swap in other chromosomes and (later) other tokenizations/encodings?

If all three land, Shape C (tokenization-fidelity + encoding-resolution panel) is justified as the research-novel extension.

## Design commitments

### Narrative — the dictionary metaphor, rigorously

- **Entry** (single token): its intrinsic regulatory machinery. The "definition."
- **Thesaurus / co-occurrence** (token ↔ token): tokens that travel together across BED files. The "context."
- **Featured entry** (one curated locus): a front-page story anchoring the narrative. **FTO / IRX3** on chr16 — rs1421085 (and neighboring variants) sit in an FTO intron but regulate IRX3 / IRX5 expression megabases away, a canonical demonstration that the nearest coding gene is not the regulatory target. This directly validates the viz's subordination of nearest-gene to intrinsic machinery. Target region roughly `chr16:53,700,000-54,200,000`; exact featured cCRE selected during data prep.

### Intrinsic vs extrinsic annotation — the cardinal commitment

A token's identity is what's **at** that location — the regulatory machinery present, not what gene it's near. Nearest-gene and GWAS-trait overlap are *downstream implications* of that machinery firing, rendered in a visually subordinate section below the primary card.

This commitment is load-bearing: it makes the dictionary metaphor internally coherent, and the featured entry at FTO/IRX3 is a proof-by-demonstration of the commitment itself — the entry shows active-enhancer chromatin marks, motif hits, and conservation (the intrinsic machinery) while the downstream-implications section correctly lists FTO as nearest gene — and yet the regulatory target is IRX3 megabases away. The reader sees why intrinsic > extrinsic without being told.

### Chromosome scope — chr16 for v1, pipeline chromosome-agnostic

Primary scope: chr16 (~90 Mb; **63,677 SCREEN V4 cCREs** in the 5-canonical-class filter). The entire data-prep pipeline takes a chromosome (or region BED file) as input — swapping to chr16 (FTO/IRX3 CAD locus), chr8 (8q24 MYC desert), or a custom region list is a single `config.yaml` edit with no code changes.

Chromosome choice revisited 2026-04-22 after stage 02 ran. Switched from an initial chr16 pick to chr16 because the FTO/IRX3 featured entry directly reinforces the viz's design thesis (intrinsic > extrinsic annotation), whereas FTO/IRX3 shows regulatory pleiotropy but doesn't specifically test the nearest-gene-is-wrong frame.

### Tissues and marks

Defaults (confirmable during data prep):
- Tissues: **K562, GM12878, HepG2** — matches the `2026-04-19-within-ccre-shape-pilot.md` and `2026-04-19-r2v-token-content-comparison.md` panels. Enables annotation re-use.
- Histone marks: standard **5-mark ENCODE core set** — H3K4me1, H3K4me3, H3K27ac, H3K27me3, H3K9me3.
- Accessibility: ATAC-seq fold-change-over-control.

## Panels

### Panel 1 — Dictionary entry

Pick a token → card view with **two clearly separated sections**:

**Intrinsic definition** (the token's identity — primary visual weight):
- SCREEN cCRE class badge (PLS / pELS / dELS / CTCF-only / DNase-H3K4me3)
- Histone-mark heatmap: rows = 5 marks, columns = 3 tissues, color = fold-change signal
- ATAC coverage thumbnail: small line chart per tissue of within-cCRE signal shape
- TF motif tags: JASPAR hits as chips / tag cloud
- ChromHMM state row: state label per tissue (Active TSS / Strong Enhancer / Bivalent / ...)
- Conservation sparkline: phastCons100way over the interval
- R2V stats: "appears in N files," nearest-neighbor cosine distances, embedding-space norm

**Downstream implications** (visually subordinate, framed explicitly as *what this machinery might regulate*):
- Nearest protein-coding gene + one-line GO function
- GWAS-catalog trait overlap (badge when present)
- eQTL overlap (badge when present)

Entry points to this panel:
- Click a point in Panel 2's UMAP
- Search-box lookup by coordinates or nearest-gene name
- "Featured entry" shortcut → the pre-curated FTO/IRX3 cCRE

### Panel 2 — Co-occurrence network / neighborhood

Three composed levels of drill-in:

1. **Overview** — UMAP of all chr16 tokens, colored by SCREEN class. Pan/zoom; hover reveals a mini-entry tooltip.
2. **Class-level co-occurrence** — aggregate tokens to SCREEN class nodes; edges weighted by mean co-occurrence across the BED corpus. Compact, readable 5-node summary network. This is where "tokens of the same regulatory type travel together" becomes visible.
3. **Token neighborhood** — click a token → force-directed graph of its ~30 nearest embedding-space neighbors, overlaid with intrinsic-annotation summary (dominant histone mark; SCREEN class; shared motifs).

Filter/highlight controls: toggle visibility by SCREEN class, dominant mark, tissue activity.

### Panel 3 — "How good is our dictionary?" (A-scoped diagnostic)

One scrollytelling frame, lightly interactive, wrapping the viz with a single verification question: does R2V-embedding proximity match intrinsic-mechanism similarity?

- Scatter: R2V cosine similarity between token pairs vs mechanism-similarity (Jaccard over motif hits; correlation over mark profiles). Sampled pairs, not full pairwise.
- One headline number: fraction of k-NN neighborhoods whose dominant SCREEN class matches the query token.
- Short interpretive paragraph with the two outcomes explicitly stated: if correlation is tight, Shape C asks *how fine can we make this before it breaks?*; if correlation is weak, Shape C asks *does finer tokenization or signal-weighted encoding fix it?*

Panel 3 is the hinge between Shape A's verification and Shape C's future work.

## Data preparation — chromosome-agnostic pipeline

All preprocessing as `snakemake` rules (or equivalent Makefile) parameterized on `CHROM=chr16` (or a region BED file). Stages:

1. **Corpus curation** (genome-wide; not chromosome-dependent):
   - Inspect `hg38_meta_t1.parquet` and `hg38_meta_t2.parquet` from `databio/bedbase-umap` to identify available metadata columns.
   - Filter the ~50k-file metadata to ~20k well-annotated files. Filter criteria set after metadata inspection; minimum bar is tissue + assay + organism confirmed hg38.
   - Retrieve the 20k BED files locally via `geniml.bbclient.BBClient`.
2. **Universe extraction**: download SCREEN cCRE registry (hg38 v4, biosample-agnostic core classification). Keep genome-wide registry as the R2V training vocabulary (cross-chromosomal co-occurrence is essential for meaningful embeddings). Tag the chr16 subset with class labels for viz-time filtering. Restrict to the 5 standard classes.
3. **Intrinsic annotations per chr16 cCRE**:
   - Histone marks: mean fold-change-over-control inside each cCRE, per tissue, per mark. 15 scalars (5 marks × 3 tissues) per cCRE.
   - ATAC shape thumbnail: 100-bin resampled coverage per cCRE per tissue.
   - TF motifs: FIMO scan of cCRE sequences against JASPAR 2024 vertebrate core, `p < 1e-4`. Top-10 hits per cCRE retained for display.
   - ChromHMM state: intersect cCREs with 15-state ChromHMM calls per tissue; majority-vote state per cCRE per tissue.
   - phastCons100way: mean score per cCRE via `bigWigAverageOverBed`.
4. **Extrinsic annotations per chr16 cCRE**:
   - Nearest protein-coding gene: `bedtools closest` against GENCODE v46. Store gene symbol + distance + top GO term.
   - GWAS Catalog overlap: `bedtools intersect`; store trait names (concatenated).
   - GTEx v8 eQTL overlap: same biosamples where available.
5. **R2V training** (custom, not the production `databio/r2v-encode-hg38`):
   - **Why retrain**: the production R2V model has its own token vocabulary baked in, not aligned 1-to-1 with SCREEN cCREs. Training a dedicated model with SCREEN cCREs as the vocabulary means every R2V token *is* a SCREEN cCRE with a clean class label and intrinsic annotation — the dictionary metaphor reads as "one entry = one cCRE" instead of "smear across overlapping cCREs." This is load-bearing for the viz's argument.
   - **Honesty note**: the viz is of a retrained-for-this-project model, not the production one. We state this explicitly in the notebook narration; a brief Panel 3 comparison against the production file-level UMAP (`databio/bedbase-umap`) can corroborate that conclusions generalize.
   - **Vocabulary**: genome-wide SCREEN cCREs (~1M tokens). Chr-restricted vocabularies would starve the embedding of cross-chromosomal context.
   - **Corpus**: the curated 20k BEDs from stage 1.
   - **Tokenization**: cCREs as atomic tokens. Signal encoding: binary active/inactive (Shape A default). Tokenization resolution and encoding become the two Shape C levers.
   - **Embedding dim**: 100.
   - **Implementation**: `geniml.region2vec.main.Region2VecExModel` training pipeline; bedboss/bbconf tokenization code as reference.
6. **Viz artifacts precompute** (chr16 subset only):
   - 2D UMAP of chr16 tokens from trained R2V embeddings.
   - k-NN nearest neighbors (k=30) per chr16 token with cosine distances.
   - Sparse co-occurrence edge list: for each pair of chr16 cCREs, count co-active BED files across the 20k corpus. Threshold by min-count.
   - Aggregated class-level co-occurrence for Panel 2 level 2 (5×5 matrix).

### Payload budget for Observable

- 45k cCREs × (UMAP 2 floats + class byte + ~30 bytes annotation summary) ≈ 2–3 MB gzipped — the always-loaded overview.
- Mark heatmap + ATAC thumbnail per cCRE: lazy-loaded per-cCRE JSON or chunked into ~100-cCRE tiles.
- Co-occurrence edges: thresholded to ~50k–200k edges → ~1 MB.

Target: <10 MB initial load, <1 MB per drill-in fetch.

## Technical stack

- Preprocessing: Python + `pyBigWig`, `pybedtools`, `gtars`, `pymemesuite` (FIMO), `polars` / `pyarrow` for metadata parquets. `uv` for environment management.
- R2V training: `geniml.region2vec` (upstream of `databio/r2v-encode-hg38`). Editable dev installs of the lab's local clones:
  - `/Users/sam/Documents/Work/ai-sandbox/workspaces/sam/bedbase/repos/bbconf`
  - `/Users/sam/Documents/Work/ai-sandbox/workspaces/sam/bedbase/repos/bedboss`
  - `geniml` (transitive; pinned version from the above)
- Corpus retrieval: `geniml.bbclient.BBClient` for BEDbase downloads.
- Notebook: Observable (cloud-hosted).
- Charting: Observable `Plot` for scatters and small multiples; `d3-force` for neighborhood network; custom HTML+CSS for the entry card.
- Data hosting: Observable file attachments for small files, GitHub Pages or S3 for lazy-loaded per-cCRE JSON if payload pushes over Observable's limits.

## Milestones — 2-week plan

**Week 1 — data + Panel 1**
- **Day 1**: Scaffold dict-viz workspace under `spatial-region-features/`. Set up `uv` env with dev installs of bbconf/bedboss/geniml. Inspect `hg38_meta_t*.parquet` to pick corpus filter criteria.
- **Day 2**: Apply filter → curate 20k BEDs → download via BBClient. Download SCREEN cCRE registry; tag chr16 subset. Smoke test preprocessing on chr22 to verify chr-agnostic parameterization.
- **Day 3–4**: Run intrinsic + extrinsic annotation pipeline on chr16. Produce all annotation tables.
- **Day 5**: Train R2V with SCREEN cCRE vocabulary on the 20k corpus. Extract chr16 token embeddings. Precompute UMAP + k-NN + co-occurrence.
- **Day 6–7**: Panel 1 (dictionary entry) live in Observable. Curate the FTO/IRX3 featured entry.

**Week 2 — Panel 2, Panel 3, polish**
- **Day 8–9**: Panel 2 overview (UMAP) + class-level co-occurrence network.
- **Day 10**: Panel 2 token-neighborhood drill-in.
- **Day 11**: Panel 3 diagnostic (embedding-mechanism correlation scatter + k-NN class purity number).
- **Day 12–13**: polish, scrollytelling copy, featured-entry walkthrough, deployment.
- **Day 14**: buffer for the unexpected.

## Open decisions

**Resolved 2026-04-22**
- Shape: A first (dictionary browser + co-occurrence network). Shape C (tokenization/encoding fidelity) as follow-on contingent on A.
- Chromosome scope: chr16 for v1; pipeline chr-agnostic. Featured entry at FTO/IRX3.
- Tissues: K562, GM12878, HepG2. Histone marks: 5-mark core set.
- R2V model: retrain with SCREEN cCRE vocabulary on curated 20k corpus (Option B). Pretrained `databio/r2v-encode-hg38` kept as reference, not primary.
- Corpus source: `databio/bedbase-umap` metadata parquets → `geniml.bbclient.BBClient` fetch.

**Still to resolve before or during Day 1**
1. **Corpus filter criteria, concrete.** Depends on what columns `hg38_meta_t*.parquet` actually has. Minimum bar: tissue + assay + organism hg38. Refine after inspection.
2. **ChromHMM source** — Roadmap 15-state vs ENCODE 18-state. Pick whichever has hg38 calls for all three tissues readily available.
3. **TF motif database** — JASPAR 2024 vertebrate core (~800 motifs), top-10 hits per cCRE displayed. Alternative: curated 100-motif "canonical TF" subset for cleaner tags.
4. **GWAS/eQTL trait filter** — show all overlaps, or curate to chr16-relevant traits (CAD / T2D / cancers) for the featured entry caption? Proposal: show all as badges, curate only for the featured-entry narration.
5. **Course-deliverable format beyond the notebook** — confirm the course's exact requirements (writeup length? video walkthrough? pre-recorded demo?).

## Bridge to Shape C — soft tokenization (aligned to PI's grant Aim 3.3)

Read of the PI's grant proposal (2026-04-22) reshaped Shape C's framing. The grant's Aim 3 names three tokenization innovations:

- **3.1 Hierarchical tokenization** — multi-level L1/L2 universe sieve (addresses out-of-vocabulary)
- **3.2 Metatokens** — group tokens by co-accessibility/co-occurrence/annotation (addresses vocabulary size)
- **3.3 Soft tokenization** — kernel-weighted overlap instead of binary present/absent (addresses token mapping / information loss)

**Shape C commits to sub-aim 3.3 — soft tokenization**, because:
- It's the most direct fix for the weakness Panel 3 exposes (binary R2V collides width with signal strength).
- Same vocabulary as Shape A (SCREEN cCREs on chr16 for viz; genome-wide for training), so the comparison is clean.
- The viz UI needs no structural change — one switch between binary and soft per-token weight surfaces.

### Concrete mechanism

- Shape A: binary active/inactive token activations in R2V training.
- Shape C: per-file, per-token kernel-weighted activation score (e.g., fraction of token overlapped by file regions, or Gaussian-weighted distance from center). Soft activation feeds R2V training directly.

### Why Shape A's Panel 3 diagnostic motivates Shape C

- Expected Shape A outcome: weak correlation between R2V proximity and intrinsic-mechanism similarity, because binary R2V can't distinguish "wide + weak" from "narrow + strong" regulatory signatures.
- Shape C's soft tokenization separates these in the activation space without changing the vocabulary. If correlation rises, soft tokenization's value is demonstrated directly on a real regulatory corpus.

### Stretch: Shape C+ (metatokens, aim 3.2)

Only if time permits. Group tokens by co-occurrence patterns (easily derivable from R2V's own output) to test whether a reduced metatoken vocabulary preserves class-separation.

### What this project contributes

Scoping the project this way positions it as a **visual demonstration of the grant's Challenge 2 (interpreting genomic region embeddings) and Aim 3.3 (soft tokenization)**, not as a tangential course exercise. The featured-entry card shows exactly what a cCRE's "definition" looks like in embedding space; the Panel 3 diagnostic shows the binary-vs-soft transition's effect on how well that definition reads as biologically coherent.

## Relationship to the broader workspace

- **Reuses within-cCRE signal-extraction infrastructure** from `2026-04-19-within-ccre-shape-pilot.md`. The mark heatmap and ATAC thumbnail apply the same within-interval extraction pattern.
- **Parallels `2026-04-19-r2v-token-content-comparison.md`** at the per-cCRE level. That pilot asks "does R2V proximity predict cCRE class better than content features?"; Panel 3 is the visual relative of that test.
- **Operationalizes the intrinsic-vs-extrinsic framing** in a way the dissertation proper can point to as a calibration artifact — an interactive way for readers to see what an embedding has learned about a real noncoding-regulatory locus.

## Status and next steps

**Draft.** Before executing:

1. Lock in open decision (1) — BED corpus source and filter.
2. Confirm default selections (tissues, marks, ChromHMM source) are acceptable.
3. Confirm course-deliverable format (open decision 6).
4. Quiz checkpoint (per user preference) on the load-bearing design decisions before starting implementation.

---

## 2026-04-26 — Major pivot away from training, toward 5-step narrative demo

After working through what R2V training would actually contribute beyond annotation-based similarity (and what a naive "dictionary of regulatory genomics" reader would expect), the project pivoted away from the train-our-own-model path. Key realizations:

1. **R2V's contribution isn't load-bearing for the dictionary metaphor.** The dictionary is carried by per-entry richness (Panel 1 entry card), which uses annotations from ENCODE/SCREEN/JASPAR/UCSC — not the embedding. R2V adds a co-occurrence-based distance metric, useful for Panel 2's UMAP layout, but a pretrained model gives us that without retraining.
2. **A naive "dictionary of regulatory genomics" reader expects a reference resource**, not a methods demo. Broad coverage, rich per-entry context, mechanism-level navigation, downstream-implication links, memorable hooks. The R2V/Shape A/C framing was methods-paper, not encyclopedia.
3. **The pretrained `databio/r2v-encode-hg38` model already has 1.06M region embeddings** that overlap SCREEN cCREs at 99.7%. Using it directly cuts ~6 hours of training + ~2 hours of tokenization off the critical path with zero loss to the dictionary's Panel 1 richness.
4. **Shape C (soft tokenization, grant Aim 3.3) drops from the project**, since it requires retraining. Future work, not this 2-week course demo.

### New direction: 5-step narrative demo

Course goal recast (2026-04-26): build a demo viz showing what *"a dictionary of regulatory genomics"* might look like for an informed-layperson audience. Five-step narrative:

1. **Here are some BED files. What's in them?** Featured intervals on chr16; show BED tracks (peaks + optional bigwig signal) at those intervals across multiple modalities (DNase, ATAC, ChIP, etc).
2. **We tokenize BED files to a common vocabulary through a universe.** Show universe regions at the same intervals; show each BED file becoming a binary presence vector against that vocabulary.
3. **What does it mean to treat tokenized regions as inputs to an NLP-style model?** Co-occurrence of regions, the idea of a "grammar" of regulatory elements; introduce the SCREEN class taxonomy and the region UMAP from the pretrained R2V embedding.
4. **Click a BED file → highlight the regions it activates** on the region UMAP. Shows how a single experiment relates to the broader vocabulary.
5. **Open-ended interpretation**: beyond ENCODE annotations, how do we hypothesize what these regions are *for*? Separate plan: `plans/2026-04-26-region-interpretation-step5.md`.

### What stages survive / change / drop

**Survive** (outputs already on disk):
- Stage 00: HF metadata inspection
- Stage 01: Tight corpus curation (84,698 files)
- Stage 02: SCREEN V4 universe + metadata
- Stage 03: 18 ENCODE bigwigs on Rivanna
- Stage 04: bigwig features per chr16 region (now keyed on pretrained regions, not cCREs)
- Stage 07 (new): load pretrained R2V universe + SCREEN-overlap tagging
- Stage 08: chr16 viz precompute (UMAP + kNN)
- Stage 09: file-level UMAP from HF coords + 20 mystery entries

**Dropped**:
- Old stage 07 (tokenize 85k corpus) — no training corpus needed
- Old stage 08 (R2V training) — pretrained model used directly
- Stage 05 (sequence features: FIMO + ChromHMM + phastCons) — deferred; entry card is rich enough for MVP without
- Stage 06 (extrinsic: nearest gene + GWAS + eQTL) — deferred

**Added (current work)**:
- Stage 10: featured-narrative data prep — curated featured intervals + featured BED files; per-file chr16 tokenization for Panel 4 highlighting; per-(file, interval) BED peak slicing for Panel 1 tracks.

### Featured items

**Featured intervals** (chr16):
- FTO / IRX3 locus (chr16:53,700,000–54,200,000) — flagship "nearest gene is wrong"
- α-globin cluster (chr16:170,000–230,000) — classic regulatory gene cluster
- MC1R locus (chr16:89,900,000–90,000,000) — pigmentation, well-studied
- CDH1 / E-cadherin region (chr16:68,700,000–68,900,000) — tumor suppressor with regulatory complexity

**Featured BED files** (~18 + 2–4 mystery):
- 18 ENCODE files: K562 / GM12878 / HepG2 × (ATAC, H3K4me1, H3K4me3, H3K27ac, H3K27me3, H3K9me3) — match between bigwig signal and peak BED comes free
- 2–4 mystery (UNKNOWN cell_line) entries from `viz_files.parquet`'s curated sample, for step 5's hypothesis-generation hook

### Observable scaffold

`genomic-dict/notebook/scaffold.md` written 2026-04-26 — copy-pasteable cells covering data loading, shared state, file UMAP (Panel A), region UMAP (Panel B), entry card placeholder (Panel C). User drives panel design from there.

### Status

Backend pipeline complete through stage 09. Stage 10 in flight. Notebook scaffold confirmed working. User implements panels and narrative copy from here.
