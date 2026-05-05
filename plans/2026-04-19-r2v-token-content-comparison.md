---
date: 2026-04-19
status: draft
description: Companion pilot to the within-cCRE shape experiment — compare R2V's learned embeddings (aggregated per cCRE) against within-cCRE content representations on a curated 3-tissue panel, with leave-one-cCRE-out attribution to understand what R2V implicitly weights when aggregating to file-level.
---

# R2V per-token content comparison

**Goal:** decompose what R2V's learned embeddings represent at per-cCRE granularity — specifically whether they carry more, less, or different biology than within-cCRE content vectors, and what cCRE-level properties R2V implicitly weights when producing file-level embeddings.

## The question

Given cCREs with known biological labels (SCREEN cCRE class), how do R2V-derived embeddings (aggregated per cCRE) compare to content-based representations computed inside each cCRE?

Three specific sub-questions:

1. Which representation — R2V-derived embedding or content vector — predicts cCRE class better on held-out cCREs?
2. Do R2V and content recover the same cluster structure among cCREs, or different structures?
3. What content properties does R2V implicitly upweight when aggregating cCREs to file-level embeddings?

## Why this is the right second experiment

1. **Direct upgrade from the retired `bed-embedding-eval` dense track** (see retrospective at `../retired/bed-embedding-eval/README.md`). That track correlated R2V file-level similarities against file-level geometric features and found ρ 0.88–0.98, but could not separate "R2V learns geometry" from "R2V learns biology with geometric signatures." This experiment moves the evaluation to the **per-cCRE level**, where the universe (SCREEN cCREs) provides biology-anchored ground truth (cCRE class per interval) and the confound resolves.
2. **Reuses infrastructure from the first pilot.** Content vectors come from the same ATAC bigwig + within-cCRE feature extraction as `./2026-04-19-within-ccre-shape-pilot.md`. The R2V model is already trained (`r2v-scatlas-hg38-v2`). No training, no new pipeline beyond what the first pilot already needs.
3. **Pairs naturally with the first pilot.** The first pilot establishes whether content alone carries cCRE-class signal. This experiment asks whether R2V adds anything above what content already provides — or matches it, or underperforms it.
4. **Each outcome is decision-relevant.** If R2V is redundant with content, the lab's token-based embedding stack has a ceiling. If R2V adds signal above content, we learn *what* it's adding (via Step C's influence-vs-content correlation). If R2V underperforms content on biology-anchored tasks, the within-universe content direction earns empirical support.

## Method

### Curated panel

- **Universe**: **ENCODE SCREEN cCRE registry (hg38, biosample-agnostic core classification)** — the same universe as the first pilot, so both plans operate on the same biology-anchored intervals. ~1 M cCREs, each with a class label (PLS, pELS, dELS, CTCF-only, DNase-H3K4me3).
- **R2V source**: the pretrained `r2v-scatlas-hg38-v2` model provides embeddings for Atacformer scATLAS tokens, not for cCREs directly. To get a per-cCRE R2V embedding, aggregate: for each cCRE, average the R2V embeddings of all scATLAS tokens whose intervals overlap the cCRE. cCREs that overlap zero scATLAS tokens are excluded from analysis.
- **Content signal source**: ENCODE ATAC bigwig fold-change-over-control for **K562, GM12878, HepG2** — the same three tissues as the first pilot, for direct comparability.
- **BED file corpus for Step B**: curated set of ~20–50 BED files across the three tissues, drawn from BEDbase (multiple assays per tissue, multiple replicates per assay where available).

### Per-cCRE representations

- **(a) R2V-derived embedding**: mean of the R2V embeddings of scATLAS tokens that overlap the cCRE. 100-dim vector per cCRE. This is an aggregation, not R2V-native, so the experiment is strictly a test of R2V's per-token structure *as surfaced at the cCRE granularity*. If scATLAS-level behavior differs meaningfully from cCRE-aggregate behavior, that's an interesting separate question — a follow-up could add a dual scATLAS-level comparison if Step A produces signal.
- **(b) Content vector**: for each cCRE, compute the 5-feature content vector using the first pilot's pipeline (summit position, summit sharpness, signal concentration, edge-to-core ratio, total signal). Average across the three tissues. 5-dim vector per cCRE. **Directly reuses the first pilot's extraction output**, no re-computation needed.

Tissue-specific content vectors (one per cCRE per tissue) are an optional extension if Step A shows signal.

### Step A: class-separation comparison

All classifiers use **stratified train/test splits with chromosome-based folds** (hold out whole chromosomes to prevent data leakage through nearby genomic intervals) and **class-weighted macro-F1** (SCREEN class distribution is not uniform — PLS is typically much smaller than dELS).

1. Train multinomial logistic regression: **R2V-derived embedding → cCRE class**. Report macro-F1.
2. Train **content vector → cCRE class** on the same split. Report macro-F1.
3. Train on **concatenation [R2V; content] → cCRE class**. Report macro-F1. Non-zero gain over each individually indicates complementarity.
4. Cluster R2V-derived embeddings (k-means with k tuned via silhouette, or UMAP + HDBSCAN). Separately cluster content vectors. Measure cluster-assignment agreement via Adjusted Rand Index and Normalized Mutual Information.

### Step B: leave-one-out cCRE attribution

For each BED file in the curated corpus:

1. Tokenize against the scATLAS universe (how R2V was actually trained) → active scATLAS-token set. Compute the file-level R2V embedding as the mean of active scATLAS-token embeddings (the standard R2V set-vector aggregation).
2. Group the active scATLAS tokens by which cCRE they overlap. For each overlapped cCRE C, the "file's C-contribution" is the subset of active tokens that fall inside C.
3. For each overlapped cCRE C, recompute the file-level embedding **with C's contributing tokens removed** — i.e., mean of active scATLAS tokens that do not overlap C. This gives a leave-one-cCRE-out file embedding. The delta vector against the full file embedding captures C's contribution.
4. Store `|delta_C|` (L2 norm) as cCRE C's influence in this file.

Aggregate across files: for each cCRE C that has at least M file-contributions (M ≈ 3 to avoid noisy estimates), compute `influence(C) = mean(|delta_C|)` over files where C contributes. cCREs with high influence are "load-bearing" for R2V's file-level representations.

### Step C: tie Steps A and B together

All correlations below are over cCREs, using the `influence(C)` values from Step B.

1. **R2V influence vs content properties.** Per-feature Spearman correlation of `influence(C)` against each content-vector feature (summit position, sharpness, etc.). Which content properties correspond to high-influence cCREs? Also inspect scatter plots per feature — the relationship may not be monotonic.
2. **R2V influence vs cCRE class.** Stratify the `influence(C)` distribution by cCRE class. Does R2V give systematically higher influence to specific classes (e.g., PLS)? Use class-weighted summary statistics (median, IQR per class) rather than raw means, since the class distribution is imbalanced.
3. **R2V influence vs embedding properties.** Correlate `influence(C)` with R2V-derived intrinsic properties: distance to the R2V-embedding centroid, local density in R2V-embedding space. Does R2V upweight "distinctive" cCREs?

### Informative outcomes (not pass/fail criteria)

This experiment decomposes R2V's cCRE-aggregated semantics; any clean decomposition is useful. Rough interpretive thresholds:

- **R2V F1 > content F1 by > 0.05**: R2V's distributional signal (even after cCRE aggregation) carries biology beyond within-cCRE content. Worth investigating what.
- **R2V F1 ≈ content F1** (within 0.03): the two representations converge — R2V at cCRE granularity is essentially a content summarizer via the distributional route.
- **R2V F1 < content F1 by > 0.05**: R2V's distributional signal is noisier than content alone on this class-separation task. Strong empirical support for the within-universe content direction.
- **ARI / NMI low**: the two representations genuinely see different structures. The disagreement cCREs are the interesting subset for case-study analysis.
- **Step C reveals implicit R2V priors**: if high-influence cCREs share a characteristic content profile, R2V systematically upweights certain regulatory-element types when producing file-level embeddings.

### What each outcome teaches

- R2V beats content: distributional signal carries real biology; worth dissecting what R2V is learning that content misses.
- R2V matches content: the lab's R2V-based tools are content summarizers with extra steps — useful for production, not fundamentally adding biology beyond what a within-universe content representation provides.
- R2V underperforms content: the within-universe spatial-features direction of this workspace earns strong empirical support against the distributional alternative.

## Tractability

- **Data**: same bigwigs as the first pilot, plus the R2V model (a few hundred MB), plus the curated BED file corpus (~1 GB).
- **Compute**: all CPU. Step A is standard logreg + clustering on ~1 M cCREs × 100 dims (minutes). Step B is O(Σ_files k_file × c_file) where c_file is the number of cCREs overlapping active tokens in each file; still cheap with the closed-form mean update (computing the mean-over-active-tokens minus the subset-over-a-cCRE is linear per file per cCRE). Step C is correlation computation.
- **Timeline**: 1 week if the first pilot's content-extraction pipeline already exists; 1–2 weeks starting from scratch. Cleanest sequence: run the first pilot first, then this one, reusing the pipeline.

## Open design decisions

1. **Curated BED file corpus for Step B.** Need a concrete list: how many files per tissue, which assays to include (ATAC only? or ATAC + H3K27ac + H3K4me3 to probe assay-invariance of R2V's cCRE weighting?), replicate handling. Propose: 5–10 ATAC files per tissue for the first pass, add other assays if results warrant.
2. **Tissue-averaged vs per-tissue content vectors.** Draft uses tissue-averaged. Per-tissue content is more informative but complicates the comparison with R2V (which has no tissue dimension). Run tissue-averaged first; add per-tissue as a follow-up.
3. **R2V aggregation method for file embeddings.** Mean of active scATLAS-token embeddings is the usual default for R2V; the retired `bed-embedding-eval` pipeline might have used TF-IDF-weighted aggregation or similar. Check its retrospective before starting.
4. **Class-imbalanced evaluation** (addressed in Step A and Step C above; noted here for completeness). Applied consistently throughout.
5. **cCRE-aggregated vs scATLAS-native R2V evaluation.** The draft aggregates R2V to cCRE granularity via mean-over-overlapping-scATLAS-tokens. A parallel scATLAS-native analysis (R2V embedding per scATLAS token, cCRE class inherited via intersection, tokens that don't overlap any cCRE excluded) tests whether the aggregation is lossy. Run cCRE-aggregated first; add scATLAS-native as a follow-up if Step A produces signal and the aggregation bias is worth characterizing.

## Relationship to the broader workspace

- **Direct sequel to the retired `bed-embedding-eval` dense track.** Moves evaluation from file-level correlation to token-level content comparison, resolving the geometry-vs-biology confound that retired that project.
- **Reuses content-extraction infrastructure from the first pilot** (`./2026-04-19-within-ccre-shape-pilot.md`). The first pilot builds content features against SCREEN cCRE classes; this experiment uses them as the baseline against which R2V is compared.
- **Addresses Q1 cluster (F) sub-questions F1 and F2** from [`../../dissertation/questions/q1-profile-predictivity.md`](../../dissertation/questions/q1-profile-predictivity.md). F1 (do per-token representations align with known regulatory categories) is answered by Step A. F2 (decompose embedding similarity into interpretable axes) is partially addressed by Step C's influence-vs-content analysis.
- **Tests the `embedding-interpretability.md` framing empirically.** That framing argues token co-occurrence is *circular as a semantics source* — two tokens are similar because they co-occur, they co-occur because they play similar roles, but the roles are exactly what the embedding does not expose. This experiment tests whether within-cCRE content representations expose class structure that R2V's distributional embedding does not — confirming, refuting, or refining the framing.

## Status and next steps

**Draft.** Before executing:

1. Pin down the BED file corpus for Step B (open decision 1) — choose accessions and confirm BEDbase access.
2. Confirm R2V aggregation method (open decision 3) by checking the retired `bed-embedding-eval` retrospective and source code.
3. Decide whether to run tissue-averaged content first (as drafted) or tissue-stratified from the start.

**Sequencing relative to the first pilot.** These experiments can run in parallel logically, but the first pilot's content-extraction pipeline is the upstream dependency. Cleanest path: build the pipeline once for the first pilot, reuse for this experiment, so the two plans share their infrastructure cost.

## 2026-04-19 update — LanceOtron as architectural precedent, sharpened hypotheses

The literature check in [`../../dissertation/background/within-peak-shape.md`](../../dissertation/background/within-peak-shape.md) surfaced a useful conceptual anchor for this pilot. Hentges et al. 2022 (LanceOtron) trains a **"wide and deep"** model that combines two parallel input branches for the same region — a CNN on continuous shape and a logistic regression on aggregate enrichment scalars — fused by an MLP. The paper shows combining the two beats either alone for peak-vs-noise discrimination.

**Transferred to the R2V setting**, this pilot's three-way comparison (R2V-derived embedding, within-cCRE content vector, concatenated `[R2V; content]`) is the same design principle applied at a different granularity:

- R2V corresponds to LanceOtron's aggregate/corpus-level branch (cross-file co-occurrence is the "wide" view — many tokens, summarized relationally).
- Within-cCRE content corresponds to LanceOtron's shape branch (within-interval continuous structure is the "deep" view — one token, fine-resolution).
- The two branches encode **orthogonal** information: R2V cannot see shape (binary BED input); a content encoder cannot see cross-file co-occurrence (single-interval input). The combination is a content-aware token embedding.

**This sharpens the predictions** in Step A (formerly "informative outcomes, not pass/fail"). Given LanceOtron's mechanism and Pundhir 2016's shape-beats-aggregate finding, we can now state directional predictions:

1. **F1_content > F1_R2V on cCRE class prediction (PLS / pELS / dELS).** Class is substantially determined by within-interval shape (promoter vs enhancer peak architecture, bimodal H3K4me3 flanking, mark deposition patterns) — exactly the axis R2V cannot see. Pundhir 2016 + LanceOtron directly support this prediction.
2. **F1_R2V > F1_content on tissue-activity prediction** (which tissues is this cCRE active in). Tissue identity is a cross-file co-occurrence property that R2V is structurally built to capture; shape is tissue-specific but the corpus-level signal carries identity more cleanly than within-interval shape alone.
3. **F1_combined > max(F1_R2V, F1_content)** on both tasks. This is the LanceOtron wide-and-deep result at embedding granularity — two orthogonal branches combine into a richer representation than either alone.

**This is where the genuine novelty lives.** LanceOtron's published work does not test tokenization-based embeddings at all; its aggregate branch is just enrichment-ratio scalars, not a distributional embedding. So comparing R2V-style co-occurrence against within-interval content is a contribution that LanceOtron does not pre-empt. By contrast, the within-cCRE shape pilot ([`./2026-04-19-within-ccre-shape-pilot.md`](./2026-04-19-within-ccre-shape-pilot.md)) is largely predicted by prior art; this pilot is where new ground is broken.

**Step B influence-attribution is untouched by the prior art**: Pundhir and LanceOtron do not attempt to characterize what an embedding method implicitly weights when aggregating to file-level. Step B's leave-one-out cCRE attribution remains novel.

**Tissue-activity task needed.** Step A currently focuses on cCRE class prediction. To test prediction #2 above, add a second classification task: **given a cCRE, predict which of K562 / GM12878 / HepG2 it is maximally active in** (or the binary active-in-tissue task per tissue). This should rely on the per-tissue content features that are already an optional extension in the current draft; promote them from optional to required under this repositioning. Use the ATAC activity thresholds from `01_select_ccres.py` as the tissue-activity ground truth.

**Revised status**: this pilot is now positioned as **the primary novelty-bearing experiment** in the workspace, with the within-cCRE shape pilot recast as an exercise/infrastructure-build that calibrates the content-extraction pipeline this pilot depends on.
