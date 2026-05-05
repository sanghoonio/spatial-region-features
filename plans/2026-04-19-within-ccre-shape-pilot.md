---
date: 2026-04-19
status: draft
description: Pilot experiment testing whether within-cCRE ATAC peak-shape features distinguish regulatory classes (PLS / pELS / dELS) reproducibly across three tissues, using hand-engineered features on continuous bigwig signal.
---

# Within-cCRE peak-shape pilot

**Goal:** produce a single concrete result that either motivates continued work on within-universe spatial features or cleanly closes the direction.

## The question

**Do within-cCRE ATAC peak-shape features distinguish SCREEN cCRE classes (PLS, pELS, dELS), reproducibly across three tissues, in a way that is not trivially reproduced by total signal alone?**

Alternative framings considered and set aside for now:

- *"Within-cCRE pattern predicts validation status (MPRA hit vs non-hit)."* More ambitious; depends on MPRA label sources and belongs to a second-stage experiment that this pilot is the falsifier for.
- *"Within-cCRE pattern distinguishes tissues within a single cCRE class."* Interesting but tests a weaker property — cross-class separation is the cleaner sanity check first.
- *"Full continuous multi-track spatial encoder (learned features)."* That is the Phase D item in the Q1 roadmap. This pilot is the gate that has to pass before investing in the learned version.

## Why this is the right first experiment

1. **Bounded scope.** Three tissues, one assay, three cCRE classes. All data public in ENCODE. Runs on a laptop or a small SLURM job.
2. **Direct test of the universe-as-scope commitment.** If within-cCRE spatial patterns do not distinguish cCRE classes even with hand-features, the commitment has a problem that the bigger program has to answer. If they do, universe-as-scope earns empirical support independent of any learned-encoder story.
3. **Avoids the circularity that bit the retired `region-cnn`.** All features are computed inside biology-anchored intervals (cCREs have known regulatory class labels). Aggregations ("mean summit sharpness across PLS cCREs in K562") have biological identity and are interpretable. Compare the retired region-cnn's feature aggregations over stitched-but-unlabeled groupings.
4. **Also partially avoids the SCREEN-definitional circularity.** SCREEN classifies cCREs using DNase + H3K4me3 signal. This pilot uses **ATAC** as the feature source. ATAC and DNase both measure accessibility so the residual overlap is real, but ATAC vs DNase is not the same signal — a pilot result using ATAC features is a real test rather than a tautology. (A DNase arm is proposed below as a positive control.)
5. **Produces the decision input we actually need.** The question informing this workspace right now is whether within-universe spatial-feature learning is exciting enough to build a longer-term program around. A clear yes / moderate / no on "within-cCRE shape features carry class signal" is that kind of result.

## Method

### Data

- **Universe**: ENCODE SCREEN cCRE registry (hg38, current stable version). Each interval has a class label among PLS, pELS, dELS, CTCF-only, DNase-H3K4me3. The plan uses the **biosample-agnostic "core" classification** (integrated signal across all biosamples) so that every cCRE has a single stable class label independent of tissue. SCREEN also publishes biosample-specific class labels where the same cCRE can be PLS in one tissue and pELS in another — those are out of scope for this pilot; a follow-up could test tissue-specific-class reproducibility using them.
- **Tissues**: three well-characterized ENCODE cell lines with deep ATAC coverage — **K562, GM12878, HepG2**.
- **Signal source**: ENCODE ATAC-seq **bigwig** files (fold-change-over-control track), not peak-call BED files. The question is about within-interval shape, which requires continuous signal.
- **Scope reduction**: for the pilot, restrict to PLS, pELS, dELS (three classes). Skip CTCF-only and DNase-H3K4me3 for simplicity and class-size balance.
- **Replicate handling**: use the ENCODE primary replicate per tissue. No cross-replicate averaging for the pilot.

### Feature extraction

For each (tissue, cCRE) pair, extract the ATAC fold-change signal vector over the cCRE interval. Compute a small hand-engineered feature set per interval:

1. **Summit position** — bp offset of the argmax position from the interval center, normalized by interval length.
2. **Summit sharpness** — FWHM of the summit peak within the interval, normalized by interval length.
3. **Signal concentration** — fraction of total signal in the top 20% of positions by value.
4. **Edge-to-core ratio** — mean signal in the outer 25% vs the inner 50% of the interval.
5. **Total signal** — sum of signal across the interval. Density-equivalent. Separated from the shape features for the ablation (see below).

Total: 4 shape features + 1 density baseline, per (tissue, cCRE) pair.

### Analysis

1. **Within-tissue class separation.** For each tissue separately, train a multinomial logistic regression to predict cCRE class from features. Report macro-F1 across the three classes. Run three model variants: (a) all 5 features, (b) shape features only (drop total signal), (c) total signal only (density baseline). The **shape contribution** is F1(a) − F1(c); reporting F1(b) alongside separates "shape adds to density" from "shape alone is sufficient."

2. **Cross-tissue reproducibility.** Train on two tissues, test on the third, for each held-out choice. If shape features transfer, cross-tissue F1 is within a reasonable margin of within-tissue F1; if the features are tissue-idiosyncratic, cross-tissue F1 degrades.

3. **Density baseline.** Model variant (c) above — classifier trained on total-signal-only. Compare to variants (a) and (b) to establish how much of the class signal is "this cCRE has high ATAC" vs "this cCRE has a specific ATAC shape."

4. **DNase positive control (optional, run if time permits).** Repeat the full pipeline with ENCODE DNase-seq signal for the same three tissues. DNase was used in SCREEN's cCRE classification, so DNase features should recover cCRE class near-trivially. Establishes the upper bound of what "signal shape distinguishes cCRE class" can look like when the signal is near-definitional. The interesting comparison is how much of that upper bound ATAC recovers.

### Success criteria

- **Strong positive**: macro-F1 > 0.70 with ATAC shape features, cross-tissue F1 within 0.05 of within-tissue F1, shape contribution > 0.10 F1 above density-only.
- **Moderate positive**: macro-F1 0.55–0.70, shape contribution 0.05–0.10 — some real signal, but bounded.
- **Negative**: macro-F1 < 0.55, OR shape contribution < 0.03 F1 — shape features as defined do not carry class-relevant signal. Either the features are too coarse (motivating a learned encoder) or within-cCRE shape does not carry much class-relevant information that density does not already carry.

### What each outcome teaches

- **Strong positive**: within-cCRE shape features carry class-relevant, tissue-transferable information. Universe-as-scope is empirically defensible. Expand: more tissues, more assays, multi-track joint features, eventually a learned encoder.
- **Moderate positive**: some signal, bounded. The ceiling question becomes the driving one — what's missing that hand-features can't capture? Motivates a learned-encoder pilot as the next step.
- **Negative**: the within-cCRE spatial story is weaker than expected at the hand-feature level. Either (a) the feature set is too coarse (worth a learned-encoder retry), (b) the cCRE classes are already fully captured by interval identity under the signal source used (worth testing with a different assay, e.g., H3K27ac), or (c) within-cCRE shape genuinely does not carry more class information than total signal. Each is a concrete update to the program.

## Tractability

- **Data**: 3 tissues × (1 ATAC bigwig per tissue) = 3 files, ~1 GB each. Optional DNase arm doubles this. ENCODE direct download, ~30 min.
- **Compute**: feature extraction is one-pass over ~1 M cCREs × 3 tissues. With `pyBigWig` (or an equivalent bigwig accessor), ~15–30 min per tissue on a laptop. Fits on a laptop; no SLURM needed for the pilot.
- **Analysis**: classifier training on ~3 × 1 M rows × 5 features is trivial via `sklearn` (seconds).
- **Timeline**: ~1–2 weeks from start to a first decisive result, assuming no data-curation surprises.

## Open design decisions before executing

1. **ATAC vs DNase choice.** ATAC as primary, DNase as positive-control arm if time permits. Alternative: run both from the start to make the circularity caveat a clean comparison rather than a disclaimer.

2. **Bigwig representation.** ENCODE publishes fold-change-over-control, pileup, and signal p-value bigwigs. Propose **fold-change-over-control** as the standard published normalization. Revisit if the features look noisy.

3. **Interval-length handling.** cCREs vary (typically 150–350 bp). Two options: (a) compute features on variable-length intervals directly, using length-normalized features where meaningful (as drafted above); (b) zero-pad every cCRE to a fixed 1 kb window centered on the interval, compute features on the padded window. Option (a) is cleaner when the features are length-invariant; option (b) is needed if we later want to feed the signal vectors into a learned encoder at fixed length. Start with (a) for the hand-feature pilot.

4. **Feature set validation.** The 5 hand-features are a reasonable sketch, not a principled biology-down set. Worth a 30-minute lit search on *"spatial signatures of promoters vs enhancers in chromatin accessibility data"* before executing — maybe there's a known-to-discriminate feature we are missing (e.g., nucleosome footprint V-plots, bimodal-vs-unimodal distributions around TSS, directional asymmetry around the summit).

## Relationship to the broader workspace

- Tests the **universe-as-scope** commitment articulated in [`../../dissertation/framing/universe-as-scope.md`](../../dissertation/framing/universe-as-scope.md). A positive result is empirical support for the methodological commitment independent of any learned-encoder work; a negative result forces the commitment to be defended on other grounds.
- Contributes to Q1 cluster (F), sub-question F1 ("do within-token multi-track pattern clusters align with known regulatory categories"), in a simplified single-track, hand-feature, single-classifier form. See [`../../dissertation/questions/q1-profile-predictivity.md`](../../dissertation/questions/q1-profile-predictivity.md).
- Explicitly **not** the Phase D self-supervised encoder from the Q1 buildup roadmap. This is the falsifier that has to pass first. If it fails, the Phase D design has to change (or the whole program direction does). If it passes, Phase D is better justified.

## Status and next steps

**Draft.** Before executing:

1. Sam to confirm the three design decisions above (ATAC-only vs ATAC+DNase, bigwig choice, variable-length vs fixed-window).
2. Sam to confirm the feature set, either as drafted or after the literature check.
3. Scaffold a minimal project directory (e.g., `within-ccre-shape/`) under `spatial-region-features/` with `src/`, `data/`, `results/` once the plan is frozen.

First implementation step after the plan is frozen: download the three ATAC bigwigs and the cCRE registry, verify the universe-to-tissue alignment works cleanly on a small subset (e.g., chr22 only), then scale to genome-wide feature extraction.

## 2026-04-19 update — prior art and repositioning

After a literature check (see [`../../dissertation/background/within-peak-shape.md`](../../dissertation/background/within-peak-shape.md)) the motivation for this pilot has shifted. Two pieces of prior art substantially overlap the question:

- **Pundhir et al. 2016** (PARE, Nucleic Acids Research) showed that within-peak shape (valley depth `nfrDip`, peak-valley-peak asymmetry) on H3K4me1/H3K4me3 discriminates active enhancer vs active promoter and active vs poised elements at 0.89 specificity vs DEEP's 0.52 on 105/81 validated HeLa enhancers. Narrow marks + narrow task + narrow universe, but a direct "shape beats aggregate" precedent.
- **Hentges et al. 2022** (LanceOtron, Bioinformatics) trains a "wide and deep" model that explicitly separates a shape branch (CNN on 2 kb bigwig coverage) from an aggregate enrichment branch (logreg on 11 enrichment-ratio scalars at multiple window sizes) and fuses them via an MLP. Beats MACS2 on F1 across ATAC, DNase, TF ChIP, histone ChIP. Task is peak-vs-noise, not cCRE class, but the architecture is directly transferable.
- **Vu et al. 2023** (RCL, Genome Research) quote: "mean coverage as a simple score already does well... but RCL learns additional signals, perhaps peak shape, that further improve the performance."

**Repositioning**: the pilot's motivation is no longer "does shape carry signal beyond aggregates" — Pundhir and LanceOtron have settled the yes at narrow scope. The sharper residual claim is "does shape carry signal across the full SCREEN 5-way cCRE taxonomy, with multi-track inputs, at registry scale, against a modern max-Z baseline." That is genuinely untested. But the answer is now strongly predicted (yes), which shrinks the novelty of this pilot considerably.

**Practical consequence**: the pilot is more usefully run as an **exercise + infrastructure build** than as a standalone paper contribution. Keep the scope tight, borrow the LanceOtron wide+deep ablation template explicitly (shape features vs aggregate features vs combined), and treat the run primarily as calibration for the R2V comparison pilot ([`./2026-04-19-r2v-token-content-comparison.md`](./2026-04-19-r2v-token-content-comparison.md)), which addresses the part LanceOtron does not touch (tokenization-based embeddings).

**Additions to the feature set**. The current 5 features (summit position, summit sharpness, signal concentration, edge-to-core ratio, total signal) should be augmented with Pundhir's validated features:

- **`nfrDip`** — valley depth between the two flanking summits (the PVP valley). The current `edge-to-core ratio` is a coarser cousin; run both.
- **PVP asymmetry** — ratio of upstream to downstream flanking-summit heights. Only meaningful for strand-aware tracks (H3K4me3 mostly). Captures TSS-directionality information that Pundhir showed correlates with GRO-seq (ρ=0.45) and CAGE (ρ=0.63).

Add these as a baseline arm, since skipping them leaves the obvious reviewer question ("did your feature set include the known-good features?") unanswered.

**BED-vs-bigwig arm justified**. The BED-vs-bigwig comparison we added (pilot of what binary peak input captures vs continuous signal) is empirically predicted by Vu 2023's finding that learned shape adds signal above mean coverage on the same interval. Keep this arm; it now has an adjacent precedent to cite.

**Remaining genuine contributions** (what is not-yet-settled):

- Multi-track shape features (jointly across ATAC + H3K27ac + H3K4me1 + H3K4me3 + H3K27me3 + DNase), not the single-track that Pundhir and LanceOtron tested.
- 5-way SCREEN classification including CTCF-only and DNase-H3K4me3, not the binary/ternary that Pundhir tested.
- Full registry scale (~2.3 M rDHSs), not the ~200 hand-validated intervals Pundhir used.
- Tissue-matching decomposition — does shape predict class consistently across tissues, or is the shape signature tissue-dependent?
