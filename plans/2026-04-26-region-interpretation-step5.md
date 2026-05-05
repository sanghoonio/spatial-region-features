---
date: 2026-04-26
status: draft
description: Plan for step 5 of the regulatory dictionary viz — how to interpret what a region is *for* or *doing* beyond the ENCODE/SCREEN annotations baked into the dictionary's entry cards. Open-ended; this plan scopes which interpretive moves are tractable for the course demo and which are research-paper directions.
---

# Step 5 — region interpretation beyond ENCODE annotations

## The question

After steps 1–4, the viz has shown a reader:

- BED files contain regions of various modalities (DNase, ATAC, ChIP).
- Files become tokenized as binary presence vectors against a fixed universe.
- The universe's regions cluster in embedding space — apparently meaningfully, since SCREEN classes separate visibly on the UMAP.
- A single BED file activates a subset of the universe; the activated subset has a shape in embedding space.

Step 5's question: **OK, but what does any specific region actually *do*?** Beyond "SCREEN says it's a dELS in K562," can the reader generate a useful hypothesis about a region's function?

This is the genuine research question of regulatory genomics interpretability. The grant's Challenge 2 ("interpreting genomic region embeddings") names it explicitly. We can't solve it; we can scope a demo move that demonstrates the *kind of question the dictionary helps us ask*.

## Why this is hard

1. **Annotations are incomplete.** SCREEN classes are coarse (5 broad bins). Many regulatory mechanisms (silencers, condition-specific enhancers, locus control regions, anti-sense regulation) aren't in standard catalogs.
2. **The embedding is bounded.** It captures co-occurrence patterns from the training corpus. Functional axes that aren't reflected in *which experiments contain which regions* won't show up.
3. **Interpretation is a generative act.** A region doesn't have a single "function" we can read off; it has effects in contexts. Asking "what is this region for?" is shorthand for "in what biological scenarios does this region's activity matter, and what does it modulate?"
4. **The mystery files we curated have UNKNOWN cell_line metadata.** We don't know what they are. The dictionary's only handle on them is their position in the file UMAP and their token activation.

## Candidate approaches

In rough order of feasibility for the 2-week demo:

### A. Hypothesis generation by neighborhood aggregation (most feasible)

For any region (especially the rare unclassed entries):
- Look up its k-NN in embedding space (already in `viz_chr16.parquet`).
- Aggregate the SCREEN classes / mark profiles of those neighbors.
- Surface a textual hypothesis: *"This region's nearest neighbors are 23/30 dELS, mostly active in K562 (high H3K27ac). Hypothesis: a K562-specific distal enhancer."*

Same trick for mystery files (step 4 highlights the regions they activate; step 5 asks what those regions are like *as a set*, suggesting what the file experiment might be).

**What we'd build**: a small Observable cell that does the kNN aggregation + builds a templated paragraph. Maybe a confidence indicator (homogeneity of the neighborhood).

**What's already in the data**: kNN graph (Stage 08 output), class labels, mark means.

**Risk**: kNN aggregations smooth toward the dominant class (dELS). The hypothesis text might be redundant with what Panel 1's class badge already says.

### B. Experiment-based context (medium feasibility)

For a region, enumerate which BED files activate it and look at the joint metadata of those files:
- "Region X is in the activated-token set of 142 of our 84,698 files. Of those, 89% are DNase-seq experiments, 67% are in tissue=blood. Hypothesis: a blood-specific accessibility region."

**What we'd build**: a per-region inverted index over the file corpus's tokenized vectors. Currently we don't have a tokenized version of the 84k corpus — we deliberately skipped that to use the pretrained embedding as-is. We'd need to tokenize a subset (say 1,000 files) to populate this.

**What's already in the data**: nothing.

**Risk**: re-introduces a chunk of the work we cut. Probably not worth it for the demo.

### C. Cross-modality validation (medium feasibility)

For a hypothesized region function from approach A or B:
- Cross-check with mark profile (do the marks match the expected pattern for this hypothesized function?).
- Cross-check with conservation (deeply conserved → suggests important function).
- Cross-check with sequence (check for the predicted TF motif in the actual DNA).

**What we'd build**: orchestration logic for the cross-checks. Conservation requires `phastCons100way` (Stage 05 work, deferred). Motif requires FIMO (Stage 05 work, deferred). Marks already in our parquet.

**What's already in the data**: mark profiles only.

**Risk**: requires Stage 05 work we deferred. Bumping it back in is feasible if the user prioritizes step 5 highly.

### D. Model probing on embeddings (research direction, not demo)

Train a small classifier on top of the pretrained embedding for a target task: "given embedding, predict tissue activity in K562," or "given embedding, predict SCREEN class." Use the classifier's confidence on novel/unclassed tokens as an interpretive output.

**What we'd build**: a sklearn `LogisticRegression` (or gradient-boosted tree) on the pretrained embeddings + a label of choice.

**Risk**: introduces an additional ML layer that's a methods-paper move, not a "dictionary" move. The reader would lose track of what the embedding alone already showed.

### E. Mystery-file walking (demo-friendly framing of A)

For the 20 mystery entries in `viz_files.parquet`:
- For each, find its file-level UMAP nearest neighbors among the labeled files.
- Aggregate those neighbors' assays / cell lines / tissues.
- Surface: *"Mystery entry X is 5 of its 10 file-UMAP neighbors are K562 DNase-seq, 3 are GM12878 ATAC-seq. Hypothesis: a primary blood cell accessibility experiment."*

Mirror of approach A but at the file level. Same code shape.

**What we'd build**: file-level kNN (we have UMAP coords; can compute kNN in Observable directly). Plus the templated text generation.

**What's already in the data**: file UMAP coords. (We don't currently have a file-level kNN graph; it's cheap in JS.)

## Recommended demo path

**Build (A) and (E) only** for the course MVP. Both reuse data we already have:

1. Region-side: for any selected region (especially unclassed ones), aggregate its 30 kNN's SCREEN classes + mark means → templated hypothesis.
2. File-side: for any selected mystery file, compute its file-UMAP kNN among labeled files → templated hypothesis.

These two cells (with simple text-generation logic) materialize step 5's hypothesis-generation hook. Reader sees: *"the dictionary lets us hypothesize about things SCREEN didn't catalog."*

Skip B, C, D for the MVP. They're real research directions; the demo doesn't need to fight that fight.

## What "good" looks like for the demo

A reader who lands on step 5 should:

- **Click a mystery file** → see "these 7 nearest neighbors are mostly DNase-seq from blood tissues, suggesting this is likely a blood DNase-seq experiment."
- **Click an unclassed region** (or any region without strong class signal) → see "its embedding-space neighbors are a mix of dELS (60%) and CA-CTCF (30%); marks show high CTCF binding in K562 → suggests an insulator-like distal element."
- **Realize the dictionary isn't just a lookup table** — it's a tool for asking "what is this thing like?" and getting a probabilistic-flavored answer grounded in the embedding's structure.

## What's NOT in scope (research-direction notes)

- Validation of the generated hypotheses (would require MPRA / actual functional assays).
- Confidence calibration of the hypothesis text (we'd need a held-out evaluation set).
- Multi-region hypothesis composition (regions in a window → joint hypothesis).
- Explicit attention-style attribution (which file in the corpus pushed this region toward this neighborhood?).
- Reverse: "give me a function description, find regions that match it" — text→region is grant Aim 2 territory.

## Open questions

1. **Should step 5 be its own panel, or a card-augmentation that fires on click?** A card augmentation reads as part of Panel 1's entry view; a separate panel reads as more deliberate. Probably the former for narrative continuity.
2. **How tight should the templated text be?** Pure data table is dry; full prose is over-claiming. Halfway: structured bullets ("dominant class," "dominant tissue," "best-guess function").
3. **Do we want any user-driven prompting?** *"What kind of regulatory element does this look like?"* with a free-text input that influences the aggregation? Probably not for the MVP — adds a chat-style affordance the demo doesn't need.

## Status

Draft. Implementation deferred to after Stage 10 lands and the Observable scaffold's main panels are working. Approach A + E are 1–2 hours of Observable cell work each.
