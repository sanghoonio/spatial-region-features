# spatial-region-features

Research workspace for exploring spatial patterns **inside** known regulatory-element intervals — the universe-as-scope methodology where a fixed universe (e.g., ENCODE cCREs, tissue-matched consensus peaks) defines *where* the model looks, while the representation inside each interval stays continuous and fine-resolution.

A "universe," in this workspace's Sheffield-Lab inherited sense, is a fixed collection of genomic regions that bounds what any downstream analysis can observe. Under **universe-as-scope**, the universe's role is to *scope* the model's domain of computation — not to tokenize and collapse intra-region structure as in Atacformer-style discrete-token foundation models.

## The motivating question

Token-based embedding models for genomic region sets (Region2Vec, scEmbed, Atacformer) produce per-region vectors whose similarity comes from **co-occurrence** across the training corpus. They are useful for retrieval, clustering, and set-level classification, but they do not tell a user *what shared tokens actually represent*. "Your query BED shares 5000 tokens with these 30 BEDs" is biologically opaque — tokens are coordinates, and their similarity is defined relationally (by what co-occurs) rather than semantically (by what they contain).

This workspace investigates whether learning spatial features **inside** known regulatory-element intervals can give universe tokens transferable semantic types, making existing embedding-based region tools interpretable. The near-term focus is one concrete experiment at a time — e.g., whether within-cCRE spatial features cluster by cCRE class reproducibly across tissues — not a multi-phase research program.

## Why universe-as-scope

Alternatives, and why they are not what this workspace does:

- **Universe-as-vocabulary** (Atacformer, ChromFound, EpiAgent): discrete tokens, one vector per universe region; intra-region structure collapsed.
- **No universe** (Enformer, AlphaGenome): genome-wide 128-bp bins; compute cost beyond a single-student budget.
- **Universe-as-scope** (this workspace): the universe tells the model where to compute features; inside each interval, representation stays continuous and multi-resolution.

Two justifications:

1. **Compute tractability.** ~1 M intervals of ~2–10 kb is a different compute regime than genome-wide modeling at full resolution.
2. **Biological interpretability.** Features anchored to known regulatory elements aggregate into statements with biological identity ("mean summit sharpness across PLS cCREs"); features over bottom-up groupings do not. This is the lesson from the retired `region-cnn` and `bed-embedding-eval` predecessors — see their READMEs.

## Relationship to the dissertation workspace

The sibling `../dissertation/` directory contains framings, question drafts, and literature notes that inform this workspace. The framings were written against a research-program-scale narrative; the empirical work here is scoped tighter. Read in order:

- `../dissertation/framing/universe-as-scope.md` — the methodological commitment.
- `../dissertation/framing/embedding-interpretability.md` — the downstream motivation.
- `../dissertation/framing/dissertation-argument.md` — aspirational through-line. (Caveat: written in the language of a committed PhD dissertation; the current reality is pre-PhD exploration.)
- `../dissertation/questions/q1-profile-predictivity.md` — the broader question this workspace's experiments contribute to.

## Current contents

- `plans/` — experiment plans, one file per planned experiment. Currently: within-cCRE shape pilot and R2V per-token content comparison.
- `retired/` — retrospectives on predecessor projects that were walked away from.
  - `retired/region-cnn/` — 1D-CNN spatial features without a universe commitment. See its README for why retired.
  - `retired/bed-embedding-eval/` — evaluation of R2V / text2bed / Atacformer embeddings via file-level geometric correlation. See its README for why retired.
- `scratch/` — a small hands-on sandbox with ENCODE chr22 bigwigs (K562 ATAC / H3K27ac / H3K4me1), the SCREEN cCRE chr22 slice, and example scripts for `pyBigWig`, `deeptools`, and `gtars`. Useful for prototyping feature extraction before scaling up to the plans. Not tracked content-wise for science; for learning-by-doing.

New experimental subdirectories will be added as specific questions get scoped and run.

## Setup

One `uv` environment for the whole workspace. Scratch scripts and future experiment subdirectories share it.

```bash
cd ~/Documents/Work/spatial-region-features
uv sync
```

`pyproject.toml` at the root declares the full dependency set (bigwig / region handling for scratch, plus numpy / scikit-learn / matplotlib / umap-learn etc. for the plans). `.venv/` is gitignored.

Python 3.13+. Optional system deps: `bedtools` (for deeptools interop in scratch), IGV (for visual inspection of bigwigs).

## Status

Exploratory. The workspace exists to help decide — through concrete empirical results on small, well-scoped experiments — whether within-universe spatial-feature learning is interesting enough to build a longer-term research program around. No experiments are underway yet.
