# ccre_examples

Curated playground for the planned pilots. Selects a small set of cCREs
stratified by class + a biologically meaningful property + validation status,
then pulls bigwigs and peak BEDs for six modalities across three tissues so
you can visually compare continuous-signal and binary-peak representations
inside real regulatory intervals before scaling up.

## What gets selected (30 cCREs total)

| class | stratum                     | n | source of label                       |
|-------|-----------------------------|---|---------------------------------------|
| PLS   | ubiquitous (3 tissues)      | 5 | ATAC activity in K562+GM12878+HepG2   |
| PLS   | K562-specific               | 5 | ATAC activity in K562 only            |
| pELS  | ubiquitous                  | 5 | ATAC activity in K562+GM12878+HepG2   |
| pELS  | K562-specific               | 5 | ATAC activity in K562 only            |
| dELS  | Fulco-2019 positive         | 5 | CRISPRi-validated enhancer in K562    |
| dELS  | Fulco-2019 negative         | 5 | tested in CRISPRi, not significant    |

Tissue-specificity is the biologically meaningful property we stratify by for
PLS and pELS. For dELS we swap in a validation contrast — both sides were
tested in the same CRISPRi screen, one scored as functional, the other did
not. This is cleaner than "validated vs. untested" because both sides have
actual experimental evidence.

## Tissues and modalities

- **Tissues**: K562, GM12878, HepG2 (matches the within-cCRE shape pilot).
- **Modalities per tissue** (each as both bigwig + peak BED):
  - ATAC (pilot primary)
  - DNase (pilot positive control; near-definitional for cCRE class)
  - H3K27ac (active mark)
  - H3K4me3 (promoter mark — should light up on PLS)
  - H3K4me1 (enhancer mark — should light up on pELS/dELS)
  - H3K27me3 (repressive mark — negative contrast)

Six modalities × three tissues × (bigwig + peak BED) = 36 source files. Pulled
once, sliced per cCRE + 2 kb flank, then the full files can be deleted.

## Layout

```
ccre_examples/
├── README.md
├── scripts/
│   ├── 01_select_ccres.py        # select 30 cCREs; writes data/ccres.bed + metadata
│   ├── 02_fetch_and_slice.py     # (planned) download 36 files, slice per cCRE
│   └── 03_plot_ccre.py           # (planned) per-cCRE track-stack PNG
├── data/                         # gitignored
│   ├── _cache/                   # downloaded ENCODE files (registry, bigwigs)
│   ├── fulco_2019_crispri.tsv    # user-provided (see below)
│   ├── ccres.bed                 # 30 selected intervals, BED6
│   ├── ccres_metadata.csv        # class, stratum, validation, per-tissue activity
│   ├── signal/                   # per-cCRE signal slices (from 02)
│   └── peaks/                    # per-cCRE peak overlaps (from 02)
└── outputs/                      # gitignored
    └── plots/                    # per-cCRE PNG stacks (from 03)
```

## Setup

Workspace-root `uv` env (same as `chr22_demo/`):

```bash
cd ~/Documents/Work/spatial-region-features
uv sync
```

## Fulco 2019 CRISPRi table (auto-fetched)

Fulco et al. 2019 (Nature Genetics, doi:10.1038/s41588-019-0538-0) published a
CRISPRi screen of ~5k enhancer-gene pairs in K562. `01_select_ccres.py` pulls
a pre-lifted, pre-harmonized copy from the Engreitz lab's ENCODE-E2G benchmark
repo (Gschwind 2023), filters to the Fulco 2019 subset (3,501 pairs on
GRCh38), and writes `data/fulco_2019_crispri.tsv` with columns `chrom, start,
end, gene, significant`.

The `significant` column is the benchmark's `Regulated` flag (significant +
negative effect on target gene — i.e., the canonical enhancer-regulation
definition, not just "any significant effect"). Of the 3,501 Fulco 2019 pairs
in the benchmark, 52 are `Regulated=True`.

**Numerical caveat**: the benchmark re-ran statistics through a common
pipeline, so effect sizes and p-values differ from the original paper
(hit/miss calls are preserved). If exact original values matter, swap the
benchmark TSV for Supplementary Table 6a from the Springer XLSX
(`https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0538-0/MediaObjects/41588_2019_538_MOESM3_ESM.xlsx`, sheet `Supplementary Table 6a`,
hg19 — requires liftOver).

If the auto-fetch fails (e.g. offline, GitHub down), the script prints a
warning and falls back to stratifying dELS by tissue-specificity (same design
as PLS/pELS), marking `validation_status = "untested"`.

## Running

```bash
cd scratch/ccre_examples
uv run python scripts/01_select_ccres.py
```

First run downloads:
- ENCODE cCRE registry (ENCFF420VPZ, ~18 MB gzipped)
- Three ATAC-seq fold-change-over-control bigwigs (~300-500 MB each)

Total ~1.5 GB to `data/_cache/`. Subsequent runs use the cache.

Outputs `data/ccres.bed` (30 rows, BED6) and `data/ccres_metadata.csv`
(30 rows, full per-cCRE metadata).

## Design notes

- **Tissue activity is a local approximation of SCREEN's biosample-specific
  labels.** SCREEN publishes per-biosample cCRE classifications, but they are
  not exposed as a simple static download. We approximate by thresholding
  ATAC fold-change-over-control at the 80th percentile per tissue. This is a
  tunable in `01_select_ccres.py` (`ACTIVITY_PERCENTILE`).

- **Fulco 2019 is K562-specific.** The validation stratum only applies to
  K562 dELS cCREs. We could extend to HepG2 with Agarwal et al. 2025 MPRA or
  to other tissues with more recent CRISPRi screens; not in this first pass.

- **PLS may under-populate K562-specific.** Promoters are usually ubiquitous
  across cell types, so the K562-specific PLS pool may be small. The script
  samples whatever is available and logs the count.

- **ENCODE experiment choices**: the script uses canonical ATAC-seq experiments
  for each tissue. IDs are hard-coded at the top of `01_select_ccres.py` and
  easy to change. The fold-change-over-control bigwig is resolved at runtime
  via the ENCODE REST API.

## Playground observations (2026-04-19)

From eyeballing the IGV sessions across the 30 cCREs:

- **ATAC dominates visually at most cCREs, and aligns tightly with the cCRE
  interval.** Partly a rendering effect (autoscale gives ATAC its own y-axis,
  and ATAC FC peaks are much larger than histone FC). Partly mechanism: ATAC
  concentrates at the ~200 bp nucleosome-depleted core while histone marks
  smear across flanking nucleosomes. Partly definitional: SCREEN built cCREs
  by anchoring on DNase summits, so any accessibility signal is centered on
  the interval by construction. Implication: **ATAC features at cCRE scale
  benefit from a shared-signal-source advantage** over H3K27ac / H3K4me*
  features when predicting SCREEN cCRE class — the H3K27ac / H3K4me*
  features are the less-circular test for the shape pilot.

- **BED peaks omit most sub-peak structure visible in the bigwigs.** Much of
  the shape seen in the continuous tracks never shows up as called peaks
  unless the signal is very strong. Three causes stack: peak-caller
  thresholding (MACS2 / IDR / replicated-peaks drop sub-q-value shoulders
  even when biologically obvious), rounded boundaries (MACS2 extends by
  fragment length and bins to ~100–200 bp), and magnitude-flattening (BED
  binarizes so FC=3 and FC=50 look identical). Implication:
  **tokenization-based models (R2V, scEmbed, Atacformer) inherit these three
  losses by construction** — their per-token embeddings cannot encode
  anything the peak-caller discarded. This is the direct empirical case for
  the universe-as-scope + continuous-inside-intervals commitment. Caveat:
  not all sub-threshold shape is real signal; some is noise. Separating the
  two is part of why a learned encoder over continuous signal might
  eventually earn its weight over hand features.
