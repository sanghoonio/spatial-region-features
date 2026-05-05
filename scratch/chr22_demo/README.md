# profile-playground

A sandbox for getting hands-on with BED files and bigwigs — the file types that
regulatory genomics is actually stored in — and with the standard tools for
visualizing them and computing region-level profile statistics on them.

The data is pinned to **chr22** to keep everything fast and portable. That's
~1.6% of the genome, a few tens of MB on disk, and enough to make meta-profiles
and heatmaps look like they do in the papers.

## What's in here

```
profile-playground/
├── pyproject.toml               uv-managed Python project
├── data/
│   ├── hg38.chrom.sizes         UCSC-hosted hg38 chromosome sizes
│   ├── hg38_cCRE_chr22.bed      ENCODE cCRE registry, chr22 only (~43k regions)
│   ├── K562_ATAC_chr22.bw       K562 ATAC-seq signal p-value, chr22 only
│   ├── K562_H3K27ac_chr22.bw    K562 H3K27ac ChIP-seq signal p-value, chr22 only
│   ├── K562_H3K4me1_chr22.bw    K562 H3K4me1 ChIP-seq signal p-value, chr22 only
│   └── _full/                   original full-genome bigwigs (can be deleted
│                                once chr22 slices are built)
├── scripts/
│   ├── fetch_chr22_bigwigs.py   downloads ENCODE files + slices to chr22
│   ├── 01_pybigwig_basics.py    read signal, region summaries, multi-track stack
│   ├── 02_deeptools_profile.sh  computeMatrix + plotProfile + plotHeatmap
│   └── 03_gtars_bed_ops.py      gtars + bedtools interval operations on cCREs
├── igv/
│   └── playground_session.xml   IGV session auto-loading the chr22 tracks
└── outputs/                     written by the example scripts
```

## Data sources

All files are from ENCODE, GRCh38 assembly:

| File                         | Source                                              |
|------------------------------|-----------------------------------------------------|
| K562 ATAC-seq                | ENCSR483RKN / ENCFF867QEW (signal p-value bigwig)   |
| K562 H3K27ac ChIP-seq        | ENCSR000AKP / ENCFF467OGB (signal p-value bigwig)   |
| K562 H3K4me1 ChIP-seq        | ENCSR000EWC / ENCFF948OCT (signal p-value bigwig)   |
| hg38 cCRE registry (agnostic)| ENCSR800VNX / ENCFF420VPZ                           |
| hg38 chromosome sizes        | hgdownload.soe.ucsc.edu/goldenPath/hg38             |

The K562 tracks give you the canonical **active enhancer triplet**: chromatin
accessibility (ATAC), the active-mark histone modification (H3K27ac), and the
general enhancer mark (H3K4me1). Running a meta-profile around cCRE centers
should reproduce the classic figure-one pattern: ATAC and H3K27ac peak sharply
at the center, H3K4me1 shows a bimodal shape flanking it.

## Setup

Environment is managed at the **workspace root** (`~/Documents/Work/spatial-region-features/`), not here. One `uv sync` at the root installs everything the scratch scripts need alongside the plan experiments' deps.

```bash
cd ~/Documents/Work/spatial-region-features
uv sync                       # once, at root
brew install bedtools         # once, for deeptools interop
```

Scripts in this dir use **relative paths** (e.g. `data/...`, `outputs/...`), so run them from inside `scratch/chr22_demo/`; `uv` walks up to the root `pyproject.toml` and uses the root env automatically:

```bash
cd scratch/chr22_demo
uv run python scripts/fetch_chr22_bigwigs.py
```

IGV 2.19.7 can be installed separately at `/Applications/IGV_2.19.7.app` — not managed by uv.

## Running the examples

### 1. pyBigWig basics — what is a bigwig, and how do you read it?

```bash
uv run python scripts/01_pybigwig_basics.py
```

Walks through five things in order:

1. Opening a bigwig and reading its header (chromosomes, lengths).
2. Reading per-base signal values in a narrow window (`bw.values(...)`).
3. Region-level summaries with `bw.stats(chrom, start, end, type=...)`.
4. The **row view** — all three tracks queried at a single position, giving
   you the "what does the regulatory stack look like here" picture.
5. Scoring every chr22 cCRE against K562 H3K27ac and ranking the top ten.

This is the conceptual foundation: the bigwig is a continuous signal value at
every genomic position, the BED file is a set of regions, and the primitive
you care about is "signal-at-regions" as a summary statistic.

### 2. deepTools meta-profile and heatmap

```bash
./scripts/02_deeptools_profile.sh
# or:  uv run bash scripts/02_deeptools_profile.sh
```

Runs the three canonical deepTools commands:

1. `computeMatrix reference-point` — scores every cCRE against all three
   bigwigs in +/- 2 kb windows at 50 bp bin resolution. Produces a 3D tensor
   of shape `(regions × bigwigs × bins)` saved as a gzipped matrix.
2. `plotProfile` — averages across regions, produces the meta-profile plot
   (signal as a function of position around cCRE centers).
3. `plotHeatmap` — plots the full tensor as a heatmap with one row per cCRE,
   sorted by H3K27ac signal.

Outputs:

- `outputs/ccre_matrix.gz`
- `outputs/ccre_profile.png`
- `outputs/ccre_heatmap.png`

Open the PNGs. The meta-profile should show ATAC and H3K27ac peaking sharply
at the cCRE center; H3K4me1 should be bimodal (peaks flanking the center, dip
at the center itself — because H3K4me1 marks the nucleosomes flanking active
regulatory elements but not the accessible core).

This is the conceptual template for the Phase A1 linkage tool. `computeMatrix`
takes a BED file and a set of bigwigs, produces a tensor. Phase A1's variant-
centered version would take a VCF/variant list instead of a BED file and
produce the same kind of tensor, keyed on variant positions with windowing
handled per variant.

### 3. gtars + bedtools interval operations

```bash
uv run python scripts/03_gtars_bed_ops.py
```

Walks through:

1. What gtars exposes as its top-level API.
2. Counting cCREs by SCREEN classification (pELS / dELS / PLS / CA-H3K4me3 /
   CA-CTCF / CA-TF / CA / TF).
3. Running `bedtools merge -d 1000` to collapse proximal-enhancer-like
   cCREs into merged clusters, then counting the before/after.
4. Scoring the merged clusters against the K562 H3K27ac bigwig with
   pyBigWig and listing the top 10 by mean signal.

The point of this script is to show the primitives layer — how you actually
call into gtars / bedtools / pyBigWig to do the small operations you'd
compose into a larger pipeline. It's deliberately not a user-facing analysis;
it's the plumbing.

## Visualizing in IGV

The data is designed to be opened in IGV with one step. Start IGV:

```bash
open -a IGV_2.19.7
```

Then load the session file via **File → Open Session...** and pick
`igv/playground_session.xml`. The three bigwig tracks and the cCRE BED will
load automatically, zoomed to `chr22:30,000,000-32,000,000` (a region with
plenty of active regulation to look at).

**Things worth doing once it's loaded:**

- Scroll around chr22 and watch how the three tracks co-vary. Peaks where
  ATAC + H3K27ac are both high are active enhancers or promoters; peaks where
  only H3K4me1 is high (with weaker H3K27ac) are poised enhancers.
- Zoom into a single peak and look at the shape. ATAC peaks are usually
  narrower and sharper than histone ChIP peaks because DNase/Tn5 sensitivity
  localizes to the ~150 bp of unprotected DNA between flanking nucleosomes.
- Right-click on a track → **Set Data Range** and flip between autoscale and
  manual scale to see what the per-track signal actually looks like vs what
  IGV is defaulting to.
- Compare a few cCREs flagged by `01_pybigwig_basics.py` as top-ranked for
  H3K27ac: jump to their coordinates via the locus bar and see whether the
  visual signal matches the ranking.

## What this playground is not

- Not a faithful analysis pipeline. The bigwigs are signal p-value tracks
  from a single replicate each; real analyses use IDR-thresholded peaks,
  multiple replicates, and proper input controls.
- Not a substitute for learning the underlying biology. A meta-profile that
  looks "right" in a heatmap doesn't tell you anything without context about
  what the regulatory elements are, what cell type they're active in, and
  what the assay is actually measuring.
- Not a complete tour. There's no Hi-C / contact data, no TF ChIP, no RNA-seq,
  no methylation. If you want to branch out, the same
  `scripts/fetch_chr22_bigwigs.py` pattern extends trivially to any other
  ENCODE experiment.

## Cleanup

To reclaim disk space once you're done:

```bash
# remove the full-genome bigwigs (~1.5 GB) but keep the chr22 slices
rm -rf data/_full
```

To nuke everything and start over:

```bash
rm -rf ~/Documents/Work/profile-playground
```
