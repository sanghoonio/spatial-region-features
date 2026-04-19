#!/usr/bin/env bash
#
# deepTools meta-profile and heatmap on K562 cCREs.
#
# Runs the canonical "compute a (regions x bigwigs x bins) matrix then plot it"
# pipeline. Produces:
#
#   outputs/ccre_matrix.gz      - the per-region per-bin signal matrix
#   outputs/ccre_profile.png    - meta-profile plot (average signal vs position)
#   outputs/ccre_heatmap.png    - heatmap with one row per cCRE
#
# Run from the playground root:
#   ./scripts/02_deeptools_profile.sh
# or under uv:
#   uv run bash scripts/02_deeptools_profile.sh

set -euo pipefail

# Resolve the playground root regardless of cwd
PLAYGROUND="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PLAYGROUND"

BIGWIGS=(
    "data/K562_ATAC_chr22.bw"
    "data/K562_H3K27ac_chr22.bw"
    "data/K562_H3K4me1_chr22.bw"
)
LABELS=("ATAC" "H3K27ac" "H3K4me1")
REGIONS="data/hg38_cCRE_chr22.bed"
MATRIX="outputs/ccre_matrix.gz"
PROFILE="outputs/ccre_profile.png"
HEATMAP="outputs/ccre_heatmap.png"

mkdir -p outputs

echo "[1/3] computeMatrix: scoring $REGIONS against ${#BIGWIGS[@]} bigwigs"
echo "      window: +/- 2 kb around cCRE centers, 50 bp bins"
computeMatrix reference-point \
    --referencePoint center \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 2000 \
    --binSize 50 \
    --scoreFileName "${BIGWIGS[@]}" \
    --samplesLabel "${LABELS[@]}" \
    --regionsFileName "$REGIONS" \
    --outFileName "$MATRIX" \
    --skipZeros \
    --numberOfProcessors 4
echo "      wrote $MATRIX"

echo
echo "[2/3] plotProfile: meta-profile (average signal vs position around center)"
plotProfile \
    --matrixFile "$MATRIX" \
    --outFileName "$PROFILE" \
    --perGroup \
    --plotTitle "K562 cCREs on chr22 (n=$(wc -l < "$REGIONS" | tr -d ' '))"
echo "      wrote $PROFILE"

echo
echo "[3/3] plotHeatmap: one row per cCRE, columns are bins, color = signal"
plotHeatmap \
    --matrixFile "$MATRIX" \
    --outFileName "$HEATMAP" \
    --colorMap Blues Reds Greens \
    --whatToShow 'heatmap and colorbar' \
    --sortRegions descend \
    --sortUsingSamples 2
echo "      wrote $HEATMAP"

echo
echo "Done. Open outputs/ccre_profile.png and outputs/ccre_heatmap.png."
echo
echo "What you're seeing:"
echo "  - profile: average signal across all cCREs as a function of position"
echo "    relative to the cCRE center. ATAC and H3K27ac should peak at the"
echo "    center; H3K4me1 should show a bimodal shape flanking the center."
echo "  - heatmap: each row is one cCRE, each column is one 50 bp bin, color"
echo "    encodes signal. Rows are sorted by H3K27ac signal (descending), so"
echo "    strong enhancer-like elements cluster at the top of the plot."
