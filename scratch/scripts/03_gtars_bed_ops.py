"""
gtars basics: interval operations on BED files.

Demonstrates the Sheffield Lab's Rust-backed gtars package for region-set
operations. gtars sits at the 'primitives layer' — it's the library you'd
import to build a pipeline on top of, rather than a user-facing analysis tool.

The script loads the chr22 ENCODE cCRE BED into a gtars RegionSet, computes
a few basic statistics, filters to proximal-enhancer-like elements (pELS),
builds a new RegionSet from the filtered rows, and scores the filtered
regions with K562 H3K27ac signal via pyBigWig.

Outputs written to outputs/ for downstream inspection.
"""

from pathlib import Path
from statistics import median
import pyBigWig
from gtars.models import RegionSet

PLAYGROUND = Path(__file__).resolve().parent.parent
DATA = PLAYGROUND / "data"
OUT = PLAYGROUND / "outputs"
OUT.mkdir(exist_ok=True)

CCRES = DATA / "hg38_cCRE_chr22.bed"
H3K27AC = DATA / "K562_H3K27ac_chr22.bw"


def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


def load_region_set() -> RegionSet:
    section("1. Load chr22 cCREs into a gtars RegionSet")
    rs = RegionSet(str(CCRES))
    print(f"path:              {CCRES.relative_to(PLAYGROUND)}")
    print(f"identifier:        {rs.identifier}")
    print(f"n regions:         {len(rs):,}")
    print(f"mean region width: {rs.mean_region_width():.1f} bp")
    print(f"nucleotide length: {rs.get_nucleotide_length():,} bp "
          f"({rs.get_nucleotide_length() / 50_818_468 * 100:.2f}% of chr22)")
    print()
    print("The 'identifier' is a content-addressable digest of the region set.")
    print("Two BED files with the same regions (regardless of order or file")
    print("format details) produce the same identifier — Sheffield Lab's core")
    print("commitment to 'region sets as first-class objects' made concrete.")
    return rs


def neighbor_distance_distribution(rs: RegionSet):
    section("2. Inter-region spacing (gtars.neighbor_distances)")
    distances = rs.neighbor_distances()
    distances_sorted = sorted(distances)
    n = len(distances_sorted)
    print(f"neighbor distances computed for {n:,} regions")
    print()
    print("percentiles of inter-region distance (bp):")
    for pct in (5, 25, 50, 75, 95):
        idx = max(0, min(n - 1, int(n * pct / 100)))
        print(f"  p{pct:<3d}  {distances_sorted[idx]:>10,.0f}")
    print()
    under_1kb = sum(1 for d in distances if d < 1_000)
    under_10kb = sum(1 for d in distances if d < 10_000)
    print(f"cCREs within 1 kb of a neighbor:  {under_1kb:,}  ({under_1kb/n*100:.1f}%)")
    print(f"cCREs within 10 kb of a neighbor: {under_10kb:,}  ({under_10kb/n*100:.1f}%)")
    print()
    print("This is a spatial statistic you'd compute to characterize how")
    print("clustered the region set is. cCREs tend to cluster tightly near")
    print("active genes — the p25 distance here should be in the 100s of bp,")
    print("not the 1000s, because cCREs in the same enhancer-cluster sit")
    print("adjacent in the BED file.")


def parse_classes_and_filter(bed_path: Path, target: str):
    """Filter the BED by cCRE class and return parallel arrays for gtars."""
    section(f"3. Filter cCREs to class = {target} (SCREEN classification)")
    counts: dict[str, int] = {}
    chrs, starts, ends = [], [], []
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 10:
                continue
            label = parts[9]
            counts[label] = counts.get(label, 0) + 1
            if label == target:
                chrs.append(parts[0])
                starts.append(int(parts[1]))
                ends.append(int(parts[2]))

    total = sum(counts.values())
    print(f"total cCREs on chr22: {total:,}")
    print()
    print(f"{'class':<15s}  {'count':>8s}  {'pct':>6s}")
    for label in sorted(counts, key=lambda k: -counts[k]):
        pct = counts[label] / total * 100
        marker = " <-- filtered" if label == target else ""
        print(f"{label:<15s}  {counts[label]:>8,d}  {pct:>5.1f}%{marker}")
    return chrs, starts, ends


def build_pels_regionset(chrs, starts, ends) -> RegionSet:
    section("4. Build a new RegionSet from vectors (gtars.RegionSet.from_vectors)")
    pels = RegionSet.from_vectors(chrs, starts, ends)
    print(f"n regions:         {len(pels):,}")
    print(f"mean region width: {pels.mean_region_width():.1f} bp")
    print(f"nucleotide length: {pels.get_nucleotide_length():,} bp")
    print()
    print("from_vectors is the zero-copy way to build a RegionSet in memory")
    print("when you've already filtered or generated the intervals yourself")
    print("— no BED file round-trip required.")
    return pels


def set_comparison(full: RegionSet, subset: RegionSet):
    section("5. Set-theoretic comparison (gtars.jaccard / overlap_coefficient)")
    j = subset.jaccard(full)
    oc = subset.overlap_coefficient(full)
    print(f"jaccard(pELS, all cCREs):            {j:.4f}")
    print(f"overlap_coefficient(pELS, all cCREs): {oc:.4f}")
    print()
    print("Jaccard is |A ∩ B| / |A ∪ B|. Since pELS is a pure subset of the")
    print("full cCRE set, jaccard = |pELS| / |all| ≈ 0.16 (pELS is ~16% of")
    print("cCREs on chr22). Overlap coefficient is |A ∩ B| / min(|A|, |B|),")
    print("which equals 1.0 when A is a subset of B — which it is here.")


def score_with_pybigwig(pels: RegionSet, bigwig_path: Path):
    section("6. Score pELS regions with K562 H3K27ac (via pyBigWig)")
    # gtars RegionSet doesn't iterate to (chrom, start, end) directly yet,
    # so re-read the BED subset by iterating parallel vectors we kept.
    out_bed = OUT / "pELS_chr22.bed"
    # Export RegionSet to BED for downstream tools.
    pels.to_bed(str(out_bed))
    print(f"[write] {out_bed.relative_to(PLAYGROUND)}")

    bw = pyBigWig.open(str(bigwig_path))
    scores = []
    with open(out_bed) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            sig = bw.stats(chrom, start, end, type="mean")[0] or 0
            scores.append((sig, chrom, start, end))
    bw.close()

    scores.sort(reverse=True)
    print()
    print(f"scored {len(scores):,} pELS regions against K562 H3K27ac")
    print(f"median signal: {median([s[0] for s in scores]):.3f}")
    print()
    print("top 10 pELS by mean H3K27ac signal:")
    print(f"  {'rank':<5s} {'coords':<32s} {'mean signal':>12s}")
    for i, (sig, chrom, start, end) in enumerate(scores[:10], 1):
        coords = f"{chrom}:{start:,}-{end:,}"
        print(f"  {i:<5d} {coords:<32s} {sig:>12.2f}")
    print()
    print("Every rank-1 enhancer candidate here is now ready to drop into IGV")
    print("for visual inspection. Copy any coordinate above and paste it into")
    print("IGV's locus bar to see the multi-track signal shape at that site.")


def main():
    rs = load_region_set()
    neighbor_distance_distribution(rs)

    chrs, starts, ends = parse_classes_and_filter(CCRES, "pELS")
    pels = build_pels_regionset(chrs, starts, ends)
    set_comparison(rs, pels)
    score_with_pybigwig(pels, H3K27AC)


if __name__ == "__main__":
    main()
