"""
pyBigWig basics: open a bigwig, read values, compute region summaries.

This script walks through the core operations you care about when thinking
about "what does a bigwig contain" and "how do I pull signal at a position
or across a set of regions."

Runs against the chr22-only K562 bigwigs produced by the fetch step.
"""

from pathlib import Path
import pyBigWig

PLAYGROUND = Path(__file__).resolve().parent.parent
DATA = PLAYGROUND / "data"

TRACKS = {
    "ATAC":    DATA / "K562_ATAC_chr22.bw",
    "H3K27ac": DATA / "K562_H3K27ac_chr22.bw",
    "H3K4me1": DATA / "K562_H3K4me1_chr22.bw",
}


def section(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


def main():
    # --- 1. Open a bigwig and inspect its header ---------------------------
    section("1. Open a bigwig and read its header")
    bw = pyBigWig.open(str(TRACKS["H3K27ac"]))
    print(f"chroms: {bw.chroms()}")
    print(f"is_bigwig: {bw.isBigWig()}")
    print(f"header: {bw.header()}")
    print()
    print("A bigwig is an indexed binary file with one continuous signal value")
    print("per genomic position (stored in compressed blocks). The header tells")
    print("you which chromosomes it covers and what their lengths are.")
    bw.close()

    # --- 2. Read signal at a single locus ---------------------------------
    section("2. Read signal values in a narrow window")
    locus_start, locus_end = 38_000_000, 38_000_200  # arbitrary chr22 locus
    print(f"locus: chr22:{locus_start:,}-{locus_end:,}  ({locus_end - locus_start} bp)")
    for name, path in TRACKS.items():
        bw = pyBigWig.open(str(path))
        values = bw.values("chr22", locus_start, locus_end)
        mean_val = sum(v for v in values if v is not None) / len(values)
        print(f"  {name:8s}  mean = {mean_val:.3f}  first 5 bp = {[round(v, 2) if v else 0 for v in values[:5]]}")
        bw.close()
    print()
    print("`bw.values(chrom, start, end)` returns a list of per-base signal")
    print("values. None means 'no data' (outside any indexed block, i.e. zero).")

    # --- 3. Region-level summaries via .stats() ---------------------------
    section("3. Region-level summaries via .stats()")
    region_start, region_end = 30_000_000, 40_000_000  # 10 Mb chunk of chr22
    print(f"region: chr22:{region_start:,}-{region_end:,}  ({(region_end - region_start)/1e6:.0f} Mb)")
    print()
    print(f"{'track':<10s}  {'mean':>8s}  {'max':>8s}  {'stdev':>8s}  {'cov_pct':>8s}")
    for name, path in TRACKS.items():
        bw = pyBigWig.open(str(path))
        mean_v = bw.stats("chr22", region_start, region_end, type="mean")[0] or 0
        max_v = bw.stats("chr22", region_start, region_end, type="max")[0] or 0
        std_v = bw.stats("chr22", region_start, region_end, type="std")[0] or 0
        cov_v = bw.stats("chr22", region_start, region_end, type="coverage")[0] or 0
        print(f"{name:<10s}  {mean_v:>8.3f}  {max_v:>8.3f}  {std_v:>8.3f}  {cov_v*100:>7.1f}%")
        bw.close()
    print()
    print(".stats() is the fast path for 'just give me the summary' — much")
    print("cheaper than pulling all base-level values and computing them yourself.")

    # --- 4. Cross-track profile at a single position ---------------------
    section("4. The row-view: all tracks at a single position (the 'stack')")
    position = 32_000_000  # arbitrary chr22 position
    window = 500
    print(f"position: chr22:{position:,}  (+/- {window} bp window)")
    print()
    print(f"{'track':<10s}  {'signal (mean over window)':>30s}")
    for name, path in TRACKS.items():
        bw = pyBigWig.open(str(path))
        signal = bw.stats("chr22", position - window, position + window, type="mean")[0] or 0
        print(f"{name:<10s}  {signal:>30.4f}")
        bw.close()
    print()
    print("This is the 'what does the regulatory stack look like at this position'")
    print("view — one number per track, same position. Pick three positions with")
    print("different multi-track patterns and you see three different regulatory")
    print("states. Pick the same position across many tissues and you see its")
    print("tissue-specificity.")

    # --- 5. Read BED regions and score each one --------------------------
    section("5. Score a BED file of regions against one bigwig")
    bed_path = DATA / "hg38_cCRE_chr22.bed"
    if not bed_path.exists():
        print(f"(skipping; {bed_path.name} not found)")
        return

    bw = pyBigWig.open(str(TRACKS["H3K27ac"]))
    n_scored = 0
    top = []  # (score, name, chrom, start, end, label)
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) > 3 else "."
            label = parts[9] if len(parts) > 9 else "."
            score = bw.stats(chrom, start, end, type="mean")[0] or 0
            top.append((score, name, chrom, start, end, label))
            n_scored += 1
    bw.close()

    top.sort(reverse=True)
    print(f"scored {n_scored:,} cCREs against K562 H3K27ac")
    print()
    print("Top 10 cCREs by mean H3K27ac signal:")
    print(f"  {'rank':<5s} {'cCRE':<15s} {'coords':<30s} {'class':<12s} {'signal':>8s}")
    for i, (score, name, chrom, start, end, label) in enumerate(top[:10], 1):
        coords = f"{chrom}:{start:,}-{end:,}"
        print(f"  {i:<5d} {name:<15s} {coords:<30s} {label:<12s} {score:>8.2f}")
    print()
    print("This is the column-view of a bigwig summarized over a region set:")
    print("one score per region for one track. It's the simplest profile")
    print("statistic and the primitive every more elaborate analysis builds on.")


if __name__ == "__main__":
    main()
