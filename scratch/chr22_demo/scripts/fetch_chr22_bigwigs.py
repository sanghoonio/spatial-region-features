"""
Fetch chr22-only bigwigs for the playground.

Strategy: download the full-genome bigwigs from ENCODE into data/_full/,
then use pyBigWig to read each one and write a chr22-only copy into data/.
The full files are kept around in case you want to branch out beyond chr22.

Source experiments (all K562, GRCh38):
  - H3K27ac:  ENCSR000AKP  ->  ENCFF467OGB.bigWig  (signal p-value)
  - H3K4me1:  ENCSR000EWC  ->  ENCFF948OCT.bigWig  (signal p-value)
  - ATAC-seq: ENCSR483RKN  ->  ENCFF867QEW.bigWig  (signal p-value)

Writes:
  data/_full/{name}_full.bw    (original full-genome files, ~300-500 MB each)
  data/{name}_chr22.bw         (chr22-only slice, ~5-10 MB each)
"""

from pathlib import Path
import subprocess
import sys
import pyBigWig

PLAYGROUND = Path(__file__).resolve().parent.parent
DATA = PLAYGROUND / "data"
FULL = DATA / "_full"

SOURCES = [
    ("K562_H3K27ac", "https://www.encodeproject.org/files/ENCFF467OGB/@@download/ENCFF467OGB.bigWig"),
    ("K562_H3K4me1", "https://www.encodeproject.org/files/ENCFF948OCT/@@download/ENCFF948OCT.bigWig"),
    ("K562_ATAC",    "https://www.encodeproject.org/files/ENCFF867QEW/@@download/ENCFF867QEW.bigWig"),
]

CHROM = "chr22"


def download_full(name: str, url: str) -> Path:
    FULL.mkdir(parents=True, exist_ok=True)
    out = FULL / f"{name}_full.bw"
    if out.exists() and out.stat().st_size > 1_000_000:
        print(f"[skip] {out.name} already present ({out.stat().st_size/1e6:.0f} MB)")
        return out
    print(f"[curl] {name} <- {url}")
    subprocess.run(["curl", "-sSL", "-o", str(out), url], check=True)
    print(f"[done] {out.name}  {out.stat().st_size/1e6:.0f} MB")
    return out


def slice_to_chr22(name: str, full_path: Path) -> Path:
    out_path = DATA / f"{name}_chr22.bw"
    if out_path.exists() and out_path.stat().st_size > 1000:
        print(f"[skip] {out_path.name} already exists")
        return out_path

    print(f"[open] {full_path.name}")
    src = pyBigWig.open(str(full_path))
    if src is None:
        raise RuntimeError(f"failed to open {full_path}")

    chroms = src.chroms()
    if CHROM not in chroms:
        src.close()
        raise RuntimeError(f"{CHROM} not in {full_path.name}")

    chr22_len = chroms[CHROM]
    print(f"[read] {CHROM} length = {chr22_len:,}")

    intervals = src.intervals(CHROM)
    src.close()

    if not intervals:
        raise RuntimeError(f"no intervals on {CHROM} in {full_path.name}")

    print(f"[read] {len(intervals):,} intervals")

    dst = pyBigWig.open(str(out_path), "w")
    dst.addHeader([(CHROM, chr22_len)])
    dst.addEntries(
        [CHROM] * len(intervals),
        [iv[0] for iv in intervals],
        ends=[iv[1] for iv in intervals],
        values=[float(iv[2]) for iv in intervals],
    )
    dst.close()

    print(f"[write] {out_path.name}  {out_path.stat().st_size/1e6:.2f} MB")
    return out_path


def main():
    for name, url in SOURCES:
        full = download_full(name, url)
        slice_to_chr22(name, full)
        print()

    print("chr22 bigwigs in data/:")
    for p in sorted(DATA.glob("*_chr22.bw")):
        print(f"  {p.name}  {p.stat().st_size/1e6:.2f} MB")


if __name__ == "__main__":
    main()
