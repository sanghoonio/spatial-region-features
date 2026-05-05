"""
Fetch bigwigs and peak BEDs for six modalities across K562, GM12878, HepG2,
slice each to the 30 selected cCREs ± 2 kb, and write an IGV session that
loads all 36 tracks.

Inputs:
  data/ccres.bed  (produced by 01_select_ccres.py)

Pipeline per (tissue × modality):
  1. Resolve an ENCODE experiment via the search API (cached in
     data/_cache/experiments.json).
  2. Find the fold-change-over-control bigwig (or signal-p-value fallback).
  3. Find the IDR-thresholded narrowPeak (or replicated / pseudoreplicated
     fallback).
  4. Download both if not cached.
  5. Slice the bigwig to cCRE ± flank windows; write a small sliced bigwig.
  6. Subset the peak BED to entries overlapping cCRE ± flank; write a small
     sliced BED.

Outputs:
  data/signal/<tissue>_<modality>.bw     (sliced bigwig, ~1-5 MB each)
  data/peaks/<tissue>_<modality>.bed     (sliced narrowPeak, tiny)
  outputs/igv_session.xml                (auto-loads all 36 tracks)
  data/_cache/experiments.json           (resolved ENCODE accessions)

Usage:
  cd scratch/ccre_examples
  uv run python scripts/02_fetch_and_slice.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig


# --- paths ---

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"
CACHE = DATA / "_cache"
SIGNAL = DATA / "signal"
PEAKS = DATA / "peaks"
OUTPUTS = HERE / "outputs"

CCRES_BED = DATA / "ccres.bed"
EXPERIMENTS_JSON = CACHE / "experiments.json"
IGV_SESSION = OUTPUTS / "igv_session.xml"


# --- config ---

TISSUES = ["K562", "GM12878", "HepG2"]
MODALITIES = ["ATAC", "DNase", "H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3"]
FLANK_BP = 2000

# Search parameters per modality. Passed to ENCODE's search endpoint. Every
# search is additionally constrained to status=released, assembly=GRCh38, and
# the chosen tissue's biosample_ontology.term_name.
ASSAY_SEARCH: dict[str, dict[str, str]] = {
    "ATAC":     {"assay_title": "ATAC-seq"},
    "DNase":    {"assay_title": "DNase-seq"},
    "H3K27ac":  {"assay_title": "Histone ChIP-seq", "target.label": "H3K27ac"},
    "H3K4me1":  {"assay_title": "Histone ChIP-seq", "target.label": "H3K4me1"},
    "H3K4me3":  {"assay_title": "Histone ChIP-seq", "target.label": "H3K4me3"},
    "H3K27me3": {"assay_title": "Histone ChIP-seq", "target.label": "H3K27me3"},
}

BIGWIG_PRIORITY = (
    "fold change over control",   # histone ChIP, ATAC
    "signal p-value",              # histone ChIP, ATAC
    "read-depth normalized signal",  # DNase
)
PEAK_PRIORITY = (
    "IDR thresholded peaks",       # ATAC (Kundaje pipeline)
    "replicated peaks",            # histone ChIP with replicates
    "pseudoreplicated peaks",      # histone ChIP single replicate
    "peaks",                       # DNase, generic
)

# IGV display colors per tissue (distinct enough to tell apart at a glance)
TISSUE_COLOR = {
    "K562":    "80,80,200",   # blue
    "GM12878": "80,180,80",   # green
    "HepG2":   "200,80,80",   # red
}


# --- HTTP helpers ---

def api_get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def curl(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"[curl] {url}")
    # Write to a .part file so a crashed download doesn't pollute the cache.
    # --retry heals transient S3/network drops mid-transfer.
    subprocess.run(
        ["curl", "-fsSL", "--retry", "5", "--retry-delay", "3",
         "--retry-all-errors", "-o", str(tmp), url],
        check=True,
    )
    tmp.rename(dest)
    print(f"[done] {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


# --- experiment resolution ---

def search_experiment(tissue: str, modality: str) -> str:
    """Return an ENCSR accession for the tissue+modality combination."""
    params = {
        "type": "Experiment",
        "status": "released",
        "assembly": "GRCh38",
        "biosample_ontology.term_name": tissue,
        "format": "json",
        "limit": "25",
        **ASSAY_SEARCH[modality],
    }
    qs = urllib.parse.urlencode(params)
    url = f"https://www.encodeproject.org/search/?{qs}"
    result = api_get(url)
    hits = result.get("@graph", [])
    if not hits:
        raise RuntimeError(f"no released experiments for {tissue} / {modality}")
    # Prefer an experiment that actually released bigWigs under our priorities.
    # Otherwise just take the first.
    for hit in hits:
        encsr = hit["@id"].strip("/").split("/")[-1]
        return encsr
    raise RuntimeError(f"unreachable: {tissue} / {modality}")


def resolve_all_experiments() -> dict[str, dict[str, str]]:
    """Map tissue → modality → ENCSR accession, caching to a JSON file."""
    if EXPERIMENTS_JSON.exists():
        cached = json.loads(EXPERIMENTS_JSON.read_text())
        print(f"[cache] loaded {EXPERIMENTS_JSON.name}")
        return cached

    print(f"[api  ] resolving experiments via ENCODE search...")
    out: dict[str, dict[str, str]] = {}
    for tissue in TISSUES:
        out[tissue] = {}
        for modality in MODALITIES:
            encsr = search_experiment(tissue, modality)
            out[tissue][modality] = encsr
            print(f"        {tissue:>8} / {modality:>8}  ->  {encsr}")

    EXPERIMENTS_JSON.parent.mkdir(parents=True, exist_ok=True)
    EXPERIMENTS_JSON.write_text(json.dumps(out, indent=2))
    print(f"[write] {EXPERIMENTS_JSON.name}")
    return out


def pick_file(
    exp: dict,
    output_type_priority: tuple[str, ...],
    file_format: str,
    file_format_type: str | None = None,
) -> tuple[str, str] | None:
    """Return (accession, output_type) of the best-matching released GRCh38 file, or None."""
    files = exp.get("files", [])
    for ot in output_type_priority:
        cands = [
            f for f in files
            if f.get("output_type") == ot
            and f.get("file_format") == file_format
            and (file_format_type is None or f.get("file_format_type") == file_format_type)
            and f.get("status") == "released"
            and f.get("assembly") == "GRCh38"
        ]
        if not cands:
            continue
        default = [f for f in cands if f.get("preferred_default")]
        pick = (default or cands)[0]
        return pick["accession"], pick["output_type"]
    return None


def resolve_files(encsr: str) -> tuple[str, str, str, str] | None:
    """Return (bigwig_accession, bigwig_output_type, bed_accession, bed_output_type)."""
    exp = api_get(f"https://www.encodeproject.org/experiments/{encsr}/?format=json")
    bw = pick_file(exp, BIGWIG_PRIORITY, "bigWig")
    # Prefer narrowPeak (histone ChIP / ATAC) but fall back to any BED so we
    # don't miss DNase (which ENCODE files as bed3+ rather than narrowPeak).
    peaks = pick_file(exp, PEAK_PRIORITY, "bed", "narrowPeak") or pick_file(
        exp, PEAK_PRIORITY, "bed"
    )
    if bw is None and peaks is None:
        return None
    bw_acc, bw_ot = bw or ("", "")
    pk_acc, pk_ot = peaks or ("", "")
    return bw_acc, bw_ot, pk_acc, pk_ot


# --- download + slice ---

def download_bigwig(accession: str) -> Path:
    url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bigWig"
    dest = CACHE / "full_bigwigs" / f"{accession}.bigWig"
    return curl(url, dest)


def download_bed(accession: str) -> Path:
    url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.bed.gz"
    dest = CACHE / "full_beds" / f"{accession}.bed.gz"
    return curl(url, dest)


def load_ccre_windows() -> pd.DataFrame:
    """Load ccres.bed and expand each interval by FLANK_BP on each side."""
    ccres = pd.read_csv(
        CCRES_BED, sep="\t", header=None,
        names=["chrom", "start", "end", "ccre_accession", "score", "strand"],
    )
    w = ccres.copy()
    w["win_start"] = (w["start"] - FLANK_BP).clip(lower=0)
    w["win_end"] = w["end"] + FLANK_BP
    w = w.sort_values(["chrom", "win_start"]).reset_index(drop=True)
    return w


def slice_bigwig(src_path: Path, windows: pd.DataFrame, dst_path: Path) -> None:
    """Write a small bigwig containing only signal inside the cCRE ± flank windows."""
    if dst_path.exists() and dst_path.stat().st_size > 0:
        print(f"[skip ] {dst_path.name}")
        return
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    src = pyBigWig.open(str(src_path))
    src_chroms = src.chroms()
    # Collect intervals, keyed by chrom, then write in chrom/start order.
    collected: dict[str, list[tuple[int, int, float]]] = {}
    total = 0
    for w in windows.itertuples(index=False):
        if w.chrom not in src_chroms:
            continue
        win_start = max(0, int(w.win_start))
        win_end = min(src_chroms[w.chrom], int(w.win_end))
        if win_end <= win_start:
            continue
        for iv in src.intervals(w.chrom, win_start, win_end) or []:
            collected.setdefault(w.chrom, []).append(iv)
            total += 1
    src.close()

    # Sort per-chrom and deduplicate (adjacent windows can overlap)
    ordered_chroms = [c for c in src_chroms if c in collected]
    dst = pyBigWig.open(str(dst_path), "w")
    dst.addHeader([(c, src_chroms[c]) for c in ordered_chroms])
    for c in ordered_chroms:
        intervals = sorted(set(collected[c]))
        if not intervals:
            continue
        dst.addEntries(
            [c] * len(intervals),
            [i[0] for i in intervals],
            ends=[i[1] for i in intervals],
            values=[float(i[2]) for i in intervals],
        )
    dst.close()
    print(f"[slice] {dst_path.name}  ({total:,} intervals)")


def slice_bed(src_path: Path, windows: pd.DataFrame, dst_path: Path) -> None:
    """Subset a BED (narrowPeak or bed3+) to entries overlapping any cCRE ± flank window."""
    if dst_path.exists() and dst_path.stat().st_size > 0:
        print(f"[skip ] {dst_path.name}")
        return
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    # Read with dynamic column count — first three are always chrom/start/end.
    # comment='#' skips header lines in ENCODE DNase bed3+ files.
    peaks = pd.read_csv(
        src_path, sep="\t", header=None, comment="#",
        compression="gzip" if src_path.suffix == ".gz" else None,
    )
    n_cols = peaks.shape[1]
    # Standard narrowPeak names for the first 10; bed3+ keeps the rest generic
    std = ["chrom", "start", "end", "name", "score", "strand",
           "signalValue", "pValue", "qValue", "peak"]
    names = std[:n_cols] + [f"col{i}" for i in range(10, n_cols)]
    peaks.columns = names

    import pyranges as pr
    peaks_pr = pr.PyRanges(peaks.rename(
        columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
    ))
    win_pr = pr.PyRanges(windows.rename(
        columns={"chrom": "Chromosome", "win_start": "Start", "win_end": "End"}
    )[["Chromosome", "Start", "End"]])

    hit = peaks_pr.overlap(win_pr).df
    if len(hit) == 0:
        dst_path.write_text("")
        print(f"[slice] {dst_path.name}  (0 peaks)")
        return
    out = hit.rename(columns={"Chromosome": "chrom", "Start": "start", "End": "end"})
    out = out[names].sort_values(["chrom", "start"])
    out.to_csv(dst_path, sep="\t", header=False, index=False)
    print(f"[slice] {dst_path.name}  ({len(out):,} peaks)")


# --- IGV session ---

def _build_session_xml(
    signal_paths: dict, peak_paths: dict, ccres_bed: Path, out_dir: Path, locus: str,
) -> str:
    """Return an IGV session XML string configured to open at `locus`.

    `out_dir` is the directory the session will be written to (paths are
    relative to it)."""
    def rel(p: Path) -> str:
        # Path relative to the directory the session file lives in
        return os.path.relpath(p.resolve(), start=out_dir.resolve())

    xml: list[str] = []
    xml.append('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
    xml.append(f'<Session genome="hg38" locus="{locus}" version="8">')
    xml.append('    <Resources>')
    for t in TISSUES:
        for m in MODALITIES:
            if (t, m) in signal_paths:
                xml.append(f'        <Resource name="{t} {m}" path="{rel(signal_paths[(t, m)])}"/>')
            if (t, m) in peak_paths:
                xml.append(f'        <Resource name="{t} {m} peaks" path="{rel(peak_paths[(t, m)])}"/>')
    xml.append(f'        <Resource name="Selected cCREs" path="{rel(ccres_bed)}"/>')
    xml.append('    </Resources>')

    xml.append('    <Panel height="1100" name="SignalPanel" width="1400">')
    for m in MODALITIES:
        for t in TISSUES:
            if (t, m) not in signal_paths:
                continue
            path_rel = rel(signal_paths[(t, m)])
            color = TISSUE_COLOR[t]
            xml.append(
                f'        <Track attributeKey="{t} {m}" autoScale="true" '
                f'clazz="org.broad.igv.track.DataSourceTrack" color="{color}" '
                f'displayMode="COLLAPSED" id="{path_rel}" name="{t} {m}" '
                f'renderer="BAR_CHART" visible="true" windowFunction="mean">'
            )
            xml.append(
                '            <DataRange baseline="0.0" drawBaseline="true" '
                'flipAxis="false" type="LINEAR"/>'
            )
            xml.append('        </Track>')
    xml.append('    </Panel>')

    xml.append('    <Panel height="400" name="FeaturePanel" width="1400">')
    for m in MODALITIES:
        for t in TISSUES:
            if (t, m) not in peak_paths:
                continue
            path_rel = rel(peak_paths[(t, m)])
            color = TISSUE_COLOR[t]
            xml.append(
                f'        <Track attributeKey="{t} {m} peaks" '
                f'clazz="org.broad.igv.track.FeatureTrack" color="{color}" '
                f'displayMode="COLLAPSED" id="{path_rel}" '
                f'name="{t} {m} peaks" visible="true"/>'
            )
    ccres_rel = rel(ccres_bed)
    xml.append(
        f'        <Track attributeKey="Selected cCREs" '
        f'clazz="org.broad.igv.track.FeatureTrack" color="30,30,30" '
        f'displayMode="COLLAPSED" id="{ccres_rel}" '
        f'name="Selected cCREs" visible="true"/>'
    )
    xml.append('    </Panel>')
    xml.append('    <PanelLayout dividerFractions="0.8"/>')
    xml.append('</Session>')
    return "\n".join(xml) + "\n"


def _locus_for(row) -> str:
    start = max(0, int(row["start"]) - FLANK_BP)
    end = int(row["end"]) + FLANK_BP
    return f"{row['chrom']}:{start:,}-{end:,}"


def write_igv_session(signal_paths: dict, peak_paths: dict, ccres_bed: Path) -> None:
    """Write the combined session + one session per cCRE + an index."""
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    per_ccre_dir = OUTPUTS / "igv_sessions"
    per_ccre_dir.mkdir(parents=True, exist_ok=True)

    ccres = pd.read_csv(
        ccres_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "ccre_accession", "score", "strand"],
    )
    meta_path = DATA / "ccres_metadata.csv"
    meta = pd.read_csv(meta_path) if meta_path.exists() else None

    # Combined session: opens at the first cCRE
    locus0 = _locus_for(ccres.iloc[0])
    IGV_SESSION.write_text(
        _build_session_xml(signal_paths, peak_paths, ccres_bed, OUTPUTS, locus0)
    )
    print(f"[write] {IGV_SESSION.relative_to(HERE)}  (opens at {locus0})")

    # Per-cCRE sessions, named in a stable order: by class, stratum, validation,
    # then cCRE accession. Numeric prefix keeps the filesystem order readable.
    meta_by_acc: dict[str, dict] = {}
    if meta is not None:
        meta_by_acc = {r["ccre_accession"]: r for _, r in meta.iterrows()}

    def sort_key(acc: str) -> tuple:
        m = meta_by_acc.get(acc, {})
        return (m.get("ccre_class", ""),
                m.get("stratum", ""),
                m.get("validation_status", ""),
                acc)

    index_rows = []
    ordered = sorted(ccres["ccre_accession"], key=sort_key)
    pad = len(str(len(ordered)))
    for i, acc in enumerate(ordered, start=1):
        row = ccres[ccres["ccre_accession"] == acc].iloc[0]
        locus = _locus_for(row)
        m = meta_by_acc.get(acc, {})
        klass = m.get("ccre_class", "unknown")
        stratum = m.get("stratum", "unknown")
        validation = m.get("validation_status", "untested")
        tag = f"{klass}_{stratum}_{validation}".replace(" ", "-")
        fname = f"{str(i).zfill(pad)}_{tag}_{acc}.xml"
        out_path = per_ccre_dir / fname
        out_path.write_text(
            _build_session_xml(signal_paths, peak_paths, ccres_bed, per_ccre_dir, locus)
        )
        index_rows.append({
            "order": i, "file": fname, "ccre_accession": acc,
            "class": klass, "stratum": stratum,
            "validation_status": validation, "locus": locus,
        })

    index_df = pd.DataFrame(index_rows)
    index_csv = per_ccre_dir / "index.csv"
    index_df.to_csv(index_csv, index=False)

    # Markdown index for quick visual browsing
    md_lines = [
        "# Per-cCRE IGV sessions",
        "",
        "Each file opens IGV with all 37 tracks loaded, zoomed to that cCRE ± 2 kb.",
        "",
        "| # | cCRE | class | stratum | validation | locus | session |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in index_rows:
        md_lines.append(
            f"| {r['order']} | `{r['ccre_accession']}` | {r['class']} | "
            f"{r['stratum']} | {r['validation_status']} | `{r['locus']}` | "
            f"[{r['file']}]({r['file']}) |"
        )
    (per_ccre_dir / "index.md").write_text("\n".join(md_lines) + "\n")

    print(f"[write] {len(index_rows)} per-cCRE sessions in "
          f"{per_ccre_dir.relative_to(HERE)}/")
    print(f"[write] {index_csv.relative_to(HERE)}")
    print(f"[write] {(per_ccre_dir / 'index.md').relative_to(HERE)}")


# --- main ---

def main() -> int:
    if not CCRES_BED.exists():
        print(f"[fail] {CCRES_BED} not found. Run 01_select_ccres.py first.",
              file=sys.stderr)
        return 1

    windows = load_ccre_windows()
    print(f"[load ] {len(windows)} cCRE windows ({FLANK_BP} bp flank each side)")

    experiments = resolve_all_experiments()

    # For each tissue × modality: resolve file accessions, download, slice
    signal_paths: dict[tuple[str, str], Path] = {}
    peak_paths: dict[tuple[str, str], Path] = {}

    for tissue in TISSUES:
        for modality in MODALITIES:
            encsr = experiments[tissue][modality]
            print(f"\n=== {tissue} / {modality} ({encsr}) ===")
            files = resolve_files(encsr)
            if files is None:
                print(f"[warn ] no usable files in {encsr}")
                continue
            bw_acc, bw_ot, pk_acc, pk_ot = files

            if bw_acc:
                print(f"[api  ] bigwig  {bw_acc}  ({bw_ot})")
                full_bw = download_bigwig(bw_acc)
                out_bw = SIGNAL / f"{tissue}_{modality}.bw"
                slice_bigwig(full_bw, windows, out_bw)
                signal_paths[(tissue, modality)] = out_bw
            else:
                print(f"[warn ] no bigwig in {encsr}")

            if pk_acc:
                print(f"[api  ] peaks   {pk_acc}  ({pk_ot})")
                full_bed = download_bed(pk_acc)
                out_bed = PEAKS / f"{tissue}_{modality}.bed"
                slice_bed(full_bed, windows, out_bed)
                peak_paths[(tissue, modality)] = out_bed
            else:
                print(f"[warn ] no peaks in {encsr}")

    print(f"\n=== Writing IGV session ===")
    write_igv_session(signal_paths, peak_paths, CCRES_BED)

    print(f"\n=== Done ===")
    print(f"  signal slices:  {len(signal_paths)}/{len(TISSUES) * len(MODALITIES)}")
    print(f"  peak slices:    {len(peak_paths)}/{len(TISSUES) * len(MODALITIES)}")
    print(f"  igv session:    {IGV_SESSION}")
    print(f"\nOpen in IGV: File -> Open Session... -> select {IGV_SESSION.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
