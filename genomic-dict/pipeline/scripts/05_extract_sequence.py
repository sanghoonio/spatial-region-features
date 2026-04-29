"""Stage 05 — sequence-based intrinsic features for chr16 regions.

Per chr16 region, derive: TF motif hits (JASPAR + FIMO), ChromHMM state per
tissue (Roadmap 15-state), phastCons mean.

Reads:
  data/annotations/pretrained_universe.parquet       (chr16 subset; coords)
  Downloaded: JASPAR 2024 MEME, hg38 chr16 FASTA, Roadmap ChromHMM BEDs,
              UCSC phastCons100way bigwig
Writes:
  data/annotations/intrinsic_sequence.parquet
  results/05_extract_sequence/summary.json
"""
from __future__ import annotations

import gzip
import shutil
import sys
import urllib.request
from pathlib import Path

import polars as pl
import requests

sys.path.insert(0, str(Path(__file__).parent))
from _common import PROJECT_ROOT, file_record, stage_start, write_summary  # noqa: E402


def download(url: str, dest: Path) -> Path:
    if dest.exists():
        print(f"  cached: {dest.name}", file=sys.stderr)
        return dest
    print(f"  downloading {url}", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=(10, 600)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as out:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    out.write(chunk)
    tmp.rename(dest)
    return dest


def fetch_chr16_fasta(work_dir: Path) -> Path:
    """Download UCSC hg38 chr16 FASTA (uncompressed, ~90 MB)."""
    fa = work_dir / "chr16.fa"
    if fa.exists():
        return fa
    gz = work_dir / "chr16.fa.gz"
    if not gz.exists():
        download("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr16.fa.gz", gz)
    print("  decompressing chr16.fa.gz", file=sys.stderr)
    with gzip.open(gz, "rb") as fin, open(fa, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    return fa


def extract_chr16_seqs(fa: Path) -> str:
    """Return chr16 sequence (uppercase) as a single string."""
    print("  loading chr16 sequence into memory", file=sys.stderr)
    seq_parts = []
    with open(fa) as f:
        next(f)  # skip header
        for line in f:
            seq_parts.append(line.strip())
    return "".join(seq_parts).upper()


def fimo_scan(meme_path: Path, regions_df: pl.DataFrame, chr16_seq: str,
              p_threshold: float, top_n: int) -> dict[int, list[tuple[str, float]]]:
    """For each region, top-N motif hits below p_threshold. Returns {token_id: [(motif_name, pvalue), ...]}."""
    from pymemesuite.common import MotifFile, Sequence
    from pymemesuite.fimo import FIMO

    print(f"  loading motifs from {meme_path.name}", file=sys.stderr)
    with MotifFile(str(meme_path)) as mf:
        motifs = list(mf)
        bg = mf.background
    print(f"    {len(motifs)} motifs loaded", file=sys.stderr)

    fimo = FIMO(threshold=p_threshold, both_strands=True)
    out: dict[int, list[tuple[str, float]]] = {}

    n = len(regions_df)
    rows = regions_df.iter_rows(named=True)
    for i, r in enumerate(rows):
        seq_str = chr16_seq[r["start"]:r["end"]]
        if not seq_str or "N" in seq_str[:10]:  # skip Ns
            out[r["token_id"]] = []
            continue
        seq = Sequence(seq_str.encode(), name=str(r["token_id"]).encode())
        hits = []
        for motif in motifs:
            res = fimo.score_motif(motif, [seq], bg)
            for m in res.matched_elements:
                hits.append((motif.name.decode(), float(m.pvalue)))
        hits.sort(key=lambda h: h[1])
        out[r["token_id"]] = hits[:top_n]
        if (i + 1) % 1000 == 0:
            print(f"    motif scan: {i+1:,}/{n:,}", file=sys.stderr)
    return out


def chromhmm_states(regions_df: pl.DataFrame, chromhmm_paths: dict[str, Path]) -> dict[str, dict[int, str]]:
    """Per tissue, return token_id -> majority-vote ChromHMM state."""
    out: dict[str, dict[int, str]] = {}
    for tissue, path in chromhmm_paths.items():
        print(f"  intersecting ChromHMM for {tissue}", file=sys.stderr)
        # Roadmap files are bed.gz with cols: chrom, start, end, state_label
        chrom = pl.read_csv(path, separator="\t", has_header=False,
                            new_columns=["chrom", "start", "end", "state"])
        chrom = chrom.filter(pl.col("chrom") == "chr16")
        # Per region, find overlapping states. Use polars range-style approach via sort + scan.
        out[tissue] = {}
        # Naive: for each region, filter chromhmm rows that overlap. Slow for large universes
        # but chr16 chromhmm is ~30k rows and regions ~36k; ~1B comparisons too slow.
        # Use pyranges instead.
        import pyranges as pr
        regions_pr = pr.PyRanges(regions_df.select(["chrom", "start", "end", "token_id"])
                                 .to_pandas().rename(columns={"chrom": "Chromosome",
                                                              "start": "Start",
                                                              "end": "End"}))
        chrom_pr = pr.PyRanges(chrom.to_pandas().rename(columns={"chrom": "Chromosome",
                                                                  "start": "Start",
                                                                  "end": "End"}))
        joined = regions_pr.join(chrom_pr, suffix="_state").df
        if len(joined) == 0:
            continue
        joined["overlap"] = (joined[["End", "End_state"]].min(axis=1)
                             - joined[["Start", "Start_state"]].max(axis=1)).clip(lower=0)
        best = (joined.sort_values(["token_id", "overlap"], ascending=[True, False])
                .drop_duplicates(subset=["token_id"], keep="first"))
        out[tissue] = dict(zip(best["token_id"], best["state"]))
    return out


def phastcons_means(regions_df: pl.DataFrame, bw_path: Path) -> dict[int, float]:
    """Mean phastCons over each region."""
    import pyBigWig
    print(f"  reading phastCons", file=sys.stderr)
    bw = pyBigWig.open(str(bw_path))
    out: dict[int, float] = {}
    for r in regions_df.iter_rows(named=True):
        try:
            v = bw.stats("chr16", int(r["start"]), int(r["end"]), type="mean", nBins=1)[0]
            out[r["token_id"]] = float(v) if v is not None else float("nan")
        except Exception:
            out[r["token_id"]] = float("nan")
    bw.close()
    return out


def main() -> None:
    ctx = stage_start("05_extract_sequence", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    chrom = sc["focus_chromosome"]
    pretrained_path = PROJECT_ROOT / paths["annotations_dir"] / "pretrained_universe.parquet"
    annot_dir = PROJECT_ROOT / paths["annotations_dir"]
    out_path = annot_dir / "intrinsic_sequence.parquet"
    work_dir = annot_dir / "stage_05_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    inputs = [file_record(pretrained_path)]

    try:
        if not pretrained_path.exists():
            raise FileNotFoundError(f"{pretrained_path} missing.")

        regions = (pl.read_parquet(pretrained_path)
                   .filter(pl.col("chrom") == chrom)
                   .select(["token_id", "chrom", "start", "end"])
                   .sort(["start", "end"]))
        print(f"  {len(regions):,} chr16 regions to annotate", file=sys.stderr)

        # --- Conservation ---
        phast_url = sc.get("phastcons_url", "")
        phast_means: dict[int, float] = {}
        if phast_url:
            phast_path = work_dir / "phastCons100way.chr16.bw"
            # Whole-genome phastCons is ~10 GB; for chr16 use UCSC's per-chr file if available.
            # Otherwise download whole and accept the size.
            download(phast_url, phast_path)
            phast_means = phastcons_means(regions, phast_path)
            inputs.append(file_record(phast_path))

        # --- ChromHMM ---
        chromhmm_states_dict: dict[str, dict[int, str]] = {}
        chromhmm_urls = sc.get("chromhmm_urls", {})
        if chromhmm_urls:
            chromhmm_paths = {}
            for tissue, url in chromhmm_urls.items():
                p = work_dir / f"chromhmm_{tissue}.bed.gz"
                download(url, p)
                chromhmm_paths[tissue] = p
                inputs.append(file_record(p))
            chromhmm_states_dict = chromhmm_states(regions, chromhmm_paths)

        # --- Motifs ---
        motif_hits: dict[int, list[tuple[str, float]]] = {}
        jaspar_url = sc.get("jaspar_motifs_url", "")
        if jaspar_url:
            meme_path = work_dir / "jaspar_2024.meme"
            download(jaspar_url, meme_path)
            inputs.append(file_record(meme_path))
            chr16_fa = fetch_chr16_fasta(work_dir)
            inputs.append(file_record(chr16_fa))
            chr16_seq = extract_chr16_seqs(chr16_fa)
            p_thresh = float(sc.get("fimo_pvalue_threshold", 1e-4))
            top_n = int(sc.get("fimo_top_n_per_ccre", 10))
            motif_hits = fimo_scan(meme_path, regions, chr16_seq, p_thresh, top_n)

        # --- Assemble output ---
        rows = regions.to_dicts()
        for r in rows:
            tid = r["token_id"]
            mh = motif_hits.get(tid, [])
            r["top_motifs"] = [m for m, _ in mh]
            r["top_motif_pvalues"] = [p for _, p in mh]
            r["phastcons_mean"] = phast_means.get(tid)
            for tissue in chromhmm_urls:
                r[f"chromhmm_{tissue}"] = chromhmm_states_dict.get(tissue, {}).get(tid)

        out_df = pl.DataFrame(rows).select([
            "token_id", "top_motifs", "top_motif_pvalues", "phastcons_mean",
            *[f"chromhmm_{t}" for t in chromhmm_urls],
        ])
        annot_dir.mkdir(parents=True, exist_ok=True)
        out_df.write_parquet(out_path)
        outputs = [file_record(out_path, record_count=len(out_df))]

        metrics = {
            "n_regions": len(out_df),
            "n_with_motifs": int(sum(1 for r in rows if r["top_motifs"])),
            "n_with_phastcons": int(sum(1 for r in rows if r["phastcons_mean"] is not None
                                        and not (isinstance(r["phastcons_mean"], float) and r["phastcons_mean"] != r["phastcons_mean"]))),
            "tissues_with_chromhmm": list(chromhmm_urls.keys()),
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(f"\n  wrote {out_path.name}", file=sys.stderr)

    except Exception as e:
        write_summary(ctx, "error", inputs=inputs,
                      error={"class": type(e).__name__, "message": str(e),
                             "suggestion": "Check URLs in config; FIMO is the slowest step (~30 min); verify pymemesuite installed.",
                             "retryable": True})
        raise


if __name__ == "__main__":
    main()
