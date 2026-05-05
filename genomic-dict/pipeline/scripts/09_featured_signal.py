"""Stage 14 — featured-interval bigwig signal slices.

For each (featured file, featured interval) pair, sample the corresponding
bigwig at N evenly-spaced bins inside the interval and store the per-bin
mean fold-change-over-control. The output drives the continuous half of
the "continuous → peaks → tokens" toggle in the genomic-regions viz Step 1.

Per-file → bigwig mapping mirrors stage 03's `target_key`:
  ATAC-seq files       → <cell_line>__ATAC
  Histone ChIP-seq     → <cell_line>__<target>      (e.g., K562__H3K27ac)

Bigwig source resolution (per-file, in order):
  1. Local file at <REPO_ROOT>/<entry.local_path>
  2. Otherwise the entry.url (read remotely via byte-range; ~kb of network
     traffic per slice, no full download needed)

Reads:
  data/precomputed/featured_files.parquet         (from stage 10)
  data/precomputed/featured_intervals.parquet     (from stage 10)
  data/annotations/bigwig_manifest.yaml           (from stage 03)
Writes:
  data/precomputed/featured_signal.parquet        — one row per (file, interval)
  results/14_featured_signal/summary.json
"""
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import polars as pl
import pyBigWig
import requests
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, REPO_ROOT, file_record, stage_start, write_summary,
)


def manifest_key(cell_line: str, assay: str, target: str | None) -> str:
    """Map a featured file's metadata to its bigwig manifest key."""
    if assay == "ATAC-seq":
        return f"{cell_line}__ATAC"
    return f"{cell_line}__{target}"


def resolve_remote_url(url: str) -> str:
    """ENCODE's /@@download/ URLs 307-redirect to time-signed S3 URLs that
    pyBigWig's libcurl doesn't follow. The encode-public S3 bucket allows
    public anonymous reads, so we resolve the redirect once and strip the
    signature query string, leaving a stable unsigned S3 URL we can hand to
    pyBigWig for byte-range reads."""
    if "encode-public.s3" in url and "?" not in url:
        return url
    r = requests.head(url, allow_redirects=False, timeout=30)
    if r.status_code in (301, 302, 303, 307, 308):
        loc = r.headers.get("location", "")
        if loc:
            return loc.split("?", 1)[0]
    return url


def open_bigwig(entry: dict[str, Any]) -> tuple["pyBigWig.pyBigWig", str]:
    """Open from local file if present, else remote URL. Return (handle, source)."""
    local_rel = entry.get("local_path")
    if local_rel:
        local = REPO_ROOT / local_rel
        if local.exists():
            return pyBigWig.open(str(local)), f"local:{local_rel}"
    url = entry.get("url")
    if not url:
        raise RuntimeError(f"manifest entry has no local file and no url: {entry}")
    resolved = resolve_remote_url(url)
    return pyBigWig.open(resolved), f"remote:{resolved}"


def main() -> None:
    ctx = stage_start("09_featured_signal", __doc__)
    sc = ctx.stage_cfg

    n_bins = int(sc.get("n_bins", 200))
    paths = sc["paths"]
    precomp = PROJECT_ROOT / paths["precomputed_dir"]
    annot = PROJECT_ROOT / paths["annotations_dir"]

    files_pq = precomp / "featured_files.parquet"
    intervals_pq = precomp / "featured_intervals.parquet"
    manifest_yaml = annot / "bigwig_manifest.yaml"
    out_pq = precomp / "featured_signal.parquet"

    inputs = [
        file_record(files_pq),
        file_record(intervals_pq),
        file_record(manifest_yaml),
    ]

    try:
        for f in (files_pq, intervals_pq, manifest_yaml):
            if not f.exists():
                raise FileNotFoundError(f"missing input: {f}")

        manifest = yaml.safe_load(manifest_yaml.read_text())["entries"]
        files_df = pl.read_parquet(files_pq)
        intervals_df = pl.read_parquet(intervals_pq).sort("interval_id")

        rows: list[dict[str, Any]] = []
        per_key_metrics: dict[str, Any] = {}
        warnings: list[str] = []

        for f_row in files_df.iter_rows(named=True):
            key = manifest_key(f_row["cell_line"], f_row["assay"], f_row["target"])
            if key not in manifest:
                msg = f"no manifest entry for {key} (file_id={f_row['file_id']})"
                print(f"  WARN: {msg}", file=sys.stderr)
                warnings.append(msg)
                continue

            entry = manifest[key]
            t0 = time.time()
            try:
                bw, source = open_bigwig(entry)
            except Exception as e:
                msg = f"failed to open {key}: {e}"
                print(f"  WARN: {msg}", file=sys.stderr)
                warnings.append(msg)
                continue

            print(f"  {key} ({source})", file=sys.stderr)
            try:
                chroms = set(bw.chroms())
                for iv in intervals_df.iter_rows(named=True):
                    chrom = iv["chrom"]
                    start = int(iv["start"])
                    end = int(iv["end"])
                    bin_size = (end - start) / n_bins
                    positions = [int(start + (i + 0.5) * bin_size) for i in range(n_bins)]
                    if chrom not in chroms:
                        values = [float("nan")] * n_bins
                    else:
                        try:
                            raw = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
                            values = [float(v) if v is not None else float("nan") for v in raw]
                        except Exception as e:
                            msg = f"stats() failed on {key}/{iv['interval_id']}: {e}"
                            print(f"  WARN: {msg}", file=sys.stderr)
                            warnings.append(msg)
                            values = [float("nan")] * n_bins
                    rows.append({
                        "file_id": f_row["file_id"],
                        "interval_id": iv["interval_id"],
                        "n_bins": n_bins,
                        "positions": positions,
                        "values": values,
                    })
            finally:
                bw.close()

            per_key_metrics[key] = {
                "source": source,
                "duration_seconds": round(time.time() - t0, 2),
                "file_accession": entry.get("file_accession"),
            }

        if not rows:
            raise RuntimeError("no signal rows produced — check manifest + featured_files mapping")

        out_df = pl.DataFrame(
            rows,
            schema={
                "file_id": pl.Utf8,
                "interval_id": pl.Utf8,
                "n_bins": pl.Int32,
                "positions": pl.List(pl.Int64),
                "values": pl.List(pl.Float32),
            },
        )
        out_df.write_parquet(out_pq)

        write_summary(
            ctx,
            "success",
            inputs=inputs,
            outputs=[file_record(out_pq, record_count=len(rows))],
            metrics={
                "n_files": len(files_df),
                "n_intervals": len(intervals_df),
                "n_rows": len(rows),
                "n_bins": n_bins,
                "bigwigs": per_key_metrics,
            },
            warnings=warnings,
        )
    except Exception as e:
        write_summary(
            ctx,
            "error",
            error={
                "class": type(e).__name__,
                "message": str(e),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
