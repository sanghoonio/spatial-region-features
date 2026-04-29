"""Stage 04 — bigwig-based intrinsic features per focus_chromosome region.

Reads the chr16 subset of pretrained_universe.parquet (from stage 07) and for
each region computes:

  * **Histone mark means**: 15 scalars (5 marks × 3 tissues) — mean
    fold-change-over-control inside the interval.
  * **ATAC thumbnail**: 3 vectors (one per tissue) of 100-bin mean
    coverage, for the Panel 1 entry card.

Reads:
  data/annotations/pretrained_universe.parquet         (from stage 07)
  data/annotations/bigwig_manifest.yaml                (from stage 03)
  data/annotations/bigwigs/*.bigWig                    (from stage 03)
Writes:
  data/annotations/intrinsic_bigwig.parquet            — one row per region
  results/04_extract_intrinsic/summary.json            — AI-ingestible summary

NOTE: stage 07 must run before this stage. The pretrained universe is the
vocabulary; this stage annotates each entry with tissue-specific bigwig signal.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
import pyBigWig
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def load_focus_universe(pretrained_path: Path, focus_chrom: str) -> pl.DataFrame:
    """Return the focus-chromosome subset of the pretrained universe."""
    df = pl.read_parquet(pretrained_path).filter(pl.col("chrom") == focus_chrom)
    # Keep coordinates + identity; drop embedding here (viz loads it separately).
    return df.select([
        "token_id", "region", "chrom", "start", "end",
        "overlaps_screen", "accession_hex", "cclass",
    ])


def load_manifest(path: Path) -> dict[str, dict]:
    data = yaml.safe_load(path.read_text()) or {}
    return data.get("entries", {}) or {}


def bigwig_path(entry: dict) -> Path:
    """Resolve entry.local_path (repo-relative) against REPO_ROOT."""
    from _common import REPO_ROOT
    return REPO_ROOT / entry["local_path"]


def compute_mean_column(bw: "pyBigWig.pyBigWig", df: pl.DataFrame) -> np.ndarray:
    """Mean bigwig signal per cCRE interval. Returns np.array of len(df), NaN for missing."""
    out = np.empty(len(df), dtype=np.float32)
    chroms = df["chrom"].to_list()
    starts = df["start"].to_list()
    ends = df["end"].to_list()
    for i, (c, s, e) in enumerate(zip(chroms, starts, ends)):
        try:
            v = bw.stats(c, int(s), int(e), type="mean", nBins=1)[0]
            out[i] = float(v) if v is not None else np.nan
        except Exception:
            out[i] = np.nan
    return out


def compute_thumbnail_column(
    bw: "pyBigWig.pyBigWig", df: pl.DataFrame, n_bins: int,
) -> list[list[float]]:
    """Per-cCRE n_bins-vector of mean signal. Returns list of lists (ragged-safe)."""
    out: list[list[float]] = []
    chroms = df["chrom"].to_list()
    starts = df["start"].to_list()
    ends = df["end"].to_list()
    for c, s, e in zip(chroms, starts, ends):
        try:
            vals = bw.stats(c, int(s), int(e), type="mean", nBins=n_bins)
            # Replace None with NaN so polars stores a uniform list[f32].
            out.append([float(v) if v is not None else float("nan") for v in vals])
        except Exception:
            out.append([float("nan")] * n_bins)
    return out


def column_name(entry: dict, feature: str) -> str:
    """Produce a stable column name from (biosample, target/ATAC) + feature type."""
    tag = entry.get("target") or "ATAC"
    return f"{entry['biosample']}__{tag}__{feature}"


def main() -> None:
    ctx = stage_start("04_extract_intrinsic", __doc__)
    sc = ctx.stage_cfg

    chrom = sc["focus_chromosome"]
    n_bins = int(sc.get("atac_thumbnail_bins", 100))

    annot_dir = PROJECT_ROOT / sc["paths"]["annotations_dir"]
    pretrained_path = annot_dir / "pretrained_universe.parquet"
    manifest_path = annot_dir / "bigwig_manifest.yaml"
    out_path = annot_dir / "intrinsic_bigwig.parquet"

    inputs = [file_record(pretrained_path), file_record(manifest_path)]

    try:
        if not pretrained_path.exists():
            raise FileNotFoundError(
                f"{pretrained_path} missing. Run stage 07 first to load the pretrained universe."
            )
        if not manifest_path.exists():
            raise FileNotFoundError(
                f"{manifest_path} missing. Run stage 03 first to fetch bigwigs."
            )

        ccres = load_focus_universe(pretrained_path, chrom).sort(["start", "end"])
        manifest = load_manifest(manifest_path)

        n_ccres = len(ccres)
        print(f"  {n_ccres:,} {chrom} cCREs; {len(manifest)} bigwigs to read", file=sys.stderr)

        # Collect columns into a dict; attach to the BED schema at the end.
        feature_columns: dict[str, Any] = {}
        per_bigwig_metrics: dict[str, dict] = {}

        for key in sorted(manifest.keys()):
            entry = manifest[key]
            bw_path = bigwig_path(entry)
            if not bw_path.exists():
                raise FileNotFoundError(f"bigwig {bw_path} missing for {key}")

            print(f"  reading {key}...", file=sys.stderr)
            t0 = time.time()
            bw = pyBigWig.open(str(bw_path))
            try:
                if chrom not in bw.chroms():
                    raise RuntimeError(f"{chrom} not in {key}'s chromosome list")

                mean_col = column_name(entry, "mean")
                feature_columns[mean_col] = compute_mean_column(bw, ccres)

                # ATAC gets an additional thumbnail column; histone marks don't.
                if entry["assay"] == "ATAC-seq":
                    thumb_col = column_name(entry, "thumb")
                    feature_columns[thumb_col] = compute_thumbnail_column(bw, ccres, n_bins)
            finally:
                bw.close()

            per_bigwig_metrics[key] = {
                "file_accession": entry.get("file_accession"),
                "duration_seconds": round(time.time() - t0, 2),
                "has_thumbnail": entry["assay"] == "ATAC-seq",
            }

        # Assemble wide parquet.
        out_df = ccres.with_columns([
            pl.Series(name=col, values=vals)
            for col, vals in feature_columns.items()
        ])

        annot_dir.mkdir(parents=True, exist_ok=True)
        out_df.write_parquet(out_path)
        outputs = [file_record(out_path, record_count=len(out_df))]

        # Surface a flat NaN-rate per feature column for quick sanity.
        nan_rates = {}
        for col in feature_columns:
            if isinstance(feature_columns[col], np.ndarray):
                nan_rates[col] = float(np.mean(np.isnan(feature_columns[col])))
            # thumbnails: skip detailed nan analysis here — just count rows with all-NaN
            else:
                all_nan = sum(1 for v in feature_columns[col] if all(np.isnan(x) for x in v))
                nan_rates[col] = round(all_nan / len(ccres), 4)

        metrics = {
            "n_ccres": n_ccres,
            "n_bigwigs_read": len(manifest),
            "n_mean_columns": sum(1 for k in feature_columns if k.endswith("__mean")),
            "n_thumbnail_columns": sum(1 for k in feature_columns if k.endswith("__thumb")),
            "atac_thumbnail_bins": n_bins,
            "nan_rates": nan_rates,
            "per_bigwig": per_bigwig_metrics,
        }

        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(
            f"\n  wrote {out_path.name} ({len(out_df):,} rows × {len(out_df.columns)} cols)",
            file=sys.stderr,
        )

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": (
                    "Check that stages 02 and 03 completed. "
                    "If focus_chromosome was changed, rerun stage 02 first."
                ),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
