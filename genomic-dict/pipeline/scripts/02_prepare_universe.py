"""Stage 02 — prepare SCREEN V4 cCRE universe.

Downloads the genome-wide SCREEN V4 cCRE registry (GRCh38) from wenglab,
filters to the configured class set, sorts deterministically, and writes:

  - BED3 files (chrom/start/end only) — gtars-compatible, used by stage 07
    tokenizer. Gtars' universe parser rejects 6-column BEDs unless col 5 is
    numeric (it is not for SCREEN; col 5 is an EH38E* accession).
  - A metadata parquet carrying the labels (accession_hex, accession_5group,
    cclass) for downstream stages that need annotations. Join on chrom/start/end.

V4 BED source format: chrom, start, end, accession_hex (EH38D*),
accession_5group (EH38E*), class_label.

V4 classes: PLS, pELS, dELS, CA-CTCF, CA-H3K4me3, CA-TF, CA, TF.
We keep the 5 that map to the original plan's taxonomy by default (see
config.yaml stages.02_prepare_universe.allowed_classes).

Reads:  (nothing upstream — pulls from wenglab)
Writes:
  data/universe/GRCh38-cCREs.v4.bed.gz           — raw download (gzipped)
  data/universe/ccre.all.bed                      — BED3, genome-wide, gtars-ready
  data/universe/ccre.<chrom>.bed                  — BED3, focus chromosome
  data/universe/ccre.metadata.parquet             — labels; join on chrom/start/end
  results/02_prepare_universe/class_breakdown.md  — per-class + per-chrom counts
  results/02_prepare_universe/summary.json        — AI-ingestible summary
"""
from __future__ import annotations

import gzip
import shutil
import sys
from pathlib import Path

import polars as pl
import requests

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def download_gz(url: str, dest: Path) -> Path:
    """Download url to dest (gzipped). Idempotent: skip if dest exists."""
    if dest.exists():
        print(f"  {dest.name} already present, skipping download", file=sys.stderr)
        return dest
    print(f"  downloading {url}...", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=(10, 300)) as resp:
        resp.raise_for_status()
        with open(tmp, "wb") as out:
            with gzip.open(out, "wb") as gz:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    if chunk:
                        gz.write(chunk)
    tmp.rename(dest)
    return dest


def load_ccre_bed(path: Path) -> pl.DataFrame:
    """Read a SCREEN V4 cCRE BED (gzipped or plain)."""
    return pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "accession_hex", "accession_5group", "cclass"],
        schema_overrides={
            "chrom": pl.Utf8, "start": pl.Int64, "end": pl.Int64,
            "accession_hex": pl.Utf8, "accession_5group": pl.Utf8, "cclass": pl.Utf8,
        },
    )


def write_bed3_sorted(df: pl.DataFrame, path: Path) -> None:
    """Write chrom/start/end only, sorted. BED3 is the only format gtars'
    universe parser accepts without requiring a numeric column-5 score."""
    out = df.sort(["chrom", "start", "end"]).select(["chrom", "start", "end"])
    out.write_csv(path, separator="\t", include_header=False)


def render_breakdown_md(df: pl.DataFrame, focus_chrom: str) -> str:
    lines = ["# cCRE universe breakdown\n"]
    total = len(df)
    lines.append(f"Total cCREs after class filter: **{total:,}**\n")

    lines.append("## Per-class counts (genome-wide)\n")
    lines.append("| class | count | share |")
    lines.append("|---|---|---|")
    for row in df.group_by("cclass").agg(pl.len().alias("n")).sort("n", descending=True).iter_rows(named=True):
        lines.append(f"| `{row['cclass']}` | {row['n']:,} | {row['n']/total:.1%} |")
    lines.append("")

    focus_df = df.filter(pl.col("chrom") == focus_chrom)
    lines.append(f"## Focus chromosome ({focus_chrom})\n")
    lines.append(f"- total on {focus_chrom}: **{len(focus_df):,}**")
    lines.append("")
    lines.append(f"### Per-class on {focus_chrom}\n")
    lines.append("| class | count |")
    lines.append("|---|---|")
    for row in focus_df.group_by("cclass").agg(pl.len().alias("n")).sort("n", descending=True).iter_rows(named=True):
        lines.append(f"| `{row['cclass']}` | {row['n']:,} |")
    lines.append("")

    lines.append("## Per-chromosome totals\n")
    lines.append("| chrom | count |")
    lines.append("|---|---|")
    for row in df.group_by("chrom").agg(pl.len().alias("n")).sort("chrom").iter_rows(named=True):
        lines.append(f"| `{row['chrom']}` | {row['n']:,} |")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    ctx = stage_start("02_prepare_universe", __doc__)
    sc = ctx.stage_cfg

    universe_dir = PROJECT_ROOT / sc["paths"]["universe_dir"]
    universe_dir.mkdir(parents=True, exist_ok=True)

    raw_path = universe_dir / "GRCh38-cCREs.v4.bed.gz"
    all_path = universe_dir / "ccre.all.bed"
    focus_chrom = sc["focus_chromosome"]
    focus_path = universe_dir / f"ccre.{focus_chrom}.bed"
    allowed = list(sc["allowed_classes"])

    try:
        download_gz(sc["source_url"], raw_path)
        inputs = [file_record(raw_path)]

        # Load raw registry.
        print("  parsing cCRE registry...", file=sys.stderr)
        all_ccres = load_ccre_bed(raw_path)
        n_raw = len(all_ccres)

        # Capture raw-class distribution before filtering, to report what we dropped.
        raw_class_counts = {
            row["cclass"]: int(row["n"])
            for row in all_ccres.group_by("cclass").agg(pl.len().alias("n")).iter_rows(named=True)
        }

        # Filter.
        kept = all_ccres.filter(pl.col("cclass").is_in(allowed))
        n_kept = len(kept)

        # Write genome-wide + focus-chrom BED3 (gtars-ready).
        write_bed3_sorted(kept, all_path)
        focus_df = kept.filter(pl.col("chrom") == focus_chrom)
        write_bed3_sorted(focus_df, focus_path)

        # Sidecar metadata parquet — same rows, all columns. Stages that need
        # labels (stage 04 etc) join on (chrom, start, end).
        metadata_path = universe_dir / "ccre.metadata.parquet"
        kept.sort(["chrom", "start", "end"]).write_parquet(metadata_path)

        # Human-readable breakdown.
        breakdown_path = ctx.results_dir / "class_breakdown.md"
        breakdown_path.write_text(render_breakdown_md(kept, focus_chrom))

        outputs = [
            file_record(all_path, record_count=n_kept),
            file_record(focus_path, record_count=len(focus_df)),
            file_record(metadata_path, record_count=n_kept),
            file_record(breakdown_path),
        ]

        kept_class_counts = {
            row["cclass"]: int(row["n"])
            for row in kept.group_by("cclass").agg(pl.len().alias("n")).iter_rows(named=True)
        }
        dropped = {k: v for k, v in raw_class_counts.items() if k not in allowed}

        metrics = {
            "n_ccres_raw": n_raw,
            "n_ccres_kept": n_kept,
            "n_ccres_focus": len(focus_df),
            "focus_chromosome": focus_chrom,
            "raw_class_counts": raw_class_counts,
            "kept_class_counts": kept_class_counts,
            "dropped_classes": dropped,
            "allowed_classes": allowed,
        }

        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(
            f"\n  kept {n_kept:,}/{n_raw:,} cCREs ({n_kept/n_raw:.1%}); "
            f"{focus_chrom}: {len(focus_df):,}",
            file=sys.stderr,
        )
        if dropped:
            print(f"  dropped V4 classes: {', '.join(sorted(dropped))}", file=sys.stderr)

    except Exception as e:
        write_summary(
            ctx, "error",
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": (
                    "If download failed, check network and URL in "
                    "config.yaml stages.02_prepare_universe.source_url. "
                    "Safe to retry; download is idempotent."
                ),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
