"""Stage 01 — curate the corpus manifest (~17k balanced-cap with stage-10 reserve).

Filters the 420k UMAP parquet to: 5 main assays + cell_line known. Left-joins
rich metadata (t1/t2). Then applies a balanced per-assay cap (default 4,000),
with a two-pass reserve so the (featured_cell_lines × featured_targets) cells
that stage 10 needs survive the cut even when individually rare.

Quality score = count of populated rich-meta fields (description, cell_type,
tissue, treatment, target, antibody). Used to break ties when sampling within
each assay — prefer files with more annotation.

Reads:
  data/hf_cache/<repo>/hg38_umap_all_3_13.parquet   (from stage 00)
  data/hf_cache/<repo>/hg38_meta_t{1,2}.parquet     (from stage 00)
Writes:
  data/corpus/manifest.parquet                       — filtered manifest
  results/01_curate_corpus/summary.json              — AI-ingestible summary
  results/01_curate_corpus/assay_distribution.md     — assay × cell_line table
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


# Fields that count toward the quality score. ATAC/DNase have no target/antibody;
# they just score lower on those, which is fine — we don't want to penalize them.
QUALITY_FIELDS = ("description", "cell_type", "tissue",
                  "treatment", "target", "antibody")


def has(c: str) -> pl.Expr:
    """Treat null, empty string, and 'UNKNOWN' as missing."""
    return (pl.col(c).is_not_null()
            & (pl.col(c) != "")
            & (pl.col(c) != "UNKNOWN"))


def apply_filter(
    df: pl.DataFrame,
    allowed_assays: list[str],
    unknown_marker: str,
) -> pl.DataFrame:
    """5-assay filter + cell_line non-empty AND ≠ UNKNOWN."""
    return (df.filter(pl.col("assay").is_in(allowed_assays))
              .filter(pl.col("cell_line").is_not_null()
                      & (pl.col("cell_line") != "")
                      & (pl.col("cell_line") != unknown_marker)))


def join_rich_metadata(
    manifest: pl.DataFrame, t1_path: Path, t2_path: Path,
) -> pl.DataFrame:
    """Left-join t1 and t2 — adds target/antibody/treatment/global_experiment_id."""
    t1 = pl.read_parquet(t1_path).select([
        "id", "target", "number_of_regions", "mean_region_width", "gc_content",
    ])
    t2 = pl.read_parquet(t2_path).select([
        "id", "treatment", "antibody", "bed_compliance", "data_format",
        "library_source", "global_sample_id", "global_experiment_id",
    ])
    return manifest.join(t1, on="id", how="left").join(t2, on="id", how="left")


def add_quality_score(df: pl.DataFrame) -> pl.DataFrame:
    """Add `quality_score` = count of populated QUALITY_FIELDS (0..6)."""
    score_expr = sum(has(f).cast(pl.Int8) for f in QUALITY_FIELDS)
    return df.with_columns(score_expr.alias("quality_score"))


def build_reserve_set(
    df: pl.DataFrame,
    cell_lines: list[str],
    targets: list[str | None],
    skip_cells: list[tuple[str, str]],
    n_per_cell: int,
) -> set[str]:
    """For each (cell_line, target) featured cell, reserve top-N by quality.

    `targets` may include None for the ATAC slot. `skip_cells` are
    (cell_line, target) pairs to drop from the grid (e.g. GM12878 H3K9me3
    has zero files in the source pool).
    """
    skip_set = {(cl, tgt) for cl, tgt in skip_cells}
    reserved: set[str] = set()
    for cl in cell_lines:
        for tgt in targets:
            if tgt is not None and (cl, tgt) in skip_set:
                continue
            if tgt is None:
                pool = df.filter((pl.col("cell_line") == cl)
                                 & (pl.col("assay") == "ATAC-seq"))
            else:
                pool = df.filter((pl.col("cell_line") == cl)
                                 & (pl.col("target") == tgt)
                                 & (pl.col("assay").str.contains("ChIP")))
            keep = (pool.sort(["quality_score", "id"], descending=[True, False])
                        .head(n_per_cell))
            reserved.update(keep["id"].to_list())
    return reserved


def balanced_cap_with_reserve(
    df: pl.DataFrame,
    per_assay_cap: int,
    reserved_ids: set[str],
) -> pl.DataFrame:
    """Per-assay top-quality cap, with reserved files always included.

    Within each assay: keep all reserved files; fill the remaining slots
    with the highest-quality non-reserved files (deterministic id tiebreak).
    Files already at or below the cap (e.g. Histone ChIP-seq) are kept entire.
    """
    parts = []
    reserved_list = list(reserved_ids)
    for assay in df["assay"].unique().to_list():
        pool = df.filter(pl.col("assay") == assay)
        in_reserve = pool.filter(pl.col("id").is_in(reserved_list))
        candidate = pool.filter(~pl.col("id").is_in(reserved_list))
        n_take = max(0, per_assay_cap - len(in_reserve))
        candidate = (candidate.sort(["quality_score", "id"],
                                    descending=[True, False])
                              .head(n_take))
        parts.append(pl.concat([in_reserve, candidate]))
    return pl.concat(parts)


def render_assay_breakdown_md(
    manifest: pl.DataFrame, stage10_coverage: list[dict],
) -> str:
    lines = ["# Corpus assay breakdown (post-filter)\n"]
    n = len(manifest)
    lines.append(f"Total surviving files: **{n:,}**\n")

    lines.append("## Assay distribution\n")
    lines.append("| assay | count | share |")
    lines.append("|---|---|---|")
    for row in (manifest.group_by("assay").agg(pl.len().alias("count"))
                .sort("count", descending=True).iter_rows(named=True)):
        lines.append(f"| `{row['assay']}` | {row['count']:,} | {row['count']/n:.1%} |")
    lines.append("")

    lines.append("## Top 20 cell lines\n")
    lines.append("| cell_line | count | share |")
    lines.append("|---|---|---|")
    for row in (manifest.group_by("cell_line").agg(pl.len().alias("count"))
                .sort("count", descending=True).head(20).iter_rows(named=True)):
        lines.append(f"| `{row['cell_line']}` | {row['count']:,} | {row['count']/n:.1%} |")
    lines.append("")

    lines.append("## Stage-10 cell coverage\n")
    lines.append("| cell_line | target | n |")
    lines.append("|---|---|---|")
    for row in stage10_coverage:
        marker = "✗" if row["n"] == 0 else "✓"
        lines.append(f"| {row['cell_line']} | {row['target']} | {marker} {row['n']} |")
    lines.append("")

    if "global_experiment_id" in manifest.columns:
        rich = manifest.filter(pl.col("global_experiment_id").is_not_null())
        lines.append("## Rich-metadata coverage\n")
        lines.append(f"- files with `global_experiment_id` non-null: **{len(rich):,}** "
                     f"/ {n:,} ({len(rich)/n:.1%})")
        lines.append(f"- distinct experiments: **{rich['global_experiment_id'].n_unique():,}**")

    return "\n".join(lines)


def stage10_coverage_report(
    manifest: pl.DataFrame,
    cell_lines: list[str],
    targets: list[str | None],
    skip_cells: list[tuple[str, str]],
) -> list[dict]:
    skip_set = {(cl, tgt) for cl, tgt in skip_cells}
    rows = []
    for cl in cell_lines:
        for tgt in targets:
            if tgt is not None and (cl, tgt) in skip_set:
                continue
            if tgt is None:
                hit = manifest.filter((pl.col("cell_line") == cl)
                                      & (pl.col("assay") == "ATAC-seq")).height
                rows.append({"cell_line": cl, "target": "ATAC", "n": int(hit)})
            else:
                hit = manifest.filter((pl.col("cell_line") == cl)
                                      & (pl.col("target") == tgt)
                                      & (pl.col("assay").str.contains("ChIP"))).height
                rows.append({"cell_line": cl, "target": tgt, "n": int(hit)})
    return rows


def main() -> None:
    ctx = stage_start("01_curate_corpus", __doc__)
    sc = ctx.stage_cfg

    hf_cache = PROJECT_ROOT / sc["paths"]["hf_cache_dir"]
    source_path = hf_cache / sc["source_parquet"]
    t1_path = hf_cache / sc["rich_metadata_t1"]
    t2_path = hf_cache / sc["rich_metadata_t2"]

    inputs = [file_record(source_path), file_record(t1_path), file_record(t2_path)]

    try:
        for p in [source_path, t1_path, t2_path]:
            if not p.exists():
                raise FileNotFoundError(f"{p} missing — run stage 00 first.")

        source = pl.read_parquet(source_path)
        n_input = len(source)

        manifest = apply_filter(
            source,
            allowed_assays=sc["allowed_assays"],
            unknown_marker=sc["cell_line_unknown_marker"],
        )
        n_after_filter = len(manifest)
        print(f"  filter: {n_input:,} → {n_after_filter:,}", file=sys.stderr)

        manifest = join_rich_metadata(manifest, t1_path, t2_path)
        manifest = add_quality_score(manifest)

        # Stage-10 reserve cells. Stage 01 reads stage 10's config block —
        # single source of truth for the featured grid.
        s10 = ctx.cfg.get("stages", {}).get("10_featured_narrative", {}) or {}
        cell_lines = s10.get("featured_cell_lines", [])
        targets = list(s10.get("featured_targets", []))
        skip_cells = [tuple(p) for p in s10.get("skip_cells", [])]
        reserve_per_cell = int(sc.get("reserve_per_featured_cell", 5))

        reserved = build_reserve_set(
            manifest, cell_lines, targets, skip_cells, reserve_per_cell,
        )
        print(f"  reserved {len(reserved)} files for stage-10 grid", file=sys.stderr)

        per_assay_cap = int(sc.get("per_assay_cap", 4000))
        manifest = balanced_cap_with_reserve(manifest, per_assay_cap, reserved)
        manifest = manifest.sort("id")
        n_final = len(manifest)
        print(f"  balanced cap: {n_after_filter:,} → {n_final:,}", file=sys.stderr)

        # Coverage check before writing.
        coverage = stage10_coverage_report(manifest, cell_lines, targets, skip_cells)
        missed = [r for r in coverage if r["n"] == 0]
        if missed:
            print(f"  WARN: {len(missed)} stage-10 cells have 0 files: {missed}",
                  file=sys.stderr)

        out_dir = PROJECT_ROOT / sc["paths"]["corpus_dir"]
        out_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = out_dir / "manifest.parquet"
        manifest.write_parquet(manifest_path)

        breakdown_path = ctx.results_dir / "assay_distribution.md"
        breakdown_path.write_text(render_assay_breakdown_md(manifest, coverage))

        outputs = [
            file_record(manifest_path, record_count=n_final),
            file_record(breakdown_path),
        ]

        assay_counts = {
            row["assay"]: int(row["count"])
            for row in (manifest.group_by("assay")
                        .agg(pl.len().alias("count"))
                        .iter_rows(named=True))
        }

        # Field-coverage metric (treats UNKNOWN/empty/null as missing).
        field_coverage = {}
        for c in ["target", "antibody", "treatment", "description",
                  "cell_type", "tissue", "global_experiment_id"]:
            if c in manifest.columns:
                field_coverage[c] = int(manifest.filter(has(c)).height)

        metrics = {
            "n_input_rows": n_input,
            "n_after_filter": n_after_filter,
            "n_reserved": len(reserved),
            "per_assay_cap": per_assay_cap,
            "n_final": n_final,
            "assay_distribution": assay_counts,
            "unique_cell_lines": int(manifest["cell_line"].n_unique()),
            "unique_targets": int(manifest.filter(has("target"))["target"].n_unique()),
            "field_coverage": field_coverage,
            "stage10_cell_coverage": coverage,
            "stage10_cells_missed": missed,
            "mean_quality_score": float(manifest["quality_score"].mean()),
        }

        write_summary(
            ctx,
            "success" if not missed else "partial",
            inputs=inputs, outputs=outputs, metrics=metrics,
            warnings=[f"stage-10 cell {r} has 0 files" for r in missed] if missed else [],
        )
        print(f"\n  wrote {manifest_path.name}: {n_final:,} files", file=sys.stderr)

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": "If source_parquet is missing, rerun stage 00. "
                              "Check stages.10_featured_narrative.featured_* config keys.",
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
