"""Stage 00 — inspect HuggingFace metadata parquets from databio/bedbase-umap.

Downloads hg38_meta_t1.parquet and hg38_meta_t2.parquet, reports column schema
and value distributions for string-valued columns, and writes a human-readable
columns.md plus a machine-readable summary.json.

Purpose: surface what fields are available so we can pick filter criteria for
stage 01 (corpus curation). Does not filter anything itself.

Reads:  (nothing upstream — pulls from HuggingFace Hub)
Writes:
  data/hf_cache/<repo>/<filename>           — downloaded parquets (gitignored)
  results/00_inspect_metadata/columns.md    — human-readable column summary
  results/00_inspect_metadata/summary.json  — AI-ingestible summary
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
from huggingface_hub import hf_hub_download

# Make sibling _common importable when run as `python 00_inspect_metadata.py`.
sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def download_parquets(repo_id: str, filenames: list[str], cache_dir: Path) -> list[Path]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for filename in filenames:
        print(f"  fetching {filename}...", file=sys.stderr)
        local = hf_hub_download(
            repo_id=repo_id,
            filename=filename,
            local_dir=str(cache_dir),
        )
        paths.append(Path(local))
    return paths


def summarize_column(s: pl.Series, max_unique: int, top_n: int) -> dict:
    """Column-level schema + nullity + top values (for string-ish columns)."""
    null_count = int(s.null_count())
    n_rows = len(s)
    info: dict = {
        "dtype": str(s.dtype),
        "null_count": null_count,
        "null_fraction": round(null_count / n_rows, 4) if n_rows else None,
    }
    # Value-count top-N for string columns with bounded cardinality.
    if s.dtype == pl.Utf8:
        n_unique = int(s.n_unique())
        info["n_unique"] = n_unique
        if 0 < n_unique <= max_unique:
            vc = s.value_counts(sort=True).head(top_n)
            # polars value_counts schema is [col_name, count]
            col_name = s.name
            info["top_values"] = [
                {"value": row[col_name], "count": int(row["count"])}
                for row in vc.iter_rows(named=True)
            ]
        else:
            info["top_values"] = None
    elif s.dtype.is_numeric():
        # Cheap numeric summary.
        non_null = s.drop_nulls()
        if len(non_null) > 0:
            info["min"] = float(non_null.min())
            info["max"] = float(non_null.max())
            info["mean"] = float(non_null.mean())
    return info


def summarize_file(path: Path, max_unique: int, top_n: int) -> dict:
    df = pl.read_parquet(path)
    columns = {col: summarize_column(df[col], max_unique, top_n) for col in df.columns}
    return {
        "n_rows": int(len(df)),
        "n_cols": len(df.columns),
        "columns": columns,
    }


def render_columns_md(per_file_summary: dict[str, dict], repo_id: str) -> str:
    lines: list[str] = [
        f"# Metadata columns — `{repo_id}`",
        "",
        "Column schema + value distributions for each parquet. Use to pick "
        "corpus filter criteria for stage 01.",
        "",
    ]
    for fname, info in per_file_summary.items():
        lines += [
            f"## `{fname}`",
            "",
            f"- rows: {info['n_rows']:,}",
            f"- cols: {info['n_cols']}",
            "",
            "| column | dtype | null % | n_unique | top values |",
            "|---|---|---|---|---|",
        ]
        for col, ci in info["columns"].items():
            null_pct = (
                f"{ci['null_fraction'] * 100:.1f}%"
                if ci.get("null_fraction") is not None else "—"
            )
            n_uniq = ci.get("n_unique", "")
            tv = ci.get("top_values") or []
            top_str = (
                ", ".join(f"`{v['value']}`={v['count']}" for v in tv[:5])
                if tv else ""
            )
            lines.append(f"| `{col}` | `{ci['dtype']}` | {null_pct} | {n_uniq} | {top_str} |")
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    ctx = stage_start("00_inspect_metadata", __doc__)
    sc = ctx.stage_cfg

    hf_repo = sc["hf_repo"]
    hf_files = sc["hf_files"]
    max_unique = int(sc.get("top_values_max_unique", 500))
    top_n = int(sc.get("top_values_n", 10))

    cache_dir = PROJECT_ROOT / sc["paths"]["hf_cache_dir"] / hf_repo.replace("/", "__")

    try:
        parquet_paths = download_parquets(hf_repo, hf_files, cache_dir)

        per_file: dict[str, dict] = {}
        outputs: list[dict] = []
        for p in parquet_paths:
            per_file[p.name] = summarize_file(p, max_unique=max_unique, top_n=top_n)
            outputs.append(file_record(p, record_count=per_file[p.name]["n_rows"]))

        md_path = ctx.results_dir / "columns.md"
        md_path.write_text(render_columns_md(per_file, hf_repo))
        outputs.append(file_record(md_path))

        # Flat metrics for summary.json — per-file row/col counts plus full column info.
        metrics: dict = {
            "per_file": {
                fname: {"n_rows": info["n_rows"], "n_cols": info["n_cols"]}
                for fname, info in per_file.items()
            },
            "total_rows": sum(info["n_rows"] for info in per_file.values()),
            "columns": per_file,
        }

        write_summary(ctx, "success", outputs=outputs, metrics=metrics)
        print(f"\n  rows: {metrics['total_rows']:,} across {len(per_file)} parquet(s)", file=sys.stderr)

    except Exception as e:
        write_summary(
            ctx, "error",
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": (
                    "Verify HuggingFace access (no auth needed for databio/bedbase-umap) "
                    "and network connectivity. Stage is safe to retry."
                ),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
