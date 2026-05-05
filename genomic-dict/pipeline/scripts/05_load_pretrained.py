"""Stage 07 — load pretrained R2V universe + embeddings, intersect with SCREEN.

Replaces the prior tokenize+train path. Instead of training a new model on an
85k-BED corpus, we use the production embedding `databio/r2v-encode-hg38`
directly. Each vocabulary region becomes a dictionary entry; SCREEN overlap
gives us a class annotation where available.

For pretrained regions that overlap multiple cCREs, we keep the cCRE with the
largest overlap length. Regions with no cCRE overlap are retained but carry
null class / accession — the viz hides these by default (we can't claim they
are regulatory; the pretrained corpus just saw them often enough to learn
an embedding).

Reads:
  databio/r2v-encode-hg38                        (HuggingFace; downloaded on first run)
  data/universe/ccre.metadata.parquet            (from stage 02, for SCREEN class tagging)
Writes:
  data/annotations/pretrained_universe.parquet   — one row per pretrained region
  results/07_load_pretrained/summary.json        — AI-ingestible summary
"""
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def parse_region_string(s: str) -> tuple[str, int, int] | None:
    """Parse `chr:start-end` into (chrom, start, end). Return None for specials like '<unk>'."""
    if s.startswith("<") and s.endswith(">"):
        return None
    if ":" not in s or "-" not in s:
        return None
    try:
        chrom, coords = s.split(":", 1)
        start_s, end_s = coords.split("-", 1)
        return chrom, int(start_s), int(end_s)
    except (ValueError, TypeError):
        return None


def intersect_with_screen(
    pretrained_df: pl.DataFrame, screen_df: pl.DataFrame,
) -> pl.DataFrame:
    """Left-join SCREEN class annotations onto pretrained regions by overlap.
    For multi-overlap, keep the cCRE with the largest overlap length."""
    import pyranges as pr

    def to_pyranges(df: pl.DataFrame, extras: list[str]) -> "pr.PyRanges":
        pdf = df.select(["chrom", "start", "end", *extras]).to_pandas()
        pdf = pdf.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
        return pr.PyRanges(pdf)

    p_pr = to_pyranges(pretrained_df, ["token_id"])
    s_pr = to_pyranges(screen_df, ["accession_hex", "cclass"])

    # Intersect gives us overlap intervals with metadata from both sides.
    inter = p_pr.intersect(s_pr, how="first")  # one row per pretrained-cCRE overlap
    # `intersect` with how='first' returns one record per pretrained region with a cCRE
    # match; if multiple cCREs overlap, PyRanges picks one. We refine to largest-overlap
    # manually via a full join.

    full = p_pr.join(s_pr, suffix="_screen")
    if len(full) == 0:
        # No overlaps at all — shouldn't happen with real data.
        return pretrained_df.with_columns([
            pl.lit(None, dtype=pl.Utf8).alias("accession_hex"),
            pl.lit(None, dtype=pl.Utf8).alias("cclass"),
            pl.lit(False).alias("overlaps_screen"),
        ])

    # join returns a DataFrame (via .df) with columns for both sides.
    joined_pdf = full.df
    # Compute overlap length per join row.
    joined_pdf["overlap_len"] = (
        joined_pdf[["End", "End_screen"]].min(axis=1)
        - joined_pdf[["Start", "Start_screen"]].max(axis=1)
    ).clip(lower=0)

    # Keep the cCRE match with the largest overlap per pretrained token_id.
    best = (
        joined_pdf.sort_values(["token_id", "overlap_len"], ascending=[True, False])
        .drop_duplicates(subset=["token_id"], keep="first")
        [["token_id", "accession_hex", "cclass"]]
    )
    best_pl = pl.from_pandas(best).with_columns(pl.lit(True).alias("overlaps_screen"))

    # Left-join back to the full pretrained table so non-overlapping tokens get nulls.
    out = pretrained_df.join(best_pl, on="token_id", how="left")
    out = out.with_columns(
        pl.col("overlaps_screen").fill_null(False),
    )
    return out


def main() -> None:
    ctx = stage_start("05_load_pretrained", __doc__)
    sc = ctx.stage_cfg

    hf_model = sc["hf_model"]
    universe_dir = PROJECT_ROOT / sc["paths"]["universe_dir"]
    annot_dir = PROJECT_ROOT / sc["paths"]["annotations_dir"]
    metadata_path = universe_dir / "ccre.metadata.parquet"
    out_path = annot_dir / "pretrained_universe.parquet"

    inputs = [file_record(metadata_path)]

    try:
        if not metadata_path.exists():
            raise FileNotFoundError(
                f"{metadata_path} missing (run stage 02 first)."
            )

        from geniml.region2vec.main import Region2VecExModel

        print(f"  loading pretrained model {hf_model}...", file=sys.stderr)
        t0 = time.time()
        model = Region2VecExModel(model_path=hf_model)
        load_time = time.time() - t0
        print(f"    loaded in {load_time:.1f}s", file=sys.stderr)

        # Extract vocab {region_str: token_id}
        vocab = dict(model.tokenizer.get_vocab())
        print(f"  vocab size: {len(vocab):,}", file=sys.stderr)

        # Parse regions + pull embeddings. The pretrained R2V is a torch nn.Embedding;
        # embeddings live in model.model.projection.weight — a (vocab_size, dim) tensor.
        # Lookup by token_id row index, not by region-string key.
        print("  parsing regions and pulling embeddings...", file=sys.stderr)
        weight = model.model.projection.weight.detach().cpu().numpy()  # (vocab, dim)
        vocab_rows = weight.shape[0]
        embedding_dim = weight.shape[1]
        print(f"    embedding matrix: {weight.shape}", file=sys.stderr)

        rows: list[dict[str, Any]] = []
        skipped_special = 0
        skipped_unparseable = 0
        skipped_no_vector = 0

        for region_str, tok_id in vocab.items():
            parsed = parse_region_string(region_str)
            if parsed is None:
                if region_str.startswith("<") and region_str.endswith(">"):
                    skipped_special += 1
                else:
                    skipped_unparseable += 1
                continue
            chrom, start, end = parsed
            tid = int(tok_id)
            if tid < 0 or tid >= vocab_rows:
                skipped_no_vector += 1
                continue
            vec = weight[tid]
            rows.append({
                "token_id": tid,
                "region": region_str,
                "chrom": chrom,
                "start": start,
                "end": end,
                "embedding": vec.astype(np.float32).tolist(),
            })

        print(
            f"    {len(rows):,} regions with embeddings; "
            f"skipped {skipped_special} special, "
            f"{skipped_unparseable} unparseable, "
            f"{skipped_no_vector} missing vector",
            file=sys.stderr,
        )

        pretrained_df = pl.DataFrame(rows)

        # Intersect with SCREEN metadata for class tagging.
        print("  intersecting with SCREEN cCREs...", file=sys.stderr)
        screen_df = pl.read_parquet(metadata_path).select([
            "chrom", "start", "end", "accession_hex", "cclass",
        ])
        tagged = intersect_with_screen(pretrained_df, screen_df)

        tagged = tagged.sort(["chrom", "start", "end"])

        annot_dir.mkdir(parents=True, exist_ok=True)
        tagged.write_parquet(out_path)
        outputs = [file_record(out_path, record_count=len(tagged))]

        # Per-class summary.
        class_counts = (
            tagged.group_by("cclass").agg(pl.len().alias("n"))
            .sort("n", descending=True).to_dicts()
        )
        per_chrom = (
            tagged.group_by("chrom").agg(pl.len().alias("n"))
            .sort("chrom").to_dicts()
        )

        n_total = len(tagged)
        n_with_screen = int(tagged.filter(pl.col("overlaps_screen")).height)

        metrics = {
            "hf_model": hf_model,
            "vocab_size_raw": len(vocab),
            "n_regions_parsed": len(rows),
            "n_skipped_special": skipped_special,
            "n_skipped_unparseable": skipped_unparseable,
            "n_skipped_no_vector": skipped_no_vector,
            "n_total": n_total,
            "n_with_screen_overlap": n_with_screen,
            "screen_overlap_fraction": round(n_with_screen / max(1, n_total), 4),
            "per_class_counts": class_counts,
            "per_chromosome_counts": per_chrom,
            "load_seconds": round(load_time, 1),
            "embedding_dim": embedding_dim,
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(
            f"\n  wrote {out_path.name}: {n_total:,} regions, "
            f"{n_with_screen:,} with SCREEN class ({n_with_screen/max(1,n_total):.1%})",
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
                    "If HF download failed, check network. "
                    "If metadata missing, rerun stage 02. "
                    "Stage is idempotent."
                ),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
