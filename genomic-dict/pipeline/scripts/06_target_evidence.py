"""Stage 06 — ENCODE cCRE V4 → target gene evidence (Tier A regulatory targets).

Replaces the previously-planned nearest-gene work (which has the well-known
"closest gene is wrong" failure mode, exemplified by FTO's chr16 enhancer
acting on IRX3 ~1 Mb away). Instead, integrates ENCODE Registry V4
cCRE-Gene Links curated from direct experimental evidence:

  - 3D Chromatin contacts (Hi-C / ChIA-PET / HiChIP / pCHi-C)
  - eQTL associations (variant → gene expression links per tissue)

CRISPR perturbation evidence is also in the V4 release but contains zero
chr16 entries (CRISPRi tile-screens have so far covered chr3/10/12/19/X
only), so this stage skips it for chr16-scoped output.

Approach:
  1. Map each viz_chr16 token → current V4 EH38E accession via coordinate
     match (94.5% exact; remainder via highest-overlap-fraction fallback).
  2. Filter the chr16-pre-extracted Gene-Links TSVs to rows whose EH38E
     accession is in the mapped set.
  3. Aggregate per (token_id, target_gene, evidence_type) with biosample /
     tissue, score, and p-value retained.
  4. Output a long-format parquet usable by the dictionary card.

Reads:
  data/precomputed/viz_chr16.parquet                     (chr16 token coords)
  data/universe/screen_v4_2024-07.bed                    (V4 master cCRE list)
  data/annotations/encode_gene_links/3d_chromatin_chr16.tsv  (chr16-filtered)
  data/annotations/encode_gene_links/eqtl_chr16.tsv          (chr16-filtered)
Writes:
  data/precomputed/region_target_evidence.parquet        (one row per token×gene×evidence)
  data/precomputed/region_target_evidence_summary.parquet (per-token aggregate)
  results/06_target_evidence/summary.json

Schema (region_target_evidence.parquet):
  token_id            Int64    chr16 viz token
  evidence_type       String   "3d_chromatin" or "eqtl"
  target_gene_id      String   Ensembl gene ID
  target_gene_name    String   HGNC symbol (cleaned)
  gene_type           String   protein_coding / lncRNA / etc.
  source              String   assay (e.g., "Hi-C", "ChIA-PET") for 3D; "GTEx_v8" for eQTL
  biosample_or_tissue String   cell-line / tissue context
  experiment_id       String   ENCODE experiment ID (3D) or NULL (eQTL)
  variant_id          String   eQTL variant (NULL for 3D)
  score               Float32  contact score (3D) or |slope| (eQTL)
  pvalue              Float64  raw p-value
  ccre_accession_v4   String   EH38E... in V4 registry
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


def coord_map_viz_to_v4(
    viz: pl.DataFrame,
    v4: pl.DataFrame,
    min_overlap_frac: float = 0.5,
) -> pl.DataFrame:
    """For each viz token, find the corresponding V4 EH38E accession.

    Strategy:
      1. Exact match on (chrom, start, end) — gets ~94% of tokens.
      2. For unmatched, fallback to highest-overlap V4 cCRE; require overlap
         fraction >= min_overlap_frac.

    Returns DataFrame with columns: token_id, ccre_accession_v4 (nullable).
    """
    # Exact match
    exact = viz.join(
        v4.select(["chrom", "start", "end", "accession"]),
        on=["chrom", "start", "end"], how="left",
    ).rename({"accession": "ccre_accession_v4"})

    n_exact = exact["ccre_accession_v4"].is_not_null().sum()
    print(f"  exact coord match: {n_exact:,}/{len(viz):,} = {n_exact/len(viz):.1%}",
          file=sys.stderr)

    # Overlap fallback for unmatched
    unmatched = exact.filter(pl.col("ccre_accession_v4").is_null()).select(
        ["token_id", "chrom", "start", "end"]
    )
    if len(unmatched) == 0:
        return exact.select(["token_id", "ccre_accession_v4"])

    print(f"  resolving {len(unmatched):,} unmatched via overlap fallback...",
          file=sys.stderr)

    # Build per-chrom V4 dicts for overlap lookup
    v4_arrays: dict[str, tuple[np.ndarray, np.ndarray, list[str]]] = {}
    for chrom in unmatched["chrom"].unique():
        sub = v4.filter(pl.col("chrom") == chrom).sort("start")
        v4_arrays[chrom] = (
            sub["start"].to_numpy(),
            sub["end"].to_numpy(),
            sub["accession"].to_list(),
        )

    fallback_acc: list[str | None] = []
    for r in unmatched.iter_rows(named=True):
        starts, ends, accs = v4_arrays[r["chrom"]]
        # Overlap = min(end, end_v4) - max(start, start_v4); positive iff overlap exists
        ov_start = np.maximum(r["start"], starts)
        ov_end = np.minimum(r["end"], ends)
        ov_len = np.maximum(0, ov_end - ov_start)
        viz_len = r["end"] - r["start"]
        if viz_len <= 0:
            fallback_acc.append(None)
            continue
        ov_frac = ov_len / viz_len
        best = int(np.argmax(ov_frac))
        if ov_frac[best] >= min_overlap_frac:
            fallback_acc.append(accs[best])
        else:
            fallback_acc.append(None)

    fallback_df = unmatched.with_columns(
        pl.Series("ccre_accession_v4", fallback_acc).cast(pl.Utf8)
    ).select(["token_id", "ccre_accession_v4"])

    # Combine: exact-matched + fallback
    matched_exact = exact.filter(pl.col("ccre_accession_v4").is_not_null()).select(
        ["token_id", "ccre_accession_v4"]
    )
    combined = pl.concat([matched_exact, fallback_df]).sort("token_id")

    n_total = combined["ccre_accession_v4"].is_not_null().sum()
    print(f"  total mapped (exact + overlap): {n_total:,}/{len(viz):,} = "
          f"{n_total/len(viz):.1%}", file=sys.stderr)
    return combined


def load_3d_chromatin(path: Path) -> pl.DataFrame:
    """3D-Chromatin schema: cCRE | Gene ID | Gene Name | Gene Type | Assay | Exp ID | Biosample | Score | P-value

    Pvalue is NA for some RNAPII-ChIAPET rows; we want it nullable Float64."""
    df = pl.read_csv(
        path, separator="\t", has_header=False,
        new_columns=[
            "ccre_accession_v4", "target_gene_id", "target_gene_name",
            "gene_type", "source", "experiment_id", "biosample_or_tissue",
            "score", "pvalue",
        ],
        schema_overrides={"score": pl.Float64, "pvalue": pl.Float64},
        null_values=["NA", ""],
        infer_schema_length=10000,
    )
    return df.with_columns([
        pl.col("target_gene_name").str.strip_chars(),
        pl.col("score").cast(pl.Float32),
        pl.lit("3d_chromatin").alias("evidence_type"),
        pl.lit(None, dtype=pl.Utf8).alias("variant_id"),
    ])


def load_eqtl(path: Path) -> pl.DataFrame:
    """eQTL schema: cCRE | Gene ID | Gene Name | Gene Type | Variant ID | Source | Tissue | Slope | P-value"""
    df = pl.read_csv(
        path, separator="\t", has_header=False,
        new_columns=[
            "ccre_accession_v4", "target_gene_id", "target_gene_name",
            "gene_type", "variant_id", "source", "biosample_or_tissue",
            "score", "pvalue",
        ],
        schema_overrides={"score": pl.Float64, "pvalue": pl.Float64},
        null_values=["NA", ""],
        infer_schema_length=10000,
    )
    return df.with_columns([
        pl.col("target_gene_name").str.strip_chars(),
        pl.col("score").cast(pl.Float32),
        pl.lit("eqtl").alias("evidence_type"),
        pl.lit(None, dtype=pl.Utf8).alias("experiment_id"),
    ])


def main() -> None:
    ctx = stage_start("06_target_evidence", __doc__)
    sc = ctx.stage_cfg
    paths = sc["paths"]

    viz_path = PROJECT_ROOT / paths["precomputed_dir"] / "viz_chr16.parquet"
    v4_bed_path = PROJECT_ROOT / paths["universe_dir"] / "screen_v4_2024-07.bed"
    annot_dir = PROJECT_ROOT / paths["annotations_dir"] / "encode_gene_links"
    chrom3d_path = annot_dir / "3d_chromatin_chr16.tsv"
    eqtl_path = annot_dir / "eqtl_chr16.tsv"
    precomp_dir = PROJECT_ROOT / paths["precomputed_dir"]
    out_path = precomp_dir / "region_target_evidence.parquet"
    out_summary_path = precomp_dir / "region_target_evidence_summary.parquet"

    inputs = [
        file_record(viz_path),
        file_record(v4_bed_path),
        file_record(chrom3d_path),
        file_record(eqtl_path),
    ]

    try:
        for p in (viz_path, v4_bed_path, chrom3d_path, eqtl_path):
            if not p.exists():
                raise FileNotFoundError(f"{p} missing")

        print("loading viz_chr16...", file=sys.stderr)
        viz = pl.read_parquet(viz_path).select(["token_id", "chrom", "start", "end"])
        print(f"  {len(viz):,} chr16 tokens", file=sys.stderr)

        print("loading V4 master cCRE BED (filtering to chr16)...", file=sys.stderr)
        v4 = pl.read_csv(
            v4_bed_path, separator="\t", has_header=False,
            new_columns=["chrom", "start", "end", "dhs_id", "accession", "cclass"],
        ).filter(pl.col("chrom") == "chr16").select(
            ["chrom", "start", "end", "accession"]
        )
        print(f"  {len(v4):,} V4 chr16 cCREs", file=sys.stderr)

        print("\nbuilding viz token → V4 EH38E accession map...", file=sys.stderr)
        token_map = coord_map_viz_to_v4(viz, v4, min_overlap_frac=0.5)
        mapped = token_map.filter(pl.col("ccre_accession_v4").is_not_null())
        n_mapped = len(mapped)
        print(f"  final: {n_mapped:,}/{len(viz):,} tokens have V4 accession",
              file=sys.stderr)

        # Build accession → token_id (one-to-one for our case; if duplicates would
        # arise we'd need to handle multimap, but coord match is unique)
        acc_to_token = {
            r["ccre_accession_v4"]: r["token_id"]
            for r in mapped.iter_rows(named=True)
        }
        print(f"  unique mapped accessions: {len(acc_to_token):,}", file=sys.stderr)

        print("\nloading 3D chromatin chr16 evidence...", file=sys.stderr)
        chr3d = load_3d_chromatin(chrom3d_path)
        print(f"  {len(chr3d):,} rows raw", file=sys.stderr)
        # Filter to mapped accessions
        chr3d_matched = chr3d.filter(
            pl.col("ccre_accession_v4").is_in(list(acc_to_token.keys()))
        )
        print(f"  {len(chr3d_matched):,} rows match mapped tokens "
              f"({len(chr3d_matched) / max(len(chr3d), 1):.1%})", file=sys.stderr)

        print("\nloading eQTL chr16 evidence...", file=sys.stderr)
        eqtl = load_eqtl(eqtl_path)
        print(f"  {len(eqtl):,} rows raw", file=sys.stderr)
        eqtl_matched = eqtl.filter(
            pl.col("ccre_accession_v4").is_in(list(acc_to_token.keys()))
        )
        print(f"  {len(eqtl_matched):,} rows match mapped tokens "
              f"({len(eqtl_matched) / max(len(eqtl), 1):.1%})", file=sys.stderr)

        print("\nbuilding combined evidence parquet...", file=sys.stderr)
        # Add token_id via accession map
        chr3d_matched = chr3d_matched.with_columns(
            pl.col("ccre_accession_v4").replace_strict(acc_to_token, return_dtype=pl.Int64).alias("token_id")
        )
        eqtl_matched = eqtl_matched.with_columns(
            pl.col("ccre_accession_v4").replace_strict(acc_to_token, return_dtype=pl.Int64).alias("token_id")
        )

        # Unify schemas
        cols = [
            "token_id", "evidence_type", "target_gene_id", "target_gene_name",
            "gene_type", "source", "biosample_or_tissue", "experiment_id",
            "variant_id", "score", "pvalue", "ccre_accession_v4",
        ]
        combined = pl.concat([
            chr3d_matched.select(cols),
            eqtl_matched.select(cols),
        ]).sort(["token_id", "evidence_type", "pvalue"])

        precomp_dir.mkdir(parents=True, exist_ok=True)
        combined.write_parquet(out_path)
        print(f"  wrote {out_path} ({out_path.stat().st_size / 1e6:.1f} MB; "
              f"{len(combined):,} rows)", file=sys.stderr)

        # Per-token aggregate summary: counts by evidence type, distinct genes, top genes
        print("\nbuilding per-token summary...", file=sys.stderr)
        summary = (
            combined.group_by("token_id")
            .agg([
                pl.len().alias("n_evidence_rows"),
                pl.col("target_gene_id").n_unique().alias("n_distinct_genes"),
                (pl.col("evidence_type") == "3d_chromatin").sum().alias("n_3d_chromatin"),
                (pl.col("evidence_type") == "eqtl").sum().alias("n_eqtl"),
                pl.col("biosample_or_tissue").n_unique().alias("n_distinct_contexts"),
                # Top-3 gene symbols by frequency (proxy for "most-supported target")
                pl.col("target_gene_name").value_counts(sort=True).head(3).alias("top_genes_struct"),
            ])
            .sort("token_id")
        )
        summary = summary.with_columns(
            pl.col("top_genes_struct").map_elements(
                lambda lst: [item["target_gene_name"] for item in lst],
                return_dtype=pl.List(pl.Utf8),
            ).alias("top_genes")
        ).drop("top_genes_struct")
        summary.write_parquet(out_summary_path)
        print(f"  wrote {out_summary_path} ({out_summary_path.stat().st_size / 1e6:.2f} MB; "
              f"{len(summary):,} tokens with evidence)", file=sys.stderr)

        outputs = [
            file_record(out_path, record_count=len(combined)),
            file_record(out_summary_path, record_count=len(summary)),
        ]
        metrics: dict[str, Any] = {
            "n_viz_tokens": int(len(viz)),
            "n_v4_chr16_ccres": int(len(v4)),
            "n_tokens_mapped_to_v4": int(n_mapped),
            "n_3d_chromatin_rows_raw": int(len(chr3d)),
            "n_3d_chromatin_rows_matched": int(len(chr3d_matched)),
            "n_eqtl_rows_raw": int(len(eqtl)),
            "n_eqtl_rows_matched": int(len(eqtl_matched)),
            "n_evidence_rows_total": int(len(combined)),
            "n_tokens_with_evidence": int(len(summary)),
            "evidence_types_present": ["3d_chromatin", "eqtl"],
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": "Ensure stage 12 has run and gene-links chr16 TSVs exist.",
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
