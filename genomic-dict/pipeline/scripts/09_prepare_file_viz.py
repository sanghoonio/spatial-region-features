"""Stage 09 — prepare file-level viz data (file UMAP).

Joins the corpus manifest (from stage 01) with per-file mean embeddings (from
stage 11), runs UMAP locally with umap-learn, writes viz-ready parquet.

The input embeddings are GENOME-WIDE means (one 100-dim vector per file,
averaged over all non-UNK tokens regardless of chromosome — see stage 11).
So this UMAP shows files clustered by their full-genome regulatory profile,
not by their chr16 profile. The chr16 restriction in mattress's narrative
applies to the region-level UMAP (stage 08), not this file-level one.

Output schema:
  - id, name, description
  - umap_x, umap_y          (computed here)
  - assay, cell_line, cell_type, tissue
  - is_unlabeled            (always False — the new corpus drops UNKNOWN cell_line)

Replaces the prior approach of joining bedbase-umap's precomputed coords
(which were biased toward a different curated subset).

Reads:
  data/corpus/manifest.parquet                (from stage 01)
  data/precomputed/file_embeddings.parquet    (from stage 11)
Writes:
  data/precomputed/viz_files.parquet          — viz-ready
  results/09_prepare_file_viz/summary.json
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def main() -> None:
    ctx = stage_start("09_prepare_file_viz", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    manifest_path = PROJECT_ROOT / paths["corpus_dir"] / "manifest.parquet"
    embeddings_relative = sc.get("embeddings_relative",
                                 "data/precomputed/file_embeddings.parquet")
    embeddings_path = PROJECT_ROOT / embeddings_relative
    precomp_dir = PROJECT_ROOT / paths["precomputed_dir"]
    out_path = precomp_dir / "viz_files.parquet"

    n_neighbors = int(sc.get("umap_n_neighbors", 15))
    min_dist = float(sc.get("umap_min_dist", 0.1))
    n_components = int(sc.get("umap_n_components", 2))
    seed = int(sc.get("umap_seed", 42))

    inputs = [file_record(manifest_path), file_record(embeddings_path)]

    try:
        for p in [manifest_path, embeddings_path]:
            if not p.exists():
                raise FileNotFoundError(f"{p} missing — check upstream stages.")

        manifest = pl.read_parquet(manifest_path)
        emb_df = pl.read_parquet(embeddings_path)
        print(f"  manifest: {len(manifest):,}, embeddings: {len(emb_df):,}",
              file=sys.stderr)

        joined = manifest.join(emb_df.select(["id", "embedding"]), on="id", how="inner")
        n_joined = len(joined)
        print(f"    joined: {n_joined:,}", file=sys.stderr)
        if n_joined == 0:
            raise RuntimeError("No id overlap between manifest and embeddings.")

        # Stack embeddings into (N, D) numpy array. polars list-of-float column
        # → list of lists → np.asarray.
        emb_matrix = np.asarray(joined["embedding"].to_list(), dtype=np.float32)
        print(f"    embedding matrix: {emb_matrix.shape}", file=sys.stderr)

        import umap
        t0 = time.time()
        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            random_state=seed,
            metric="cosine",
        )
        coords = reducer.fit_transform(emb_matrix)
        print(f"    umap fit: {time.time() - t0:.1f}s", file=sys.stderr)

        out = (joined.with_columns([
                    pl.Series("umap_x", coords[:, 0].astype(np.float32)),
                    pl.Series("umap_y", coords[:, 1].astype(np.float32)),
                    pl.lit(False).alias("is_unlabeled"),
                ])
               .select([
                    "id", "name",
                    pl.col("description").fill_null("").alias("description"),
                    "umap_x", "umap_y",
                    "assay", "cell_line",
                    pl.col("cell_type").fill_null("").alias("cell_type"),
                    pl.col("tissue").fill_null("").alias("tissue"),
                    "is_unlabeled",
                ])
               .sort("id"))

        precomp_dir.mkdir(parents=True, exist_ok=True)
        out.write_parquet(out_path)
        outputs = [file_record(out_path, record_count=len(out))]

        # Distributional metrics for the summary.
        assay_counts = (out.group_by("assay").agg(pl.len().alias("n"))
                        .sort("n", descending=True).to_dicts())
        cell_line_counts = (out.group_by("cell_line").agg(pl.len().alias("n"))
                            .sort("n", descending=True).head(15).to_dicts())

        metrics = {
            "n_manifest": len(manifest),
            "n_embeddings": len(emb_df),
            "n_joined": n_joined,
            "n_final": len(out),
            "umap_n_neighbors": n_neighbors,
            "umap_min_dist": min_dist,
            "umap_seed": seed,
            "umap_metric": "cosine",
            "assay_counts": assay_counts,
            "top_cell_lines": cell_line_counts,
            "file_size_bytes": out_path.stat().st_size,
            "embedding_dim": int(emb_matrix.shape[1]),
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(f"\n  wrote {out_path.name}: {len(out):,} rows "
              f"({out_path.stat().st_size / 1e6:.1f} MB)", file=sys.stderr)

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": "If file_embeddings.parquet missing, run stage 11 on Rivanna first.",
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
