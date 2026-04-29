"""Stage 08 — precompute visualization artifacts for the focus chromosome.

Reads the pretrained universe (with SCREEN class tagging) and the stage-04
bigwig features, filters to the focus chromosome, computes:

  * 2D UMAP on the 100-dim pretrained embeddings (deterministic seed).
  * k-NN graph on the same embeddings (cosine distance, k from config).

Writes one consolidated parquet that Observable can load directly for all
three panels (entry card, neighborhood network, diagnostic).

Reads:
  data/annotations/pretrained_universe.parquet      (from stage 07)
  data/annotations/intrinsic_bigwig.parquet         (from stage 04)
Writes:
  data/precomputed/viz_<chrom>.parquet              — one row per region, joined
  results/08_precompute_viz/summary.json            — AI-ingestible summary
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
    ctx = stage_start("08_precompute_viz", __doc__)
    sc = ctx.stage_cfg

    chrom = sc["focus_chromosome"]
    annot_dir = PROJECT_ROOT / sc["paths"]["annotations_dir"]
    precomp_dir = PROJECT_ROOT / sc["paths"]["precomputed_dir"]
    pretrained_path = annot_dir / "pretrained_universe.parquet"
    intrinsic_path = annot_dir / "intrinsic_bigwig.parquet"
    out_path = precomp_dir / f"viz_{chrom}.parquet"

    umap_n_neighbors = int(sc.get("umap_n_neighbors", 15))
    umap_min_dist = float(sc.get("umap_min_dist", 0.1))
    umap_seed = int(sc.get("umap_seed", 42))
    knn_k = int(sc.get("knn_k", 30))

    inputs = [file_record(pretrained_path), file_record(intrinsic_path)]

    try:
        if not pretrained_path.exists():
            raise FileNotFoundError(f"{pretrained_path} missing (run stage 07).")
        if not intrinsic_path.exists():
            raise FileNotFoundError(f"{intrinsic_path} missing (run stage 04).")

        print(f"  loading pretrained universe...", file=sys.stderr)
        universe = pl.read_parquet(pretrained_path).filter(pl.col("chrom") == chrom)
        n = len(universe)
        print(f"    {n:,} {chrom} regions", file=sys.stderr)

        print(f"  loading stage-04 bigwig features...", file=sys.stderr)
        intrinsic = pl.read_parquet(intrinsic_path)

        # Join on token_id (both tables have it; intrinsic kept it via stage 04 load_focus_universe).
        joined = universe.join(intrinsic, on="token_id", how="inner", suffix="_intr")

        # Drop duplicate coordinate columns that come from the right side.
        drop_cols = [c for c in joined.columns if c.endswith("_intr") and c.replace("_intr", "") in joined.columns]
        if drop_cols:
            joined = joined.drop(drop_cols)
        print(f"    joined shape: {len(joined):,} × {len(joined.columns)} cols", file=sys.stderr)

        # Pull embeddings as a numpy array for UMAP + kNN.
        emb_matrix = np.array(joined["embedding"].to_list(), dtype=np.float32)
        print(f"    embedding matrix: {emb_matrix.shape}", file=sys.stderr)

        # UMAP.
        print(f"  computing UMAP (seed={umap_seed}, n_neighbors={umap_n_neighbors}, min_dist={umap_min_dist})...", file=sys.stderr)
        t0 = time.time()
        import umap
        reducer = umap.UMAP(
            n_neighbors=umap_n_neighbors,
            min_dist=umap_min_dist,
            n_components=2,
            random_state=umap_seed,
            metric="cosine",
        )
        coords = reducer.fit_transform(emb_matrix)
        umap_seconds = round(time.time() - t0, 1)
        print(f"    umap {umap_seconds}s; coords shape {coords.shape}", file=sys.stderr)

        # k-NN (exact, cosine) on embeddings.
        print(f"  computing k-NN (k={knn_k})...", file=sys.stderr)
        t0 = time.time()
        from sklearn.neighbors import NearestNeighbors
        nn = NearestNeighbors(n_neighbors=knn_k + 1, metric="cosine", algorithm="brute")
        nn.fit(emb_matrix)
        distances, indices = nn.kneighbors(emb_matrix)
        # Drop self (first column).
        knn_indices = indices[:, 1:]       # (n, k)
        knn_distances = distances[:, 1:]   # (n, k)
        knn_seconds = round(time.time() - t0, 1)
        print(f"    knn {knn_seconds}s; graph shape {knn_indices.shape}", file=sys.stderr)

        # Map kNN local row indices → token_id for cross-table lookups in the viz.
        token_ids = joined["token_id"].to_numpy()
        knn_token_ids = token_ids[knn_indices]  # (n, k)

        # Build output dataframe. Drop the 100-dim embedding to keep Observable payload small;
        # embedding is preserved in pretrained_universe.parquet if needed.
        out = joined.drop("embedding").with_columns([
            pl.Series("umap_x", coords[:, 0].astype(np.float32)),
            pl.Series("umap_y", coords[:, 1].astype(np.float32)),
            pl.Series("knn_token_ids", [list(map(int, row)) for row in knn_token_ids]),
            pl.Series("knn_distances", [[float(d) for d in row] for row in knn_distances]),
        ])

        precomp_dir.mkdir(parents=True, exist_ok=True)
        out.write_parquet(out_path)
        outputs = [file_record(out_path, record_count=len(out))]

        # Per-class summary on the focus-chrom subset.
        class_counts = (
            out.group_by("cclass").agg(pl.len().alias("n"))
            .sort("n", descending=True).to_dicts()
        )

        metrics = {
            "focus_chromosome": chrom,
            "n_regions": len(out),
            "n_columns": len(out.columns),
            "umap_seconds": umap_seconds,
            "umap_params": {
                "n_neighbors": umap_n_neighbors,
                "min_dist": umap_min_dist,
                "seed": umap_seed,
                "metric": "cosine",
            },
            "knn_seconds": knn_seconds,
            "knn_k": knn_k,
            "per_class_counts": class_counts,
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(
            f"\n  wrote {out_path.name}: {len(out):,} rows × {len(out.columns)} cols "
            f"({out_path.stat().st_size / 1e6:.1f} MB)",
            file=sys.stderr,
        )

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": "Prerequisites: stages 07 (pretrained universe) and 04 (bigwig features). Rerun either if missing.",
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
