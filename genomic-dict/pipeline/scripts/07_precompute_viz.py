"""Stage 08 — precompute region visualization artifacts for the focus chromosome.

This stage is the final region-side aggregator. It absorbs what used to be
stage 04 (intrinsic bigwig features) and the standalone embedding_features.py,
and now derives concept axes from BED peak-call counts rather than bigwig
means — keeping the embedding-side computation grounded in the same data
substrate the R2V model was trained on.

Reads:
  data/annotations/pretrained_universe.parquet      (from stage 07; has 100-dim embeddings)
  data/precomputed/tokenized_corpus_chr16.parquet   (from stage 11; per-file active token lists)
  data/corpus/manifest.parquet                       (from stage 01; file metadata)

Writes:
  data/precomputed/viz_<chrom>.parquet              — token + class + UMAP + kNN
  data/precomputed/region_concept_axes.parquet     — per-token projections onto 5 axes
  results/08_precompute_viz/summary.json            — AI-ingestible summary

Stage 04 is gone; the old `intrinsic_bigwig.parquet` + per-(cell, mark) bigwig
mean columns are no longer part of the viz substrate. Bigwigs are sampled
only by stage 14 for the featured intervals' continuous-mode plot.
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


TOP_FRAC = 0.25
ACTIVE_MARKS = ["H3K27ac", "H3K4me1", "H3K4me3"]
POLYCOMB_MARK = "H3K27me3"
SPEC_CELLS = ["K562", "GM12878", "HepG2"]


def n_files_active_per_token(
    tokenized: pl.DataFrame,
    manifest: pl.DataFrame,
    *,
    cell_line: str | None,
    targets: list[str] | None,
    assay_substring: str | None,
    chr16_token_to_idx: dict[int, int],
) -> np.ndarray:
    """For each chr16 token, count files matching (cell_line, targets, assay)
    that have the token peak-called. Returns an int64 array indexed by the
    chr16 universe order (parallel to chr16_token_to_idx values).

    Any None argument relaxes that filter dimension.
    """
    expr = pl.lit(True)
    if cell_line is not None:
        expr = expr & (pl.col("cell_line") == cell_line)
    if targets is not None:
        expr = expr & pl.col("target").is_in(targets)
    if assay_substring is not None:
        expr = expr & pl.col("assay").str.contains(assay_substring)
    file_ids = manifest.filter(expr)["id"].to_list()

    counts = np.zeros(len(chr16_token_to_idx), dtype=np.int64)
    if not file_ids:
        return counts

    sub = (
        tokenized.filter(pl.col("id").is_in(file_ids))
        .explode("chr16_active_token_ids")
        .rename({"chr16_active_token_ids": "token_id"})
        .group_by("token_id")
        .agg(pl.len().alias("n"))
    )
    for row in sub.iter_rows(named=True):
        idx = chr16_token_to_idx.get(int(row["token_id"]))
        if idx is not None:
            counts[idx] = int(row["n"])
    return counts


def project(query: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Cosine of angle between each query row and an axis direction."""
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        return np.zeros(query.shape[0], dtype=np.float32)
    a = axis / norm
    q = query / (np.linalg.norm(query, axis=1, keepdims=True) + 1e-12)
    return q @ a


def mean_emb(embed: np.ndarray, mask: np.ndarray) -> np.ndarray:
    if not mask.any():
        return np.zeros(embed.shape[1], dtype=np.float32)
    return embed[mask].mean(axis=0)


def main() -> None:
    ctx = stage_start("07_precompute_viz", __doc__)
    sc = ctx.stage_cfg

    chrom = sc["focus_chromosome"]
    annot_dir = PROJECT_ROOT / sc["paths"]["annotations_dir"]
    precomp_dir = PROJECT_ROOT / sc["paths"]["precomputed_dir"]
    corpus_dir = PROJECT_ROOT / sc["paths"]["corpus_dir"]

    pretrained_path = annot_dir / "pretrained_universe.parquet"
    tokenized_path = precomp_dir / "tokenized_corpus_chr16.parquet"
    manifest_path = corpus_dir / "manifest.parquet"

    out_path = precomp_dir / f"viz_{chrom}.parquet"
    out_axes_path = precomp_dir / "region_concept_axes.parquet"

    umap_n_neighbors = int(sc.get("umap_n_neighbors", 15))
    umap_min_dist = float(sc.get("umap_min_dist", 0.1))
    umap_seed = int(sc.get("umap_seed", 42))
    knn_k = int(sc.get("knn_k", 30))

    inputs = [
        file_record(pretrained_path),
        file_record(tokenized_path),
        file_record(manifest_path),
    ]

    try:
        for p in (pretrained_path, tokenized_path, manifest_path):
            if not p.exists():
                raise FileNotFoundError(f"{p} missing — check upstream stage.")

        # ---- Load pretrained universe (filter to focus chrom) ----
        print(f"  loading pretrained universe...", file=sys.stderr)
        universe = (
            pl.read_parquet(pretrained_path)
            .filter(pl.col("chrom") == chrom)
            .sort(["start", "end"])
        )
        n = len(universe)
        print(f"    {n:,} {chrom} regions", file=sys.stderr)

        emb_matrix = np.array(universe["embedding"].to_list(), dtype=np.float32)
        print(f"    embedding matrix: {emb_matrix.shape}", file=sys.stderr)

        chr16_token_ids = universe["token_id"].to_list()
        chr16_token_to_idx = {int(t): i for i, t in enumerate(chr16_token_ids)}

        # ---- UMAP ----
        print(
            f"  computing UMAP (seed={umap_seed}, n_neighbors={umap_n_neighbors}, "
            f"min_dist={umap_min_dist})...",
            file=sys.stderr,
        )
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

        # ---- k-NN (exact, cosine) ----
        print(f"  computing k-NN (k={knn_k})...", file=sys.stderr)
        t0 = time.time()
        from sklearn.neighbors import NearestNeighbors
        nn = NearestNeighbors(
            n_neighbors=knn_k + 1, metric="cosine", algorithm="brute"
        )
        nn.fit(emb_matrix)
        distances, indices = nn.kneighbors(emb_matrix)
        knn_indices = indices[:, 1:]      # drop self
        knn_distances = distances[:, 1:]
        knn_seconds = round(time.time() - t0, 1)
        print(f"    knn {knn_seconds}s; graph shape {knn_indices.shape}", file=sys.stderr)
        token_id_arr = np.asarray(chr16_token_ids, dtype=np.int64)
        knn_token_ids = token_id_arr[knn_indices]

        # ---- BED-derived case/ctrl masks for concept axes ----
        # (replaces the bigwig-signal-based masks the old embedding_features.py used)
        print("  loading tokenized corpus + manifest for BED-derived axes...", file=sys.stderr)
        tokenized = pl.read_parquet(tokenized_path).select(["id", "chr16_active_token_ids"])
        manifest = pl.read_parquet(manifest_path).select(["id", "cell_line", "assay", "target"])
        # Inner join so we only count files with both manifest metadata and tokenized output.
        manifest = manifest.join(
            tokenized.select(["id"]), on="id", how="inner",
        )

        # Per-token activity: peak-called in K562 active mark files vs K562 H3K27me3 files.
        n_active = n_files_active_per_token(
            tokenized, manifest,
            cell_line="K562", targets=ACTIVE_MARKS, assay_substring="ChIP",
            chr16_token_to_idx=chr16_token_to_idx,
        )
        n_polycomb = n_files_active_per_token(
            tokenized, manifest,
            cell_line="K562", targets=[POLYCOMB_MARK], assay_substring="ChIP",
            chr16_token_to_idx=chr16_token_to_idx,
        )

        def top_mask(arr: np.ndarray, frac: float = TOP_FRAC) -> np.ndarray:
            if (arr > 0).sum() < int(frac * len(arr)):
                return arr > 0
            thresh = np.percentile(arr, 100 * (1 - frac))
            return arr >= thresh

        case_act = top_mask(n_active)
        ctrl_act = top_mask(n_polycomb)
        activity_axis = mean_emb(emb_matrix, case_act) - mean_emb(emb_matrix, ctrl_act)
        print(
            f"    activity axis: case={int(case_act.sum()):,}, "
            f"ctrl={int(ctrl_act.sum()):,}, norm={np.linalg.norm(activity_axis):.3f}",
            file=sys.stderr,
        )

        # Cell-line specificity: case = high in this_cell active marks AND low elsewhere;
        # ctrl = high in other_cells active marks AND low in this_cell.
        per_cell_active: dict[str, np.ndarray] = {}
        for cell in SPEC_CELLS:
            per_cell_active[cell] = n_files_active_per_token(
                tokenized, manifest,
                cell_line=cell, targets=ACTIVE_MARKS, assay_substring="ChIP",
                chr16_token_to_idx=chr16_token_to_idx,
            )

        spec_axes: dict[str, np.ndarray] = {}
        for cell in SPEC_CELLS:
            this_n = per_cell_active[cell]
            others_n = np.maximum.reduce(
                [per_cell_active[c] for c in SPEC_CELLS if c != cell]
            )
            this_thresh = np.percentile(this_n, 75)
            others_med = np.median(others_n)
            case = (this_n >= this_thresh) & (others_n < others_med)
            others_thresh = np.percentile(others_n, 75)
            this_med = np.median(this_n)
            ctrl = (others_n >= others_thresh) & (this_n < this_med)
            spec_axes[cell] = mean_emb(emb_matrix, case) - mean_emb(emb_matrix, ctrl)
            print(
                f"    {cell} specificity: case={int(case.sum()):,}, "
                f"ctrl={int(ctrl.sum()):,}, norm={np.linalg.norm(spec_axes[cell]):.3f}",
                file=sys.stderr,
            )

        # Anchor axis: PLS centroid - dELS centroid (no BED needed; class-driven).
        cclasses = (
            universe.filter(pl.col("cclass").is_not_null())
            ["cclass"].unique().sort().to_list()
        )
        centroids: dict[str, np.ndarray] = {}
        for cc in cclasses:
            mask = (universe["cclass"] == cc).fill_null(False).to_numpy().astype(bool)
            if mask.any():
                centroids[cc] = emb_matrix[mask].mean(axis=0)

        if "PLS" in centroids and "dELS" in centroids:
            anchor_axis = centroids["PLS"] - centroids["dELS"]
        else:
            anchor_axis = np.zeros(emb_matrix.shape[1], dtype=np.float32)
        print(
            f"    anchor axis (PLS-dELS): norm={np.linalg.norm(anchor_axis):.3f}",
            file=sys.stderr,
        )

        # Project all chr16 tokens onto each axis.
        anchor_score = project(emb_matrix, anchor_axis).astype(np.float32)
        activity_score = project(emb_matrix, activity_axis).astype(np.float32)
        K562_score = project(emb_matrix, spec_axes["K562"]).astype(np.float32)
        GM12878_score = project(emb_matrix, spec_axes["GM12878"]).astype(np.float32)
        HepG2_score = project(emb_matrix, spec_axes["HepG2"]).astype(np.float32)

        # ---- Build outputs ----
        # viz_chr16: identity + class + UMAP + kNN. Drops embedding (kept in
        # pretrained_universe). No bigwig means.
        out = (
            universe.drop("embedding").with_columns([
                pl.Series("umap_x", coords[:, 0].astype(np.float32)),
                pl.Series("umap_y", coords[:, 1].astype(np.float32)),
                pl.Series("knn_token_ids", [list(map(int, row)) for row in knn_token_ids]),
                pl.Series("knn_distances", [[float(d) for d in row] for row in knn_distances]),
            ])
        )

        axes_df = pl.DataFrame({
            "token_id": universe["token_id"],
            "anchor_score": anchor_score,
            "activity_score": activity_score,
            "K562_specificity_score": K562_score,
            "GM12878_specificity_score": GM12878_score,
            "HepG2_specificity_score": HepG2_score,
        })

        precomp_dir.mkdir(parents=True, exist_ok=True)
        out.write_parquet(out_path)
        axes_df.write_parquet(out_axes_path)

        outputs = [
            file_record(out_path, record_count=len(out)),
            file_record(out_axes_path, record_count=len(axes_df)),
        ]

        class_counts = (
            out.group_by("cclass").agg(pl.len().alias("n"))
            .sort("n", descending=True).to_dicts()
        )

        # Sanity: anchor_score median per class (should run PLS > pELS > … > dELS).
        anchor_medians: list[dict[str, Any]] = []
        joined = universe.select(["token_id", "cclass"]).join(axes_df, on="token_id")
        for cc in cclasses:
            med = joined.filter(pl.col("cclass") == cc)["anchor_score"].median()
            anchor_medians.append({"cclass": cc, "anchor_score_median": float(med or 0.0)})

        metrics = {
            "focus_chromosome": chrom,
            "n_regions": len(out),
            "n_columns_viz": len(out.columns),
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
            "concept_axes": [
                "anchor", "activity",
                "K562_specificity", "GM12878_specificity", "HepG2_specificity",
            ],
            "anchor_score_medians_by_class": anchor_medians,
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(
            f"\n  wrote {out_path.name} ({out_path.stat().st_size / 1e6:.1f} MB) "
            f"+ {out_axes_path.name} ({out_axes_path.stat().st_size / 1e3:.1f} KB)",
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
                    "Prerequisites: stage 07 (pretrained universe), stage 11 "
                    "(tokenized corpus), stage 01 (manifest). Stage 11 must run "
                    "before stage 08 since concept axes are now BED-derived."
                ),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
