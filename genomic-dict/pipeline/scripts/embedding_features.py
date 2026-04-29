"""Embedding-derived features for the dictionary entry — prototype distances + concept axes.

For each chr16 token, compute:
  1. Cosine distance to each SCREEN-class centroid (a soft-classification view).
     Output: data/precomputed/region_class_prototypes.parquet
  2. Projection onto a small set of named "concept axes" defined as
     direction vectors in R2V embedding space:
       - activity_axis:    mean(K562 H3K27ac top 25%) − mean(K562 H3K27me3 top 25%)
       - K562_specificity: mean(K562-only H3K27ac high) − mean(non-K562 H3K27ac high)
       - GM12878/HepG2_specificity: analogous
       - anchor_axis:      mean(PLS centroid) − mean(dELS centroid)
     Output: data/precomputed/region_concept_axes.parquet

These ride alongside the corpus-derived stages 12 / 13 — the R2V embedding
already in viz_chr16's geometry plus a bit of derivation, with no separate
stage in the canonical pipeline. Per the plan
(plans/2026-04-28-grammar-of-regions-pipeline.md).

Reads:
  databio/r2v-encode-hg38                      (HuggingFace; uses local cache)
  data/precomputed/viz_chr16.parquet           (chr16 token list + cclass + signal cols)
Writes:
  data/precomputed/region_class_prototypes.parquet
  data/precomputed/region_concept_axes.parquet
  results/embedding_features/summary.json
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from _common import PROJECT_ROOT  # noqa: E402


HF_MODEL = "databio/r2v-encode-hg38"
TOP_FRAC = 0.25  # top-25% threshold for "active" / "specific" reference sets
EMBEDDING_DIM = 100


def cosine_distances(query: np.ndarray, references: np.ndarray) -> np.ndarray:
    """Cosine distance from each query row to each reference row.

    Returns matrix of shape (n_query, n_ref). cos_dist = 1 - cos_sim."""
    q = query / (np.linalg.norm(query, axis=1, keepdims=True) + 1e-12)
    r = references / (np.linalg.norm(references, axis=1, keepdims=True) + 1e-12)
    sims = q @ r.T
    return 1.0 - sims


def project(query: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """Cosine of angle between each query row and the axis direction.

    Bounded [-1, 1]. +1 = aligned with axis; -1 = anti-aligned."""
    a = axis / (np.linalg.norm(axis) + 1e-12)
    q = query / (np.linalg.norm(query, axis=1, keepdims=True) + 1e-12)
    return q @ a


def main() -> None:
    t0 = time.time()
    viz_path = PROJECT_ROOT / "data/precomputed/viz_chr16.parquet"
    precomp_dir = PROJECT_ROOT / "data/precomputed"
    out_proto = precomp_dir / "region_class_prototypes.parquet"
    out_axes = precomp_dir / "region_concept_axes.parquet"
    results_dir = PROJECT_ROOT / "results/embedding_features"
    results_dir.mkdir(parents=True, exist_ok=True)

    print("=== embedding_features ===", file=sys.stderr)
    if not viz_path.exists():
        raise FileNotFoundError(f"{viz_path} missing")

    print(f"loading R2V model ({HF_MODEL}) from cache...", file=sys.stderr)
    from geniml.region2vec.main import Region2VecExModel
    t1 = time.time()
    model = Region2VecExModel(model_path=HF_MODEL)
    print(f"  loaded in {time.time() - t1:.1f}s", file=sys.stderr)

    weight = model.model.projection.weight.detach().cpu().numpy()  # (vocab, dim)
    print(f"  embedding matrix: {weight.shape}", file=sys.stderr)
    if weight.shape[1] != EMBEDDING_DIM:
        print(f"  WARNING: expected dim={EMBEDDING_DIM}, got {weight.shape[1]}", file=sys.stderr)

    print("loading viz_chr16...", file=sys.stderr)
    viz = pl.read_parquet(viz_path).select([
        "token_id", "region", "cclass",
        "K562__ATAC__mean", "K562__H3K27ac__mean", "K562__H3K27me3__mean",
        "K562__H3K4me1__mean", "K562__H3K4me3__mean", "K562__H3K9me3__mean",
        "GM12878__ATAC__mean", "GM12878__H3K27ac__mean", "GM12878__H3K27me3__mean",
        "GM12878__H3K4me1__mean", "GM12878__H3K4me3__mean", "GM12878__H3K9me3__mean",
        "HepG2__ATAC__mean", "HepG2__H3K27ac__mean", "HepG2__H3K27me3__mean",
        "HepG2__H3K4me1__mean", "HepG2__H3K4me3__mean", "HepG2__H3K9me3__mean",
    ])
    n_tokens = len(viz)
    print(f"  {n_tokens:,} chr16 tokens", file=sys.stderr)

    # Slice the weight matrix to chr16 tokens (in viz_chr16 token_id order).
    token_ids = viz["token_id"].to_numpy()
    if (token_ids.max() >= weight.shape[0]) or (token_ids.min() < 0):
        raise ValueError(
            f"token_id out of bounds for embedding matrix "
            f"(min={token_ids.min()}, max={token_ids.max()}, vocab={weight.shape[0]})"
        )
    embed = weight[token_ids]  # (n_tokens, dim)
    print(f"  chr16 embeddings: {embed.shape}", file=sys.stderr)

    # ---------- (1) Class prototypes ----------
    print("\ncomputing class prototypes...", file=sys.stderr)
    cclasses = (
        viz.filter(pl.col("cclass").is_not_null())
        .select("cclass").unique().sort("cclass")["cclass"].to_list()
    )
    print(f"  classes: {cclasses}", file=sys.stderr)

    centroids: dict[str, np.ndarray] = {}
    for cc in cclasses:
        rows_mask = (viz["cclass"] == cc).fill_null(False).to_numpy().astype(bool)
        n_in_class = int(rows_mask.sum())
        if n_in_class == 0:
            continue
        centroid = embed[rows_mask].mean(axis=0)
        centroids[cc] = centroid
        print(f"    {cc:>11}: n={n_in_class:>6,} centroid_norm={np.linalg.norm(centroid):.3f}",
              file=sys.stderr)

    cclass_order = sorted(centroids.keys())
    centroid_matrix = np.stack([centroids[cc] for cc in cclass_order], axis=0)
    distances = cosine_distances(embed, centroid_matrix)  # (n_tokens, n_classes)
    proto_cols = {f"distance_{cc}": distances[:, i].astype(np.float32) for i, cc in enumerate(cclass_order)}

    proto_df = pl.DataFrame({
        "token_id": viz["token_id"],
        **{k: pl.Series(k, v) for k, v in proto_cols.items()},
    })
    proto_df.write_parquet(out_proto)
    print(f"  wrote {out_proto} ({out_proto.stat().st_size / 1e3:.1f} KB)", file=sys.stderr)

    # ---------- (2) Concept axes ----------
    print("\ncomputing concept axes...", file=sys.stderr)

    # Helper: top-fraction mask by a numeric column (returns boolean numpy array)
    def top_mask(col: str, frac: float = TOP_FRAC) -> np.ndarray:
        vals = viz[col].to_numpy()
        if (vals > 0).sum() < int(frac * n_tokens):
            # Column too sparse — relax to "any positive" for the case set.
            return vals > 0
        thresh = np.percentile(vals, 100 * (1 - frac))
        return vals >= thresh

    def mean_emb(mask: np.ndarray) -> np.ndarray:
        if not mask.any():
            return np.zeros(weight.shape[1], dtype=np.float32)
        return embed[mask].mean(axis=0)

    # activity_axis: K562 active mark high vs K562 repressive mark high.
    # Combines H3K27ac+H3K4me1+H3K4me3 for case; H3K27me3 for control.
    k562_active_max = np.maximum.reduce([
        viz["K562__H3K27ac__mean"].to_numpy(),
        viz["K562__H3K4me1__mean"].to_numpy(),
        viz["K562__H3K4me3__mean"].to_numpy(),
    ])
    case_mask = k562_active_max >= np.percentile(k562_active_max, 75)
    ctrl_mask = top_mask("K562__H3K27me3__mean")
    activity_axis = mean_emb(case_mask) - mean_emb(ctrl_mask)
    print(f"  activity_axis: case={int(case_mask.sum()):,}, ctrl={int(ctrl_mask.sum()):,}, "
          f"norm={np.linalg.norm(activity_axis):.3f}", file=sys.stderr)

    # Cell-line specificity axes: high in this cell vs high in others.
    def cell_specificity(this: str, others: list[str]) -> np.ndarray:
        this_max = np.maximum.reduce([
            viz[f"{this}__H3K27ac__mean"].to_numpy(),
            viz[f"{this}__H3K4me1__mean"].to_numpy(),
            viz[f"{this}__H3K4me3__mean"].to_numpy(),
        ])
        others_max = np.maximum.reduce([
            np.maximum.reduce([
                viz[f"{o}__H3K27ac__mean"].to_numpy(),
                viz[f"{o}__H3K4me1__mean"].to_numpy(),
                viz[f"{o}__H3K4me3__mean"].to_numpy(),
            ]) for o in others
        ])
        # Case: top 25% by this_max where others_max is below median (specifically high).
        this_thresh = np.percentile(this_max, 75)
        others_med = np.median(others_max)
        case = (this_max >= this_thresh) & (others_max < others_med)
        # Control: top 25% by others_max where this_max is below median.
        others_thresh = np.percentile(others_max, 75)
        this_med = np.median(this_max)
        ctrl = (others_max >= others_thresh) & (this_max < this_med)
        return mean_emb(case) - mean_emb(ctrl), int(case.sum()), int(ctrl.sum())

    K562_axis, k562_n, k562_c = cell_specificity("K562", ["GM12878", "HepG2"])
    GM12878_axis, gm_n, gm_c = cell_specificity("GM12878", ["K562", "HepG2"])
    HepG2_axis, hep_n, hep_c = cell_specificity("HepG2", ["K562", "GM12878"])
    print(f"  K562_specificity:    case={k562_n:,} ctrl={k562_c:,} norm={np.linalg.norm(K562_axis):.3f}",
          file=sys.stderr)
    print(f"  GM12878_specificity: case={gm_n:,} ctrl={gm_c:,} norm={np.linalg.norm(GM12878_axis):.3f}",
          file=sys.stderr)
    print(f"  HepG2_specificity:   case={hep_n:,} ctrl={hep_c:,} norm={np.linalg.norm(HepG2_axis):.3f}",
          file=sys.stderr)

    # anchor_axis: PLS centroid - dELS centroid. (promoter-like vs distal-enhancer-like)
    if "PLS" in centroids and "dELS" in centroids:
        anchor_axis = centroids["PLS"] - centroids["dELS"]
        print(f"  anchor_axis (PLS-dELS): norm={np.linalg.norm(anchor_axis):.3f}", file=sys.stderr)
    else:
        anchor_axis = np.zeros(weight.shape[1], dtype=np.float32)
        print(f"  anchor_axis: missing PLS or dELS centroid; zeros", file=sys.stderr)

    # Project all chr16 tokens onto each axis.
    activity_score = project(embed, activity_axis).astype(np.float32)
    K562_score = project(embed, K562_axis).astype(np.float32)
    GM12878_score = project(embed, GM12878_axis).astype(np.float32)
    HepG2_score = project(embed, HepG2_axis).astype(np.float32)
    anchor_score = project(embed, anchor_axis).astype(np.float32)

    axes_df = pl.DataFrame({
        "token_id": viz["token_id"],
        "activity_score": activity_score,
        "K562_specificity_score": K562_score,
        "GM12878_specificity_score": GM12878_score,
        "HepG2_specificity_score": HepG2_score,
        "anchor_score": anchor_score,
    })
    axes_df.write_parquet(out_axes)
    print(f"\n  wrote {out_axes} ({out_axes.stat().st_size / 1e3:.1f} KB)", file=sys.stderr)

    # Sanity check: classes should split sensibly along anchor_axis.
    print("\nsanity: per-class median anchor_score (should run PLS > pELS > CA-H3K4me3 > CA-CTCF > dELS):",
          file=sys.stderr)
    joined = viz.select(["token_id", "cclass"]).join(axes_df, on="token_id")
    for cc in cclass_order:
        med = joined.filter(pl.col("cclass") == cc)["anchor_score"].median()
        n = (joined["cclass"] == cc).sum()
        print(f"    {cc:>11}: n={int(n):>6,} median anchor_score = {med:>7.3f}", file=sys.stderr)

    summary = {
        "n_chr16_tokens": int(n_tokens),
        "embedding_dim": int(weight.shape[1]),
        "n_cclasses": len(cclass_order),
        "cclasses": cclass_order,
        "concept_axes": ["activity", "K562_specificity", "GM12878_specificity", "HepG2_specificity", "anchor"],
        "duration_seconds": round(time.time() - t0, 2),
        "outputs": {
            "region_class_prototypes": str(out_proto.relative_to(PROJECT_ROOT)),
            "region_concept_axes": str(out_axes.relative_to(PROJECT_ROOT)),
        },
    }
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {results_dir / 'summary.json'}", file=sys.stderr)
    print(f"total time: {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
