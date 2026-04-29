"""Stage 13 — Leiden community detection on NPMI cooccurrence graphs.

For each (stratum, gamma) pair, build an undirected weighted graph from
stage 12's region_cooccurrence_pmi.parquet (edges = top-K NPMI partners
above threshold; weights = NPMI), run Leiden community detection, and
produce per-token module assignments + per-module summaries.

Per the plan (plans/2026-04-28-grammar-of-regions-pipeline.md, log entry
2026-04-28 step 2):
  - NPMI is the primary edge metric (PPMI ties are pathological)
  - gamma = 1.0 hardcoded for v1 (Locked-in decision #3)
  - Edge threshold default 0.3 (permissive; AG cluster edges are NPMI 0.4-0.99)

Reads:
  data/precomputed/region_cooccurrence_pmi.parquet  (from stage 12)
  data/precomputed/viz_chr16.parquet                (cclass for module labels)
Writes:
  data/precomputed/region_modules.parquet           (per-token × stratum × gamma)
  data/precomputed/module_summary.parquet           (per-module catalogue)
  results/13_modules/summary.json
"""
from __future__ import annotations

import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import igraph as ig
import leidenalg as la
import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def build_graph_for_stratum(
    cooc_stratum: pl.DataFrame,
    npmi_threshold: float,
) -> tuple[ig.Graph, list[int]]:
    """Build an undirected weighted graph from a stratum slice of the cooc parquet.

    Edges are symmetrized: if both (a, b) and (b, a) appear in their respective
    top-K lists with NPMI above threshold, take the max. Returns the graph and
    the parallel list of token_ids (vertex ids index into this list).
    """
    edges: dict[tuple[int, int], float] = {}
    for r in cooc_stratum.iter_rows(named=True):
        a = r["token_id"]
        partners = list(r["partner_token_ids"] or [])
        weights = list(r["weights_npmi"] or [])
        for b, w in zip(partners, weights):
            if w < npmi_threshold:
                continue
            key = (a, b) if a < b else (b, a)
            prev = edges.get(key)
            if prev is None or w > prev:
                edges[key] = float(w)

    if not edges:
        # empty graph
        return ig.Graph(directed=False), []

    nodes = sorted({n for pair in edges for n in pair})
    node_to_idx = {n: i for i, n in enumerate(nodes)}
    edge_list = [(node_to_idx[a], node_to_idx[b]) for (a, b) in edges]
    weight_list = [edges[(a, b)] for (a, b) in edges]

    g = ig.Graph(n=len(nodes), edges=edge_list, directed=False)
    g.es["weight"] = weight_list
    return g, nodes


def detect_communities(
    g: ig.Graph,
    gamma: float,
    seed: int,
) -> list[int]:
    """Leiden community detection. Returns membership list parallel to g.vs."""
    if g.vcount() == 0:
        return []
    partition = la.find_partition(
        g,
        la.RBConfigurationVertexPartition,
        weights="weight",
        resolution_parameter=gamma,
        seed=seed,
    )
    return list(partition.membership)


def eigenvector_centrality(
    g: ig.Graph,
    member_idx: list[int],
) -> np.ndarray:
    """Per-vertex eigenvector centrality on the within-module subgraph."""
    if not member_idx:
        return np.array([])
    sub = g.subgraph(member_idx)
    if sub.vcount() <= 1:
        return np.ones(sub.vcount(), dtype=np.float64)
    try:
        return np.asarray(sub.eigenvector_centrality(weights="weight", scale=False))
    except Exception:
        # Fall back to degree centrality if eigenvector fails to converge
        return np.asarray(sub.strength(weights="weight"), dtype=np.float64)


def main() -> None:
    ctx = stage_start("13_modules", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    cooc_path = PROJECT_ROOT / paths["precomputed_dir"] / "region_cooccurrence_pmi.parquet"
    viz_path = PROJECT_ROOT / paths["precomputed_dir"] / "viz_chr16.parquet"
    precomp_dir = PROJECT_ROOT / paths["precomputed_dir"]
    out_modules = precomp_dir / "region_modules.parquet"
    out_summary = precomp_dir / "module_summary.parquet"

    inputs = [file_record(cooc_path), file_record(viz_path)]

    try:
        if not cooc_path.exists():
            raise FileNotFoundError(f"{cooc_path} missing — run stage 12 first.")
        if not viz_path.exists():
            raise FileNotFoundError(f"{viz_path} missing.")

        npmi_thresh = float(sc["npmi_threshold"])
        resolutions = list(sc["resolutions"])
        seed = int(sc["random_seed"])

        print(f"loading cooc parquet...", file=sys.stderr)
        cooc = pl.read_parquet(cooc_path)
        print(f"  {len(cooc):,} rows; columns: {cooc.columns}", file=sys.stderr)

        print(f"loading viz_chr16 (for cclass)...", file=sys.stderr)
        viz = pl.read_parquet(viz_path).select(["token_id", "region", "cclass"])
        token_meta = {
            r["token_id"]: (r["region"], r["cclass"]) for r in viz.iter_rows(named=True)
        }

        strata = cooc["stratum"].unique().sort().to_list()
        print(f"strata to process: {len(strata)}  resolutions: {resolutions}\n", file=sys.stderr)

        modules_rows: list[dict[str, Any]] = []
        summary_rows: list[dict[str, Any]] = []
        per_stratum_metrics: dict[str, dict[str, Any]] = {}

        for stratum in strata:
            sub = cooc.filter(pl.col("stratum") == stratum)
            print(f"--- stratum: {stratum}  rows={len(sub):,} ---", file=sys.stderr)

            t0 = time.time()
            g, nodes = build_graph_for_stratum(sub, npmi_threshold=npmi_thresh)
            t_build = time.time() - t0

            if g.vcount() == 0:
                print(f"  empty graph (no edges above NPMI {npmi_thresh}) — skipping", file=sys.stderr)
                per_stratum_metrics[stratum] = {"n_nodes": 0, "n_edges": 0, "n_modules_per_gamma": {}}
                continue

            print(f"  graph: {g.vcount():,} nodes, {g.ecount():,} edges  (built in {t_build:.1f}s)",
                  file=sys.stderr)

            n_modules_per_gamma: dict[float, int] = {}
            for gamma in resolutions:
                t1 = time.time()
                membership = detect_communities(g, gamma, seed=seed)
                t_leiden = time.time() - t1

                # group token indices by module
                groups: dict[int, list[int]] = defaultdict(list)
                for v_idx, mid in enumerate(membership):
                    groups[mid].append(v_idx)

                n_mods = len(groups)
                n_modules_per_gamma[gamma] = n_mods
                largest = max(len(g_) for g_ in groups.values())
                singletons = sum(1 for g_ in groups.values() if len(g_) == 1)
                print(f"  γ={gamma}: {n_mods:,} modules  (largest={largest}, "
                      f"singletons={singletons})  Leiden took {t_leiden:.1f}s",
                      file=sys.stderr)

                # Per-module: centrality, anchor, summary
                for mid, member_idx in groups.items():
                    cent = eigenvector_centrality(g, member_idx)
                    if cent.size == 0:
                        continue
                    anchor_local = int(np.argmax(cent))
                    anchor_token = nodes[member_idx[anchor_local]]

                    member_tokens = [nodes[i] for i in member_idx]
                    member_classes = [
                        (token_meta.get(t, ("?", None))[1] or "unclassed") for t in member_tokens
                    ]
                    class_counts = Counter(member_classes)
                    dominant_class, _ = class_counts.most_common(1)[0]

                    anchor_region, anchor_cls = token_meta.get(anchor_token, ("?", None))
                    auto_label = (
                        f"{stratum}/γ{gamma}/m{mid}: "
                        f"{dominant_class}-dominant ({len(member_tokens)} regions, "
                        f"anchor {anchor_token} {anchor_region})"
                    )

                    for local_i, vert_i in enumerate(member_idx):
                        modules_rows.append({
                            "token_id": int(nodes[vert_i]),
                            "stratum": stratum,
                            "gamma": float(gamma),
                            "module_id": int(mid),
                            "within_module_centrality": float(cent[local_i]),
                            "is_anchor": bool(local_i == anchor_local),
                        })
                    summary_rows.append({
                        "stratum": stratum,
                        "gamma": float(gamma),
                        "module_id": int(mid),
                        "n_tokens": int(len(member_tokens)),
                        "anchor_token_id": int(anchor_token),
                        "anchor_region": anchor_region,
                        "anchor_cclass": anchor_cls or "unclassed",
                        "dominant_class": dominant_class,
                        "class_counts": dict(class_counts),
                        "member_token_ids": [int(t) for t in member_tokens],
                        "auto_label": auto_label,
                    })

            per_stratum_metrics[stratum] = {
                "n_nodes": int(g.vcount()),
                "n_edges": int(g.ecount()),
                "n_modules_per_gamma": {str(k): v for k, v in n_modules_per_gamma.items()},
            }

        if not modules_rows:
            raise RuntimeError("no modules emitted — every stratum was empty after threshold")

        print(f"\nwriting {len(modules_rows):,} module-membership rows...", file=sys.stderr)
        modules_df = pl.DataFrame(
            modules_rows,
            schema={
                "token_id": pl.Int64,
                "stratum": pl.Utf8,
                "gamma": pl.Float32,
                "module_id": pl.Int64,
                "within_module_centrality": pl.Float32,
                "is_anchor": pl.Boolean,
            },
        )
        modules_df.write_parquet(out_modules)
        print(f"  wrote {out_modules} ({out_modules.stat().st_size / 1e6:.1f} MB)", file=sys.stderr)

        print(f"writing {len(summary_rows):,} module-summary rows...", file=sys.stderr)
        summary_df = pl.DataFrame(
            summary_rows,
            schema={
                "stratum": pl.Utf8,
                "gamma": pl.Float32,
                "module_id": pl.Int64,
                "n_tokens": pl.Int64,
                "anchor_token_id": pl.Int64,
                "anchor_region": pl.Utf8,
                "anchor_cclass": pl.Utf8,
                "dominant_class": pl.Utf8,
                "class_counts": pl.Object,
                "member_token_ids": pl.List(pl.Int64),
                "auto_label": pl.Utf8,
            },
        )
        # class_counts is a dict; polars stores as Object — drop on write or stringify
        summary_df = summary_df.with_columns(
            pl.col("class_counts").map_elements(
                lambda d: ", ".join(f"{k}={v}" for k, v in sorted(d.items())),
                return_dtype=pl.Utf8,
            )
        )
        summary_df.write_parquet(out_summary)
        print(f"  wrote {out_summary} ({out_summary.stat().st_size / 1e6:.2f} MB)", file=sys.stderr)

        outputs = [
            file_record(out_modules, record_count=len(modules_df)),
            file_record(out_summary, record_count=len(summary_df)),
        ]
        metrics: dict[str, Any] = {
            "n_strata_processed": len([s for s, m in per_stratum_metrics.items() if m["n_nodes"] > 0]),
            "n_strata_empty": len([s for s, m in per_stratum_metrics.items() if m["n_nodes"] == 0]),
            "npmi_threshold": npmi_thresh,
            "resolutions": resolutions,
            "n_module_rows": len(modules_df),
            "n_summary_rows": len(summary_df),
            "per_stratum": per_stratum_metrics,
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": "Run stage 12 (region_cooccurrence_pmi.parquet) first.",
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
