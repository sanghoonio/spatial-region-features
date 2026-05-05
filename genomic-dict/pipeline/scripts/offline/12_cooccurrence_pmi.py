"""Stage 12 — curated-strata PMI cooccurrence on chr16 tokens.

For each named biological-question stratum (defined as a manifest filter rule
in config.yaml's stages.12_cooccurrence_pmi.strata), compute the top-K PPMI
partners per chr16 token. PPMI discounts popularity-bias; the curated strata
keep each cooccurrence cohort biology-coherent (each stratum asks one named
question — see plans/2026-04-28-grammar-of-regions-pipeline.md).

Algorithm per stratum:
  1. Filter manifest by `filters` block (logical AND across keys, OR within
     a key; missing key = no filter; empty filters = full corpus).
  2. Slice the (file × token) sparse activation matrix to those rows.
  3. Recompute per-token marginals; apply min_files_active floor.
  4. Chunked sparse matmul: for tokens [c0, c1), slab = X.T @ X[:, c0:c1].
  5. Per column in the slab: PPMI = max(0, log(c · N / (n_a · n_b))); also
     compute Jaccard = c / (n_a + n_b - c). Take top-K by PPMI.
  6. Append rows to the long-format output parquet.

Reads:
  data/corpus/manifest.parquet                    (from stage 01; cell_line, assay, target)
  data/precomputed/tokenized_corpus_chr16.parquet (from stage 11; chr16 token lists)
Writes:
  data/precomputed/region_cooccurrence_pmi.parquet
    schema: (token_id, stratum, n_files_in_stratum, n_files_active,
             partner_token_ids, weights_ppmi, weights_jaccard, counts)
  results/12_cooccurrence_pmi/summary.json
"""
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def _and_block_mask(
    manifest: pl.DataFrame,
    filters: dict[str, list[str]],
) -> np.ndarray:
    """AND across keys, OR within a key. Empty = all True. Nulls excluded."""
    if not filters:
        return np.ones(len(manifest), dtype=bool)
    mask = np.ones(len(manifest), dtype=bool)
    for col, values in filters.items():
        if col not in manifest.columns:
            raise KeyError(f"stratum filter references column '{col}' missing from manifest")
        mask &= manifest[col].is_in(values).fill_null(False).to_numpy().astype(bool)
    return mask


def stratum_mask(
    manifest: pl.DataFrame,
    stratum_def: dict[str, Any],
) -> np.ndarray:
    """Boolean mask of files matching this stratum's filter rule.

    Two schemas supported:
      - `filters: {col: [values], ...}` — single AND-block (legacy / simple).
      - `filter_blocks: [{...}, {...}, ...]` — list of AND-blocks, OR'd.
        Used for case-control / disjunctive contrast strata.

    Both empty/missing means "no filter" (full corpus)."""
    blocks = stratum_def.get("filter_blocks")
    if blocks is not None:
        if not isinstance(blocks, list):
            raise ValueError("filter_blocks must be a list of filter dicts")
        mask = np.zeros(len(manifest), dtype=bool)
        for blk in blocks:
            mask |= _and_block_mask(manifest, blk or {})
        return mask
    return _and_block_mask(manifest, stratum_def.get("filters") or {})


def build_activation_matrix(
    tokenized: pl.DataFrame,
    chr16_token_ids: list[int],
) -> tuple[csr_matrix, list[int]]:
    """Build sparse (file × token) activation matrix and return ordered file ids."""
    token_to_idx = {t: i for i, t in enumerate(chr16_token_ids)}
    rows: list[int] = []
    cols: list[int] = []
    file_ids: list[int] = []
    for i, (fid, tok_list) in enumerate(zip(
        tokenized["id"].to_list(),
        tokenized["chr16_active_token_ids"].to_list(),
    )):
        file_ids.append(fid)
        if tok_list is None:
            continue
        for t in tok_list:
            j = token_to_idx.get(t)
            if j is not None:
                rows.append(i)
                cols.append(j)
    data = np.ones(len(rows), dtype=np.float32)
    X = csr_matrix(
        (data, (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32))),
        shape=(len(tokenized), len(chr16_token_ids)),
        dtype=np.float32,
    )
    return X, file_ids


def compute_stratum_top_k(
    X_sub: csr_matrix,
    chr16_tokens_arr: np.ndarray,
    *,
    top_k: int,
    ppmi_threshold: float,
    min_files_active: int,
    chunk_size: int,
    stratum_name: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Compute top-K PPMI partners per token within this stratum.

    Returns (rows for output parquet, per-stratum metrics)."""
    n_files = X_sub.shape[0]
    n_tokens = X_sub.shape[1]
    if n_files == 0:
        print(f"  [{stratum_name}] EMPTY stratum — skipping", file=sys.stderr)
        return [], {"n_files": 0, "n_tokens_emitted": 0}

    n_per_token = np.asarray(X_sub.sum(axis=0)).flatten()
    keep_mask = n_per_token >= min_files_active
    n_keep = int(keep_mask.sum())
    print(
        f"  [{stratum_name}] n_files={n_files:,}  n_tokens_active={int((n_per_token > 0).sum()):,}  "
        f"passing floor (>={min_files_active}): {n_keep:,}",
        file=sys.stderr,
    )

    if n_keep == 0:
        return [], {"n_files": n_files, "n_tokens_emitted": 0}

    log_N = np.log(float(n_files))
    X_csr = X_sub.tocsr()
    X_csc = X_sub.tocsc()

    rows_out: list[dict[str, Any]] = []
    n_chunks = (n_tokens + chunk_size - 1) // chunk_size
    t_start = time.time()

    for ci in range(n_chunks):
        c0 = ci * chunk_size
        c1 = min(c0 + chunk_size, n_tokens)
        slab = X_csr.T @ X_csc[:, c0:c1]
        slab_csc = slab.tocsc()
        for j in range(c1 - c0):
            tok_idx = c0 + j
            if not keep_mask[tok_idx]:
                continue
            col = slab_csc.getcol(j).tocoo()
            partners_idx = col.row
            joint = col.data
            self_mask = partners_idx != tok_idx
            partners_idx = partners_idx[self_mask]
            joint = joint[self_mask]
            floor = keep_mask[partners_idx]
            partners_idx = partners_idx[floor]
            joint = joint[floor]
            if partners_idx.size == 0:
                continue

            n_a = float(n_per_token[tok_idx])
            n_b = n_per_token[partners_idx].astype(np.float64)
            joint_f = joint.astype(np.float64)
            ppmi = np.log(joint_f) + log_N - np.log(n_a) - np.log(n_b)
            ppmi = np.maximum(ppmi, 0.0)
            keep = ppmi > ppmi_threshold
            partners_idx = partners_idx[keep]
            joint_f = joint_f[keep]
            ppmi = ppmi[keep]
            n_b = n_b[keep]
            if partners_idx.size == 0:
                continue

            # NPMI = PPMI / -log(P(a,b)). Normalizes by joint probability so
            # high-mass perfectly-correlated partners outrank low-mass ones.
            # Without this, raw PPMI ties at log(N/n_a) for all fully-contained
            # partners (often hundreds), and argpartition picks arbitrarily.
            log_p_ab = np.log(joint_f) - log_N  # negative
            npmi = ppmi / -log_p_ab

            jaccard = joint_f / (n_a + n_b - joint_f)
            # Rank by NPMI desc, tie-break by joint count desc (more empirical
            # evidence wins ties). lexsort uses last key as primary sort.
            order = np.lexsort((-joint_f, -npmi))
            if partners_idx.size > top_k:
                order = order[:top_k]
            partners_idx = partners_idx[order]
            ppmi = ppmi[order]
            npmi = npmi[order]
            jaccard = jaccard[order]
            joint_f = joint_f[order]

            rows_out.append({
                "token_id": int(chr16_tokens_arr[tok_idx]),
                "stratum": stratum_name,
                "n_files_in_stratum": int(n_files),
                "n_files_active": int(n_per_token[tok_idx]),
                "partner_token_ids": [int(chr16_tokens_arr[p]) for p in partners_idx],
                "weights_npmi": [float(n) for n in npmi],
                "weights_ppmi": [float(p) for p in ppmi],
                "weights_jaccard": [float(j) for j in jaccard],
                "counts": [int(c) for c in joint_f],
            })

        if (ci + 1) % 40 == 0 or ci + 1 == n_chunks:
            elapsed = time.time() - t_start
            rate = (ci + 1) / elapsed if elapsed > 0 else 0
            eta = (n_chunks - ci - 1) / rate if rate > 0 else 0
            print(
                f"    chunk {ci+1:4d}/{n_chunks} "
                f"({(ci+1)/n_chunks*100:5.1f}%)  "
                f"elapsed={elapsed:6.1f}s  eta={eta:6.1f}s  "
                f"rows={len(rows_out):,}",
                file=sys.stderr,
            )

    # Per-token marginals — one row per chr16 token regardless of partner
    # ranking. Powers stratum-population brushing in the viz (canvas2):
    # define "active in stratum" by frac_active threshold, then subset the
    # UMAP to that population without needing a per-anchor query.
    marginal_rows: list[dict[str, Any]] = []
    for tok_idx in range(n_tokens):
        n_active = int(n_per_token[tok_idx])
        if n_active == 0:
            continue
        marginal_rows.append({
            "token_id": int(chr16_tokens_arr[tok_idx]),
            "stratum": stratum_name,
            "n_files_in_stratum": int(n_files),
            "n_files_active": n_active,
            "frac_active": float(n_active) / float(n_files),
        })

    metrics = {
        "n_files": int(n_files),
        "n_tokens_passing_floor": int(n_keep),
        "n_tokens_emitted": len(rows_out),
        "n_marginal_rows": len(marginal_rows),
        "elapsed_seconds": round(time.time() - t_start, 2),
    }
    return rows_out, marginal_rows, metrics


def main() -> None:
    ctx = stage_start("12_cooccurrence_pmi", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    manifest_path = PROJECT_ROOT / paths["corpus_dir"] / "manifest.parquet"
    tokenized_path = PROJECT_ROOT / paths["precomputed_dir"] / "tokenized_corpus_chr16.parquet"
    viz_chr16_path = PROJECT_ROOT / paths["precomputed_dir"] / "viz_chr16.parquet"
    precomp_dir = PROJECT_ROOT / paths["precomputed_dir"]
    out_path = precomp_dir / "region_cooccurrence_pmi.parquet"
    marginals_path = precomp_dir / "region_stratum_marginals.parquet"

    inputs = [
        file_record(manifest_path),
        file_record(tokenized_path),
        file_record(viz_chr16_path),
    ]

    try:
        for p in (manifest_path, tokenized_path, viz_chr16_path):
            if not p.exists():
                raise FileNotFoundError(f"{p} missing.")

        top_k = int(sc["top_k_partners"])
        ppmi_threshold = float(sc["ppmi_threshold"])
        min_files_active = int(sc["min_files_active"])
        chunk_size = int(sc["chunk_size"])
        strata_cfg = sc["strata"]
        if not isinstance(strata_cfg, dict) or not strata_cfg:
            raise ValueError("config.yaml: stages.12_cooccurrence_pmi.strata must be a non-empty mapping")

        print("loading manifest...", file=sys.stderr)
        manifest_full = pl.read_parquet(manifest_path)
        print(f"  {len(manifest_full):,} rows; cols={manifest_full.columns}", file=sys.stderr)

        print("loading viz_chr16 (token universe)...", file=sys.stderr)
        viz = pl.read_parquet(viz_chr16_path).select(["token_id"])
        chr16_tokens = sorted(viz["token_id"].to_list())
        chr16_tokens_arr = np.asarray(chr16_tokens, dtype=np.int64)
        print(f"  {len(chr16_tokens):,} chr16 tokens", file=sys.stderr)

        print("loading tokenized corpus...", file=sys.stderr)
        tokenized = pl.read_parquet(tokenized_path)
        print(f"  {len(tokenized):,} files in tokenized parquet", file=sys.stderr)

        # Inner-join manifest onto tokenized so file rows align with manifest rows
        # and we only keep files present in both. Keeps `id` ordering deterministic
        # (sort by id after join).
        joined = (
            tokenized.select(["id", "chr16_active_token_ids"])
            .join(
                manifest_full.select(["id", "cell_line", "assay", "target"]),
                on="id",
                how="inner",
            )
            .sort("id")
        )
        n_joined = len(joined)
        print(f"  manifest ∩ tokenized: {n_joined:,} files", file=sys.stderr)

        manifest_aligned = joined.select(["id", "cell_line", "assay", "target"])
        tokenized_aligned = joined.select(["id", "chr16_active_token_ids"])

        print("\nbuilding sparse (file × token) activation matrix...", file=sys.stderr)
        X, file_ids = build_activation_matrix(tokenized_aligned, chr16_tokens)
        print(f"  X shape={X.shape}, nnz={X.nnz:,}", file=sys.stderr)

        all_rows: list[dict[str, Any]] = []
        all_marginals: list[dict[str, Any]] = []
        per_stratum_metrics: dict[str, dict[str, Any]] = {}

        for stratum_name, stratum_def in strata_cfg.items():
            stratum_def = stratum_def or {}
            if "filter_blocks" in stratum_def:
                desc = f"filter_blocks×{len(stratum_def['filter_blocks'])}"
            else:
                desc = f"filters={dict(stratum_def.get('filters') or {})}"
            print(f"\n--- stratum: {stratum_name}  {desc} ---", file=sys.stderr)
            mask = stratum_mask(manifest_aligned, stratum_def)
            row_idx = np.where(mask)[0]
            if row_idx.size == 0:
                print(f"  [{stratum_name}] no files match — skipping", file=sys.stderr)
                per_stratum_metrics[stratum_name] = {"n_files": 0, "n_tokens_emitted": 0}
                continue
            X_sub = X[row_idx, :]
            rows, marginal_rows, metrics = compute_stratum_top_k(
                X_sub, chr16_tokens_arr,
                top_k=top_k,
                ppmi_threshold=ppmi_threshold,
                min_files_active=min_files_active,
                chunk_size=chunk_size,
                stratum_name=stratum_name,
            )
            all_rows.extend(rows)
            all_marginals.extend(marginal_rows)
            per_stratum_metrics[stratum_name] = metrics

        print(f"\nwriting output: {len(all_rows):,} rows total", file=sys.stderr)
        if not all_rows:
            raise RuntimeError("no rows emitted — every stratum either empty or all tokens fail floor")

        out_df = pl.DataFrame(
            all_rows,
            schema={
                "token_id": pl.Int64,
                "stratum": pl.Utf8,
                "n_files_in_stratum": pl.Int64,
                "n_files_active": pl.Int64,
                "partner_token_ids": pl.List(pl.Int64),
                "weights_npmi": pl.List(pl.Float32),
                "weights_ppmi": pl.List(pl.Float32),
                "weights_jaccard": pl.List(pl.Float32),
                "counts": pl.List(pl.Int64),
            },
        )
        precomp_dir.mkdir(parents=True, exist_ok=True)
        out_df.write_parquet(out_path)
        print(f"  wrote {out_path} ({out_path.stat().st_size / 1e6:.1f} MB)", file=sys.stderr)

        # Per-token marginals: one row per (token, stratum) where the token
        # is active in at least one file. Used by canvas2 to subset the UMAP
        # to a stratum-active population for brushing.
        marginals_df = pl.DataFrame(
            all_marginals,
            schema={
                "token_id": pl.Int64,
                "stratum": pl.Utf8,
                "n_files_in_stratum": pl.Int64,
                "n_files_active": pl.Int64,
                "frac_active": pl.Float64,
            },
        )
        marginals_df.write_parquet(marginals_path)
        print(
            f"  wrote {marginals_path} ({marginals_path.stat().st_size / 1e6:.1f} MB, "
            f"{len(marginals_df):,} rows)",
            file=sys.stderr,
        )

        outputs = [
            file_record(out_path, record_count=len(out_df)),
            file_record(marginals_path, record_count=len(marginals_df)),
        ]

        metrics: dict[str, Any] = {
            "n_files_joined": int(n_joined),
            "n_chr16_tokens": int(len(chr16_tokens)),
            "n_strata": len(strata_cfg),
            "top_k_partners": top_k,
            "min_files_active": min_files_active,
            "ppmi_threshold": ppmi_threshold,
            "chunk_size": chunk_size,
            "n_rows_total": len(all_rows),
            "per_stratum": per_stratum_metrics,
            "output_size_bytes": int(out_path.stat().st_size),
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": "Run stages 01 (manifest), 08 (viz_chr16), 11 (tokenize) first.",
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
