"""Pre-flight: R2V kNN vs corpus-wide PPMI agreement.

For each chr16 token, compute:
  - top-30 partners by R2V kNN (from viz_chr16.knn_token_ids)
  - top-30 partners by corpus-wide PPMI (computed here on tokenized_corpus_chr16)
  - Jaccard between the two sets

Report distribution of per-token Jaccards. Tells us whether R2V faithfully
captured our 16,799-file corpus's distributional structure, or whether
per-stratum PMI is essential rather than value-add.

Per the plan (plans/2026-04-28-grammar-of-regions-pipeline.md):
  - High agreement (mean Jaccard > 0.5) → R2V is doing most of the work;
    per-stratum PMI is value-add narrative richness.
  - Low agreement (mean Jaccard < 0.3) → R2V and our curated corpus diverge
    significantly; per-stratum PMI is essential, not optional.

Memory note: chr16 cooccurrence is dense — ~36k nonzero pairs per token →
~1.3B total. We never materialize the full X.T @ X. Instead we slab-process:
~CHUNK columns at a time, compute PPMI + top-K inline, discard the slab.

Inputs:
  data/annotations/tokenized_corpus_chr16.parquet   per-file chr16 token lists
  data/precomputed/viz_chr16.parquet                knn_token_ids per chr16 token
Writes:
  results/preflight_r2v_faithfulness/summary.json   distribution stats + verdict
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).parent))
from _common import PROJECT_ROOT  # noqa: E402


TOP_K = 30
MIN_FILES_ACTIVE = 5
CHUNK = 200  # tokens per slab; balances memory vs matmul overhead


def main() -> None:
    t0 = time.time()
    tokenized_path = PROJECT_ROOT / "data/annotations/tokenized_corpus_chr16.parquet"
    viz_chr16_path = PROJECT_ROOT / "data/precomputed/viz_chr16.parquet"
    out_dir = PROJECT_ROOT / "results/preflight_r2v_faithfulness"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "summary.json"

    print("=== preflight_r2v_faithfulness ===", file=sys.stderr)

    print("loading tokenized corpus...", file=sys.stderr)
    tokenized = pl.read_parquet(tokenized_path)
    n_files = len(tokenized)
    print(f"  {n_files:,} files", file=sys.stderr)

    print("loading viz_chr16 (for kNN partners)...", file=sys.stderr)
    viz = pl.read_parquet(viz_chr16_path).select(["token_id", "knn_token_ids"])
    print(f"  {len(viz):,} chr16 tokens", file=sys.stderr)

    chr16_tokens = sorted(viz["token_id"].to_list())
    token_to_idx = {t: i for i, t in enumerate(chr16_tokens)}
    n_tokens = len(chr16_tokens)
    chr16_tokens_arr = np.asarray(chr16_tokens, dtype=np.int64)

    print("building sparse (file × token) activation matrix...", file=sys.stderr)
    rows: list[int] = []
    cols: list[int] = []
    for i, tok_list in enumerate(tokenized["chr16_active_token_ids"]):
        for t in tok_list:
            j = token_to_idx.get(t)
            if j is not None:
                rows.append(i)
                cols.append(j)
    data = np.ones(len(rows), dtype=np.float32)
    X = csr_matrix(
        (data, (np.asarray(rows, dtype=np.int32), np.asarray(cols, dtype=np.int32))),
        shape=(n_files, n_tokens),
        dtype=np.float32,
    )
    del rows, cols, data
    print(f"  X shape={X.shape}, nnz={X.nnz:,}", file=sys.stderr)

    n_per_token = np.asarray(X.sum(axis=0)).flatten()
    print(
        f"  per-token marginals: min={int(n_per_token.min())}, "
        f"median={int(np.median(n_per_token))}, "
        f"max={int(n_per_token.max())}",
        file=sys.stderr,
    )
    keep_mask = n_per_token >= MIN_FILES_ACTIVE
    n_keep = int(keep_mask.sum())
    print(
        f"  tokens passing floor (>= {MIN_FILES_ACTIVE} files): {n_keep:,}/{n_tokens:,}",
        file=sys.stderr,
    )

    print(f"\nchunked PPMI top-{TOP_K} extraction (CHUNK={CHUNK})...", file=sys.stderr)
    X_csr = X.tocsr()  # (files × tokens)
    X_csc = X.tocsc()  # for column-slab extraction
    log_N = np.log(float(n_files))

    ppmi_top: dict[int, set[int]] = {}
    n_chunks = (n_tokens + CHUNK - 1) // CHUNK
    t_start = time.time()

    for ci in range(n_chunks):
        c0 = ci * CHUNK
        c1 = min(c0 + CHUNK, n_tokens)
        # Slab of cooc counts for tokens [c0, c1): shape (n_tokens, c1-c0)
        slab = X_csr.T @ X_csc[:, c0:c1]  # sparse
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
            # PPMI = max(0, log(c * N / (n_a * n_b)))
            n_a = n_per_token[tok_idx]
            n_b = n_per_token[partners_idx]
            ppmi = np.log(joint) + log_N - np.log(n_a) - np.log(n_b)
            ppmi = np.maximum(ppmi, 0.0)
            pos = ppmi > 0
            partners_idx = partners_idx[pos]
            ppmi = ppmi[pos]
            if partners_idx.size == 0:
                continue
            if partners_idx.size > TOP_K:
                top_local = np.argpartition(-ppmi, TOP_K - 1)[:TOP_K]
                partners_idx = partners_idx[top_local]
            ppmi_top[int(chr16_tokens_arr[tok_idx])] = {
                int(chr16_tokens_arr[p]) for p in partners_idx
            }

        if (ci + 1) % 20 == 0 or ci + 1 == n_chunks:
            elapsed = time.time() - t_start
            rate = (ci + 1) / elapsed
            eta = (n_chunks - ci - 1) / rate if rate > 0 else 0
            print(
                f"  chunk {ci+1:4d}/{n_chunks} "
                f"({(ci+1)/n_chunks*100:5.1f}%) "
                f"elapsed={elapsed:6.1f}s eta={eta:6.1f}s "
                f"ppmi_top={len(ppmi_top):,}",
                file=sys.stderr,
            )

    print(
        f"  PPMI extraction done; tokens with at least one partner: {len(ppmi_top):,}",
        file=sys.stderr,
    )

    print(f"\nloading R2V kNN top-{TOP_K} partners...", file=sys.stderr)
    knn_top: dict[int, set[int]] = {}
    for row in viz.iter_rows(named=True):
        tok = int(row["token_id"])
        knn_ids = row["knn_token_ids"]
        if knn_ids is None:
            continue
        knn_top[tok] = {int(t) for t in list(knn_ids)[:TOP_K]}

    common = set(ppmi_top.keys()) & set(knn_top.keys())
    print(f"  tokens with both PPMI and kNN sets: {len(common):,}", file=sys.stderr)

    jaccards = []
    for t in common:
        s_p = ppmi_top[t]
        s_k = knn_top[t]
        if not s_p or not s_k:
            continue
        inter = len(s_p & s_k)
        union = len(s_p | s_k)
        jaccards.append(inter / union if union > 0 else 0.0)
    jacc = np.asarray(jaccards, dtype=np.float64)
    print(f"  Jaccards computed for {len(jacc):,} tokens", file=sys.stderr)

    print("\nJaccard distribution (R2V kNN top-30 vs corpus-wide PPMI top-30):", file=sys.stderr)
    print(f"  mean   = {jacc.mean():.3f}", file=sys.stderr)
    print(f"  median = {float(np.median(jacc)):.3f}", file=sys.stderr)
    print(f"  q25    = {float(np.percentile(jacc, 25)):.3f}", file=sys.stderr)
    print(f"  q75    = {float(np.percentile(jacc, 75)):.3f}", file=sys.stderr)
    print(f"  > 0.5  = {(jacc > 0.5).mean():.3f}", file=sys.stderr)
    print(f"  < 0.3  = {(jacc < 0.3).mean():.3f}", file=sys.stderr)
    print(f"  == 0   = {(jacc == 0).mean():.3f}", file=sys.stderr)

    bins = np.arange(0, 1.01, 0.1)
    hist, _ = np.histogram(jacc, bins=bins)
    print("\n  histogram:", file=sys.stderr)
    h_max = max(hist) if hist.size else 1
    for i in range(len(bins) - 1):
        bar = "#" * int(hist[i] / h_max * 50)
        print(f"  {bins[i]:.1f}-{bins[i+1]:.1f}  {hist[i]:6d}  {bar}", file=sys.stderr)

    if jacc.mean() > 0.5:
        verdict = (
            "R2V is doing most of the work; per-stratum PMI is value-add "
            "narrative richness."
        )
    elif jacc.mean() < 0.3:
        verdict = (
            "R2V and our corpus diverge significantly; per-stratum PMI "
            "is essential, not optional."
        )
    else:
        verdict = (
            "Intermediate result; per-stratum PMI is a meaningful complement "
            "to R2V kNN, neither redundant nor essential."
        )
    print(f"\nverdict: {verdict}", file=sys.stderr)

    summary = {
        "n_files": int(n_files),
        "n_chr16_tokens": int(n_tokens),
        "n_tokens_passing_floor": int(n_keep),
        "min_files_active": MIN_FILES_ACTIVE,
        "top_k": TOP_K,
        "chunk_size": CHUNK,
        "n_tokens_with_ppmi_partners": len(ppmi_top),
        "n_tokens_with_knn_partners": len(knn_top),
        "n_tokens_compared": int(len(jacc)),
        "jaccard_mean": float(jacc.mean()),
        "jaccard_median": float(np.median(jacc)),
        "jaccard_q25": float(np.percentile(jacc, 25)),
        "jaccard_q75": float(np.percentile(jacc, 75)),
        "fraction_above_0_5": float((jacc > 0.5).mean()),
        "fraction_below_0_3": float((jacc < 0.3).mean()),
        "fraction_zero": float((jacc == 0).mean()),
        "histogram_bins": bins.tolist(),
        "histogram_counts": [int(h) for h in hist],
        "verdict": verdict,
        "duration_seconds": round(time.time() - t0, 2),
    }
    out_path.write_text(json.dumps(summary, indent=2))
    print(f"\nwrote {out_path}", file=sys.stderr)
    print(f"total time: {time.time() - t0:.1f}s", file=sys.stderr)


if __name__ == "__main__":
    main()
