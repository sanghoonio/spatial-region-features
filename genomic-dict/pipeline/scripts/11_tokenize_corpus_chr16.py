"""Stage 11 — tokenize the curated corpus, emit chr16 tokens AND file embeddings.

Single-pass tokenization against the pretrained R2V universe (~1.06M tokens,
all chromosomes). For each BED file we emit two artifacts with DIFFERENT scope:

  - 100-dim mean embedding — GENOME-WIDE.
      Mean of per-token embeddings over every non-UNK token in the file,
      regardless of chromosome. Looked up from pretrained_universe.parquet.
      Drives stage 09's file UMAP.
      → data/precomputed/file_embeddings.parquet

  - chr16-active token list — CHR16-ONLY.
      Tokens filtered to focus_chromosome (chr16). Drives the chr16-region
      viz and downstream cooccurrence / module stages.
      → data/annotations/tokenized_corpus_chr16.parquet

R2V is Word2Vec-derived (no doc vector), so the genome-wide mean is
mathematically equivalent to model.encode(rs, pooling="mean").mean(axis=0).

Streaming writers, batch_size to bound memory.

Reads:
  data/corpus/manifest.parquet                       (from stage 01)
  data/annotations/pretrained_universe.parquet       (from stage 07; token_id + embedding + chrom)
  databio/r2v-encode-hg38                            (HuggingFace; for tokenizer)
  BBCLIENT_CACHE env                                 (Rivanna only)
Writes:
  data/annotations/tokenized_corpus_chr16.parquet    — id, chr16_active_token_ids, n_chr16_active
  data/precomputed/file_embeddings.parquet           — id, embedding (list[float], 100-dim)
  results/11_tokenize_corpus_chr16/summary.json
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, file_record, stage_start, write_summary,
)


def main() -> None:
    ctx = stage_start("11_tokenize_corpus_chr16", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    manifest_path = PROJECT_ROOT / paths["corpus_dir"] / "manifest.parquet"
    pretrained_path = PROJECT_ROOT / paths["annotations_dir"] / "pretrained_universe.parquet"
    annot_dir = PROJECT_ROOT / paths["annotations_dir"]
    precomp_dir = PROJECT_ROOT / paths["precomputed_dir"]
    out_chr16_path = annot_dir / "tokenized_corpus_chr16.parquet"
    out_emb_path = precomp_dir / "file_embeddings.parquet"

    hf_model = sc["hf_model"]
    focus_chrom = sc["focus_chromosome"]
    batch_size = int(sc.get("batch_size", 5000))
    drop_unk_only = bool(sc.get("drop_unk_only", True))

    bbclient_cache = os.environ.get("BBCLIENT_CACHE")
    inputs = [file_record(manifest_path), file_record(pretrained_path)]

    try:
        if not bbclient_cache:
            raise RuntimeError("BBCLIENT_CACHE env var not set. source the lab env.sh.")
        for p in [manifest_path, pretrained_path]:
            if not p.exists():
                raise FileNotFoundError(f"{p} missing.")

        from geniml.bbclient import BBClient
        from geniml.region2vec.main import Region2VecExModel

        print(f"  loading pretrained tokenizer from {hf_model}...", file=sys.stderr)
        model = Region2VecExModel(model_path=hf_model)
        tokenizer = model.tokenizer
        unk_id = tokenizer.unk_token_id

        # Pretrained universe — token embeddings + chrom for chr16 filter.
        # Build a numpy lookup table indexed by token_id (bound by max id) so
        # the per-file inner loop is just slicing + mean.
        print(f"  loading pretrained universe + embeddings...", file=sys.stderr)
        universe = pl.read_parquet(pretrained_path)
        if "embedding" not in universe.columns:
            raise RuntimeError(
                "pretrained_universe.parquet missing 'embedding' column — "
                "stage 07 needs to be rerun with embedding extraction."
            )
        token_ids = universe["token_id"].to_numpy()
        emb_list = universe["embedding"].to_list()
        emb_dim = len(emb_list[0])
        max_token_id = int(token_ids.max())
        emb_table = np.zeros((max_token_id + 1, emb_dim), dtype=np.float32)
        for tid, vec in zip(token_ids.tolist(), emb_list):
            emb_table[int(tid)] = np.asarray(vec, dtype=np.float32)
        chr16_set = set(
            universe.filter(pl.col("chrom") == focus_chrom)["token_id"].to_list()
        )
        print(f"    {len(token_ids):,} tokens, dim={emb_dim}, "
              f"chr16={len(chr16_set):,}", file=sys.stderr)

        bbc = BBClient(cache_folder=bbclient_cache)
        manifest = pl.read_parquet(manifest_path).select("id")
        ids = manifest["id"].to_list()
        n = len(ids)
        print(f"  tokenizing {n:,} files...", file=sys.stderr)

        annot_dir.mkdir(parents=True, exist_ok=True)
        precomp_dir.mkdir(parents=True, exist_ok=True)
        chr16_schema = pa.schema([
            pa.field("id", pa.string()),
            pa.field("chr16_active_token_ids", pa.list_(pa.int64())),
            pa.field("n_chr16_active", pa.int64()),
        ])
        emb_schema = pa.schema([
            pa.field("id", pa.string()),
            pa.field("n_tokens", pa.int64()),
            pa.field("embedding", pa.list_(pa.float32())),
        ])
        chr16_writer = pq.ParquetWriter(str(out_chr16_path), chr16_schema, compression="zstd")
        emb_writer = pq.ParquetWriter(str(out_emb_path), emb_schema, compression="zstd")

        chr16_buf_ids: list[str] = []
        chr16_buf_tokens: list[list[int]] = []
        chr16_buf_counts: list[int] = []
        emb_buf_ids: list[str] = []
        emb_buf_n: list[int] = []
        emb_buf_vec: list[list[float]] = []

        n_kept_chr16 = 0
        n_kept_emb = 0
        n_unk_only = 0
        n_no_tokens = 0
        n_failed = 0
        unk_only_samples: list[str] = []
        failed_samples: list[tuple[str, str]] = []

        def flush():
            if chr16_buf_ids:
                chr16_writer.write_table(pa.table({
                    "id": chr16_buf_ids,
                    "chr16_active_token_ids": chr16_buf_tokens,
                    "n_chr16_active": chr16_buf_counts,
                }, schema=chr16_schema))
                chr16_buf_ids.clear()
                chr16_buf_tokens.clear()
                chr16_buf_counts.clear()
            if emb_buf_ids:
                emb_writer.write_table(pa.table({
                    "id": emb_buf_ids,
                    "n_tokens": emb_buf_n,
                    "embedding": emb_buf_vec,
                }, schema=emb_schema))
                emb_buf_ids.clear()
                emb_buf_n.clear()
                emb_buf_vec.clear()

        report_every = max(1, n // 40)
        t0 = time.time()
        for i, bed_id in enumerate(ids):
            try:
                rs = bbc.load_bed(bed_id)
                toks = list(tokenizer.encode(tokenizer.tokenize(rs)))
                non_unk = [t for t in toks if t != unk_id]

                is_unk_only = len(toks) > 0 and len(non_unk) == 0
                if is_unk_only:
                    n_unk_only += 1
                    if len(unk_only_samples) < 20:
                        unk_only_samples.append(bed_id)
                    if drop_unk_only:
                        continue
                if not non_unk:
                    n_no_tokens += 1
                    continue

                # Per-file mean embedding over all non-UNK tokens (genome-wide).
                vec = emb_table[np.asarray(non_unk, dtype=np.int64)].mean(axis=0)
                emb_buf_ids.append(bed_id)
                emb_buf_n.append(len(non_unk))
                emb_buf_vec.append(vec.astype(np.float32).tolist())
                n_kept_emb += 1

                # chr16-active tokens (deduped, sorted) for stages 12/13.
                chr16_toks = sorted(set(t for t in non_unk if t in chr16_set))
                if chr16_toks:
                    chr16_buf_ids.append(bed_id)
                    chr16_buf_tokens.append(chr16_toks)
                    chr16_buf_counts.append(len(chr16_toks))
                    n_kept_chr16 += 1
            except Exception as e:
                n_failed += 1
                if len(failed_samples) < 10:
                    failed_samples.append((bed_id, f"{type(e).__name__}: {e}"))

            if len(chr16_buf_ids) >= batch_size or len(emb_buf_ids) >= batch_size:
                flush()

            if (i + 1) % report_every == 0 or i == n - 1:
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed if elapsed > 0 else 0
                eta = (n - i - 1) / rate if rate > 0 else 0
                print(
                    f"    {i+1:,}/{n:,} ({(i+1)/n:.1%}) "
                    f"rate={rate:.1f} files/s eta={eta/60:.1f}m "
                    f"emb={n_kept_emb} chr16={n_kept_chr16} "
                    f"unk={n_unk_only} failed={n_failed}",
                    file=sys.stderr,
                )

        flush()
        chr16_writer.close()
        emb_writer.close()
        outputs = [
            file_record(out_chr16_path, record_count=n_kept_chr16),
            file_record(out_emb_path, record_count=n_kept_emb),
        ]

        metrics = {
            "n_files_in_manifest": n,
            "n_files_kept_embeddings": n_kept_emb,
            "n_files_kept_chr16": n_kept_chr16,
            "n_files_unk_only": n_unk_only,
            "n_files_no_tokens": n_no_tokens,
            "n_files_failed": n_failed,
            "chr16_universe_size": len(chr16_set),
            "embedding_dim": emb_dim,
            "batch_size": batch_size,
            "unk_only_samples": unk_only_samples[:10],
            "failed_samples": [f"{bid}: {msg}" for bid, msg in failed_samples],
        }
        write_summary(ctx, "success" if n_failed == 0 else "partial",
                      inputs=inputs, outputs=outputs, metrics=metrics)
        print(f"\n  emb: {n_kept_emb:,}/{n:,}, chr16: {n_kept_chr16:,}/{n:,}",
              file=sys.stderr)

    except Exception as e:
        write_summary(ctx, "error", inputs=inputs,
                      error={"class": type(e).__name__, "message": str(e),
                             "suggestion": "Check BBCLIENT_CACHE; pretrained_universe.parquet "
                                           "must have an 'embedding' column (rerun stage 07 if missing).",
                             "retryable": True})
        raise


if __name__ == "__main__":
    main()
