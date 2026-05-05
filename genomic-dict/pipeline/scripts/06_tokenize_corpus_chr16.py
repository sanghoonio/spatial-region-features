"""Stage 11 — tokenize the curated corpus + compute the file UMAP.

Single-pass tokenization against the pretrained R2V universe (~1.06M tokens,
all chromosomes). For each BED file emits three artifacts:

  - chr16-active token list — CHR16-ONLY.
      Tokens filtered to focus_chromosome (chr16). Drives the chr16-region
      viz, on-demand NPMI, and offline cooccurrence stages.
      → data/precomputed/tokenized_corpus_chr16.parquet

  - 100-dim mean embedding — GENOME-WIDE.
      Mean of per-token embeddings over every non-UNK token in the file.
      Used internally by this stage to compute the file UMAP; persisted as
      an intermediate for offline analysis.
      → data/precomputed/file_embeddings.parquet

  - File UMAP — file-level 2D coords + manifest metadata.
      Final viz-ready table joining the manifest with file embeddings + UMAP
      coords. Replaces the prior stage 09 (which is now folded in).
      → data/precomputed/viz_files.parquet

R2V is Word2Vec-derived (no doc vector), so the genome-wide mean is
mathematically equivalent to model.encode(rs, pooling="mean").mean(axis=0).

Streaming writers for the per-file artifacts, batch_size to bound memory.
File UMAP is computed once at the end (~17k × 100-dim, trivial cost).

Reads:
  data/corpus/manifest.parquet                       (from stage 01)
  data/annotations/pretrained_universe.parquet       (from stage 07)
  databio/r2v-encode-hg38                            (HuggingFace; for tokenizer)
  BBCLIENT_CACHE env                                 (Rivanna only)
Writes:
  data/precomputed/tokenized_corpus_chr16.parquet    — id, chr16_active_token_ids, n_chr16_active
  data/precomputed/file_embeddings.parquet           — id, embedding (list[float], 100-dim)
  data/precomputed/viz_files.parquet                 — id + manifest + umap_x/y
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
    ctx = stage_start("06_tokenize_corpus_chr16", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    manifest_path = PROJECT_ROOT / paths["corpus_dir"] / "manifest.parquet"
    pretrained_path = PROJECT_ROOT / paths["annotations_dir"] / "pretrained_universe.parquet"
    precomp_dir = PROJECT_ROOT / paths["precomputed_dir"]
    out_chr16_path = precomp_dir / "tokenized_corpus_chr16.parquet"
    out_emb_path = precomp_dir / "file_embeddings.parquet"
    out_files_path = precomp_dir / "viz_files.parquet"

    # File UMAP knobs (folded in from former stage 09 — same defaults).
    files_n_neighbors = int(sc.get("files_umap_n_neighbors", 15))
    files_min_dist = float(sc.get("files_umap_min_dist", 0.1))
    files_seed = int(sc.get("files_umap_seed", 42))

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

        # ---- File UMAP (folded from former stage 09) ----
        # Read back the just-written file_embeddings, join manifest, run UMAP,
        # write viz_files.
        print("\n  computing file UMAP...", file=sys.stderr)
        manifest_df = pl.read_parquet(manifest_path)
        emb_df = pl.read_parquet(out_emb_path)
        joined = manifest_df.join(
            emb_df.select(["id", "embedding"]), on="id", how="inner",
        )
        n_joined = len(joined)
        if n_joined == 0:
            raise RuntimeError("No id overlap between manifest and file embeddings.")
        print(f"    joined {n_joined:,} files", file=sys.stderr)

        emb_matrix = np.asarray(joined["embedding"].to_list(), dtype=np.float32)
        import umap
        t_umap0 = time.time()
        reducer = umap.UMAP(
            n_neighbors=files_n_neighbors,
            min_dist=files_min_dist,
            n_components=2,
            random_state=files_seed,
            metric="cosine",
        )
        coords = reducer.fit_transform(emb_matrix)
        files_umap_seconds = round(time.time() - t_umap0, 1)
        print(f"    file umap {files_umap_seconds}s; coords shape {coords.shape}", file=sys.stderr)

        out_files = (
            joined.with_columns([
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
            .sort("id")
        )
        out_files.write_parquet(out_files_path)
        print(
            f"    wrote {out_files_path.name} ({out_files_path.stat().st_size / 1e6:.1f} MB)",
            file=sys.stderr,
        )

        outputs = [
            file_record(out_chr16_path, record_count=n_kept_chr16),
            file_record(out_emb_path, record_count=n_kept_emb),
            file_record(out_files_path, record_count=len(out_files)),
        ]

        # Distributional metrics for the file UMAP — same shape as old stage 09.
        assay_counts = (
            out_files.group_by("assay").agg(pl.len().alias("n"))
            .sort("n", descending=True).to_dicts()
        )
        cell_line_counts = (
            out_files.group_by("cell_line").agg(pl.len().alias("n"))
            .sort("n", descending=True).head(15).to_dicts()
        )

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
            # File UMAP (folded stage 09):
            "n_files_umap": len(out_files),
            "files_umap_seconds": files_umap_seconds,
            "files_umap_params": {
                "n_neighbors": files_n_neighbors,
                "min_dist": files_min_dist,
                "seed": files_seed,
                "metric": "cosine",
            },
            "assay_counts": assay_counts,
            "top_cell_lines": cell_line_counts,
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
