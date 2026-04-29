"""Stage 10 — featured-narrative data prep.

Supports steps 1/2/4 of the viz narrative arc:

  Step 1 — featured BED files at featured intervals (peaks per file per interval)
  Step 2 — universe regions per interval; per-(file, interval) binary activation
  Step 4 — per-file chr16 token activation, for highlighting a clicked file's
            regions on the region UMAP

Curates ~18 BED files across 3 viz tissues × 6 marks (one per pair, drawn from
the rich-metadata Tight subset), plus a few mystery (UNKNOWN cell_line) entries.

Tokenizes each BED against the **pretrained R2V model's universe** (same one
stage 07 unpacks into pretrained_universe.parquet) — so token_ids returned
here align with the parquet, and Step 2's universe band lines up correctly
with file activations.

Reads:
  data/corpus/manifest.parquet                       (from stage 01)
  data/precomputed/viz_files.parquet                 (from stage 09; for mystery sample)
  data/annotations/pretrained_universe.parquet       (from stage 07)
  databio/r2v-encode-hg38                            (HuggingFace; same as stage 07)
  BBCLIENT_CACHE env                                 (Rivanna only)
Writes:
  data/precomputed/featured_intervals.parquet        — per-interval narrative + universe slice
  data/precomputed/featured_files.parquet            — per-file chr16 token activation
  data/precomputed/featured_tracks.parquet           — per-(file, interval) peaks + active tokens
  results/10_featured_narrative/summary.json
"""
from __future__ import annotations

import os
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


def select_featured_files(
    manifest: pl.DataFrame,
    cell_lines: list[str],
    targets: list[str | None],
    skip_cells: list[tuple[str, str]],
    seed: int,
) -> list[dict]:
    """Pick one BED file per (cell_line, target) with rich metadata.

    Targets: list of values like "H3K4me1" / "H3K27ac" / None (= ATAC-seq).
    ATAC-seq files have target = NULL; histone ChIP files have target set.
    skip_cells: (cell_line, target) pairs known to have zero files in the
    source pool (e.g., GM12878 H3K9me3) — drop from the grid silently.
    Sampling deterministically: pick the file with the smallest id within each group.
    """
    skip_set = {(cl, tgt) for cl, tgt in skip_cells}
    selected: list[dict] = []
    for cl in cell_lines:
        for tgt in targets:
            if tgt is not None and (cl, tgt) in skip_set:
                continue
            if tgt is None:
                # ATAC-seq: filter by assay, ignore target
                pool = manifest.filter(
                    (pl.col("cell_line") == cl)
                    & (pl.col("assay") == "ATAC-seq")
                )
            else:
                pool = manifest.filter(
                    (pl.col("cell_line") == cl)
                    & (pl.col("target") == tgt)
                    & (pl.col("assay").str.contains("ChIP"))
                )
            if len(pool) == 0:
                print(f"  WARN: no file found for ({cl}, {tgt or 'ATAC'})", file=sys.stderr)
                continue
            row = pool.sort("id").head(1).to_dicts()[0]
            row["_target_label"] = tgt or "ATAC-seq"
            row["_role"] = "featured"
            selected.append(row)
    return selected


def select_mystery_files(viz_files: pl.DataFrame, n: int, seed: int) -> list[dict]:
    """Pick n curated mystery (UNKNOWN cell_line) entries from viz_files.parquet."""
    pool = viz_files.filter(pl.col("is_unlabeled"))
    if len(pool) == 0:
        return []
    if len(pool) <= n:
        sample = pool
    else:
        sample = pool.sample(n=n, seed=seed)
    rows = []
    for r in sample.to_dicts():
        r["_target_label"] = r.get("assay") or "UNKNOWN"
        r["_role"] = "mystery"
        rows.append(r)
    return rows


def tokens_in_interval(
    universe: pl.DataFrame, chrom: str, start: int, end: int,
) -> list[int]:
    """Return token_ids whose intervals overlap [start, end) on chrom."""
    sub = universe.filter(
        (pl.col("chrom") == chrom)
        & (pl.col("end") > start)
        & (pl.col("start") < end)
    )
    return sub["token_id"].to_list()


def main() -> None:
    ctx = stage_start("10_featured_narrative", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    manifest_path = PROJECT_ROOT / paths["corpus_dir"] / "manifest.parquet"
    viz_files_path = PROJECT_ROOT / paths["precomputed_dir"] / "viz_files.parquet"
    pretrained_path = PROJECT_ROOT / paths["annotations_dir"] / "pretrained_universe.parquet"
    precomp_dir = PROJECT_ROOT / paths["precomputed_dir"]

    out_intervals = precomp_dir / "featured_intervals.parquet"
    out_files = precomp_dir / "featured_files.parquet"
    out_tracks = precomp_dir / "featured_tracks.parquet"

    # Use the same pretrained model as stage 07 so tokenizer vocabulary IDs
    # align with the IDs in pretrained_universe.parquet.
    hf_model = sc.get("hf_model", "databio/r2v-encode-hg38")

    inputs = [
        file_record(manifest_path),
        file_record(viz_files_path),
        file_record(pretrained_path),
    ]

    bbclient_cache = os.environ.get("BBCLIENT_CACHE")
    try:
        if not bbclient_cache:
            raise RuntimeError("BBCLIENT_CACHE env var not set. source the lab env.sh first.")
        for p in [manifest_path, viz_files_path, pretrained_path]:
            if not p.exists():
                raise FileNotFoundError(f"{p} missing — check upstream stages.")

        intervals_cfg = sc["featured_intervals"]
        cell_lines = sc["featured_cell_lines"]
        targets = list(sc["featured_targets"])  # may contain None
        skip_cells = [tuple(p) for p in sc.get("skip_cells", [])]
        n_mystery = int(sc.get("n_mystery_files", 4))
        seed = int(sc.get("mystery_seed", 42))

        # --- Pick featured files -------------------------------------------------
        manifest = pl.read_parquet(manifest_path)
        viz_files = pl.read_parquet(viz_files_path)
        featured = select_featured_files(manifest, cell_lines, targets, skip_cells, seed)
        mystery = select_mystery_files(viz_files, n_mystery, seed)
        all_files = featured + mystery
        print(
            f"  selected {len(featured)} featured + {len(mystery)} mystery files",
            file=sys.stderr,
        )

        # --- Universe lookup table for slicing -----------------------------------
        universe = pl.read_parquet(pretrained_path).select([
            "token_id", "chrom", "start", "end", "cclass", "overlaps_screen",
        ])

        # --- Per-interval token slice -------------------------------------------
        interval_rows = []
        interval_token_lookup: dict[str, set[int]] = {}
        for iv in intervals_cfg:
            tids = tokens_in_interval(universe, iv["chrom"], iv["start"], iv["end"])
            interval_token_lookup[iv["id"]] = set(tids)
            interval_rows.append({
                "interval_id": iv["id"],
                "chrom": iv["chrom"],
                "start": iv["start"],
                "end": iv["end"],
                "label": iv["label"],
                "narrative_caption": iv["narrative_caption"],
                "n_universe_tokens": len(tids),
                "universe_token_ids": tids,
            })
        intervals_df = pl.DataFrame(interval_rows)
        precomp_dir.mkdir(parents=True, exist_ok=True)
        intervals_df.write_parquet(out_intervals)

        # --- Tokenize each featured/mystery file against the PRETRAINED universe.
        # Critical: must use model.tokenizer (the pretrained vocab), not a
        # separate Tokenizer loaded from our SCREEN cCRE BED — those would be
        # different ID spaces and Step 2's universe band wouldn't align with
        # file activations.
        from geniml.bbclient import BBClient
        from geniml.region2vec.main import Region2VecExModel

        print(f"  loading pretrained tokenizer from {hf_model}...", file=sys.stderr)
        model = Region2VecExModel(model_path=hf_model)
        tokenizer = model.tokenizer
        unk_id = tokenizer.unk_token_id
        bbc = BBClient(cache_folder=bbclient_cache)

        # Pre-build a chr16 universe lookup: token_id → (start, end) so we can
        # filter activation lists to chr16 only.
        chr16_token_set = set(
            universe.filter(pl.col("chrom") == "chr16")["token_id"].to_list()
        )

        file_rows = []
        track_rows = []
        for f in all_files:
            bed_id = f["id"]
            print(f"  tokenizing {f.get('_role'):8s} {bed_id[:12]}...", file=sys.stderr)
            t0 = time.time()
            try:
                rs = bbc.load_bed(bed_id)
            except Exception as e:
                print(f"    LOAD FAILED: {e}", file=sys.stderr)
                continue

            try:
                # Tokenize against universe → chr16 activation set
                token_ids = list(tokenizer.encode(tokenizer.tokenize(rs)))
                token_ids = [t for t in token_ids if t != unk_id]
                chr16_active = sorted(set(t for t in token_ids if t in chr16_token_set))
            except Exception as e:
                print(f"    TOKENIZE FAILED: {e}", file=sys.stderr)
                continue

            # Pull peak coords as chrom/start/end triples (only chr16) for tracks.
            peak_chr16 = []
            try:
                for region in rs:
                    chrom = getattr(region, "chr", None) or getattr(region, "chrom", None)
                    if chrom == "chr16":
                        peak_chr16.append((int(region.start), int(region.end)))
            except Exception as e:
                print(f"    WARN: failed to read regions for tracks ({e})", file=sys.stderr)

            file_rows.append({
                "file_id": bed_id,
                "name": f.get("name") or "",
                "assay": f.get("assay") or "",
                "cell_line": f.get("cell_line") or "",
                "target": f.get("target") or f.get("_target_label") or "",
                "role": f["_role"],
                "n_chr16_active_tokens": len(chr16_active),
                "chr16_active_token_ids": chr16_active,
            })

            # Per-interval slicing.
            for iv in intervals_cfg:
                iv_token_set = interval_token_lookup[iv["id"]]
                active_in_iv = sorted(set(t for t in chr16_active if t in iv_token_set))
                peaks_in_iv = [
                    (s, e) for (s, e) in peak_chr16
                    if e > iv["start"] and s < iv["end"]
                ]
                if not active_in_iv and not peaks_in_iv:
                    continue
                track_rows.append({
                    "file_id": bed_id,
                    "interval_id": iv["id"],
                    "n_active_tokens": len(active_in_iv),
                    "active_token_ids": active_in_iv,
                    "n_peaks": len(peaks_in_iv),
                    "peak_starts": [s for (s, _) in peaks_in_iv],
                    "peak_ends": [e for (_, e) in peaks_in_iv],
                })

            print(
                f"    chr16: {len(chr16_active):,} active tokens, {len(peak_chr16):,} peaks "
                f"({time.time() - t0:.1f}s)",
                file=sys.stderr,
            )

        files_df = pl.DataFrame(file_rows)
        tracks_df = pl.DataFrame(track_rows)
        files_df.write_parquet(out_files)
        tracks_df.write_parquet(out_tracks)

        outputs = [
            file_record(out_intervals, record_count=len(intervals_df)),
            file_record(out_files, record_count=len(files_df)),
            file_record(out_tracks, record_count=len(tracks_df)),
        ]

        metrics = {
            "n_featured_intervals": len(intervals_df),
            "n_featured_files_selected": len(featured),
            "n_mystery_files_selected": len(mystery),
            "n_files_tokenized": len(files_df),
            "n_track_records": len(tracks_df),
            "interval_token_counts": {
                r["interval_id"]: int(r["n_universe_tokens"])
                for r in interval_rows
            },
            "per_file_chr16_token_counts": files_df.select([
                "file_id", "role", "n_chr16_active_tokens"
            ]).to_dicts(),
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(f"\n  wrote 3 parquets to {precomp_dir}", file=sys.stderr)

    except Exception as e:
        write_summary(
            ctx, "error",
            inputs=inputs,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": (
                    "If BBCLIENT_CACHE missing, source /project/shefflab/rivanna_config/env.sh. "
                    "If a featured file fails to load, check rich-metadata coverage in the manifest "
                    "(some (cell_line, target) pairs may not have a representative)."
                ),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
