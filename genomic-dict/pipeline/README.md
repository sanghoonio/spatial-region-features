# Pipeline runner notes

The 2026-04-29 cleanup pruned the production pipeline to **10 stages** with
sequential numbering 00–09 reflecting actual execution order. Bigwig
dependencies were dropped from embedding-side computations; analysis-only
stages moved to `scripts/offline/`. See
`plans/2026-04-29-bed-route-precompute-strip.md` for design rationale.

## Production stages (sequential execution order)

| # | Script | Slurm wrapper | Purpose |
|---|---|---|---|
| 00 | `scripts/00_inspect_metadata.py` | `slurm/00_02_prep.sh` | Download HF metadata parquets + render distribution summary. Required: stage 01 reads the downloaded parquets. |
| 01 | `scripts/01_curate_corpus.py` | `slurm/00_02_prep.sh` | Filter ENCODE → 17k balanced manifest. |
| 02 | `scripts/02_prepare_universe.py` | `slurm/00_02_prep.sh` | Download SCREEN V4 cCRE bed + class metadata. |
| 03 | `scripts/03_fetch_bigwigs.py` | `slurm/03_fetch_bigwigs.sh` | Download 21 bigwigs (3 cells × 7 marks). Bigwigs are now consumed only by stage 09. |
| 04 | `scripts/04_target_evidence.py` | `slurm/04_target_evidence.sh` | ENCODE V4 cCRE-Gene Links → target evidence parquet. |
| 05 | `scripts/05_load_pretrained.py` | `slurm/05_load_pretrained.sh` | Load pretrained R2V model + intersect universe with SCREEN. |
| 06 | `scripts/06_tokenize_corpus_chr16.py` | `slurm/06_tokenize_corpus_chr16.sh` | Tokenize the 17k corpus + emit chr16 token lists, file embeddings, **and the file UMAP** (former stage 09 folded in). |
| 07 | `scripts/07_precompute_viz.py` | `slurm/07_precompute_viz.sh` | Region UMAP + kNN + **BED-derived concept axes** (former `embedding_features.py`). Must run **after** stage 06. |
| 08 | `scripts/08_featured_narrative.py` | `slurm/08_featured_narrative.sh` | Featured intervals + featured files + featured tracks. |
| 09 | `scripts/09_featured_signal.py` | `slurm/09_featured_signal.sh` | Bigwig sampling at the 4 featured intervals (Section 1 continuous mode). |

### Order constraints

- 00 → 01 → 02 → 05 (the universe + class join chain)
- 03 (bigwigs) → 09
- 01 → 06 (tokenizes against manifest using R2V tokenizer from 05)
- 05 + 06 → 07 (07 needs both the universe-with-embeddings and the tokenized corpus for BED-derived concept axes)
- 01 → 08 → 09 (featured narrative drives bigwig sampling)
- 04 is independent of 05–07 — can run in parallel.

A reasonable submission order on Rivanna (using `--dependency=afterok:<jobid>` to chain):

```
sbatch slurm/00_02_prep.sh                      # 00 + 01 + 02
sbatch slurm/03_fetch_bigwigs.sh
sbatch slurm/04_target_evidence.sh              # independent of 05+
sbatch slurm/05_load_pretrained.sh
sbatch slurm/06_tokenize_corpus_chr16.sh        # ~30-45 min, after 05
sbatch slurm/07_precompute_viz.sh               # ~5 min, after 06
sbatch slurm/08_featured_narrative.sh           # after 01
sbatch slurm/09_featured_signal.sh              # after 03 + 08
```

## Offline stages (analysis only; outputs not shipped to UI)

`scripts/offline/`:

- `12_cooccurrence_pmi.py` — NPMI partners per stratum (slurm/offline/12_cooccurrence_pmi.sh).
- `13_modules.py` — Leiden communities per stratum.

These are kept for offline analysis (canvas / canvas2 hub-discovery work uses
them), but neither output ships to the React demo's
`public/data/dictionary/`. The in-browser DuckDB-WASM coordinator computes
NPMI on demand from the shipped `tokenized_corpus_chr16.parquet`. Numbering
inside `offline/` is left at 12/13 for continuity with prior analysis logs.

## Dev / probe scripts

`scripts/dev/`:

- `preflight_r2v_faithfulness.py`
- `probe_token_cache.py`, `probe_token_cache_v2.py`

Diagnostic tools, not part of the production flow.

## Removed stages (2026-04-29 cleanup)

- **04 extract_intrinsic** (legacy number) — produced per-(token × bigwig) signal means; only consumer was old stage 08, which now reads `pretrained_universe` directly.
- **05 extract_sequence** (legacy number) — sequence features (FIMO, ChromHMM, phastCons). Never wired up downstream.
- **09 prepare_file_viz** (legacy number) — folded into the new stage 06. The per-file embeddings + manifest join + UMAP all happen in one pass now.
- **`embedding_features.py`** — folded into the new stage 07. Concept axes case/ctrl masks now use BED peak counts instead of bigwig signal means, keeping embedding-side computations grounded in R2V's training substrate.

## Files shipped to the React UI

After running the pipeline above, copy the following from
`data/precomputed/` into `genomic-regions/public/data/dictionary/`:

- `viz_chr16.parquet`
- `viz_files.parquet`
- `tokenized_corpus_chr16.parquet`
- `region_concept_axes.parquet`
- `region_target_evidence.parquet`
- `featured_intervals.parquet`
- `featured_files.parquet`
- `featured_tracks.parquet`
- `featured_signal.parquet`

**Do not ship**: `region_cooccurrence_pmi.parquet`, `region_stratum_marginals.parquet`, `region_modules.parquet`, `module_summary.parquet`, `region_class_prototypes.parquet`, `region_target_evidence_summary.parquet`, `region_stats.parquet`, `file_embeddings.parquet` (intermediate). These are either offline-only outputs or replaced by on-demand DuckDB queries.
