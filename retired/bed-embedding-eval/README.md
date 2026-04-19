# bed-embedding-eval (retired 2026-04-19)

Evaluation framework for BED file embeddings (dense / sparse / retrieval tracks). **Retired because the dense-track finding that motivated the broader spatial-learning program is methodologically confounded: R2V embedding similarity correlates with file-level geometric features, but the correlation cannot be separated cleanly from biological similarity. Files that share a regulatory function co-occur in training *and* share characteristic geometric signatures, so the analysis cannot distinguish "R2V learns geometry" from "R2V learns biology that has geometric signatures."** Superseded by within-universe content analysis in the parent `spatial-region-features/` workspace.

## Why we ran this

Sheffield Lab [discussion #77](https://github.com/databio/lab.databio.org/discussions/77) (dense and sparse embeddings for enrichment) and [discussion #69](https://github.com/databio/lab.databio.org/discussions/69) (sparse retrieval for genomic intervals) asked whether BED file embeddings could support enrichment, retrieval, and visualization for [BEDbase](https://bedbase.org). The project was scoped as three tracks:

- **Dense** — what do Region2Vec embeddings encode? When are they reliable?
- **Sparse** — can BM25 sparse embeddings approximate LOLA enrichment without a universe?
- **Retrieval** — which method ranks best for bed2bed similarity search?

## What got done

**Dense track (complete).** R2V embedding similarities across BEDbase hg38 files were correlated against per-file geometric features — region count, mean region width, GC content, inter-peak spacing — computed from the raw BED files. Neighbor-correlation analyses were run on train / test / val splits with hand-built ground-truth sets.

**Sparse track (partial).** BM25 sparse-vector embeddings were produced via `gtars` (from [PR #239](https://github.com/databio/gtars/pull/239) at the time of the work). The comparison against LOLA enrichment was scaffolded — run scripts and SLURM orchestration existed — but end-to-end results were not produced.

**Retrieval track (scaffolded).** Framework and `tourney_eval/` suite were set up; no comparative retrieval experiments were completed end-to-end.

## What we found

**Dense track.** Neighbor-correlation ρ 0.88–0.98 between R2V embedding similarity and file-level geometric features (region count, mean width, GC). Conclusion at the time: "physical distribution dominates R2V's signal." This result motivated the original Chapter 2 framing about "sub-token spatial features beyond what R2V captures."

**Sparse and retrieval tracks** did not produce final results.

## Why we retired it

The dense-track conclusion is methodologically confounded. The inference chain was:

1. R2V embeddings cluster similar BED files.
2. BED files that cluster together also share file-level geometric statistics.
3. Therefore R2V is primarily learning geometry.

Step 3 does not follow. BED files that share a regulatory function (e.g., multiple ChIP-seq experiments targeting promoters) both co-occur in training through shared tokens *and* share characteristic geometric signatures (promoters are dense, GC-rich, narrow-peaked). The observed correlation is equally compatible with "R2V learns geometry" and with "R2V learns biology that has geometric signatures." The analysis design cannot separate these.

The interpretability question this project was motivated by — *what do R2V and Atacformer embeddings actually represent?* — is real and still unresolved. But the right approach is to characterize the **content** of what lives inside known regulatory-element intervals (the universe-as-scope direction), not to correlate file-level embedding similarity against file-level geometric statistics.

## What to carry forward, what not to

**Carry forward:**
- The embedding-interpretability motivation itself. "Token co-occurrence gives similarity but not semantics" is the pain point, and it is real.
- The data-pipeline shape (BED downloads from BEDbase, metadata joins, yoke/brickyard sync). Useful scaffolding, modulo rework.

**Do not carry forward:**
- The file-level geometric-correlation analysis as evidence for what embeddings encode. The conclusion is confounded.
- The three-track design (dense / sparse / retrieval). The sparse and retrieval tracks were never completed and reflect an older framing where the sharp question had not yet been articulated.
- The notion that correlating embeddings against per-file features tells you what the embedding represents. It doesn't — it tells you what shares a variance-explaining factor, which could be biology, geometry, or both.
