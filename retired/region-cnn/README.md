# region-cnn (retired 2026-04-19)

Exploration of learned 1D spatial features from genomic region data. **Retired because the methodology was bottom-up rather than biology-anchored — the spatial features were computed over genome-wide groupings (peak clusters at various radii, density windows) that do not correspond to known regulatory units, and aggregations over those groupings cannot be interpreted biologically.** Superseded by the universe-anchored approach in the parent `spatial-region-features/` workspace.

## Why we ran this

Sheffield Lab [discussion #79](https://github.com/databio/lab.databio.org/discussions/79) (Learned Spatial Features for Genomic Region Data) asked whether a 1D CNN could learn interpretable "region motifs" from binary BED data, without a predefined universe. Question 1 from that discussion: *what are the "region motifs" of the genome?*

The project was scoped in two stages:

- **Stage 0 — non-DL spatial-features baseline** as a falsification gate. If hand-crafted spatial features did not beat density-only features on classification tasks, the CNN was not worth building.
- **Stage 1 — self-supervised CNN encoder-decoder** with span masked-region modeling as the main investigation.

## What got done

**Stage 0 (complete).** A `gtars-genomicdist` feature set was extracted per BED file: `inter_peak_spacing`, `peak_clusters` at multiple radii, `density_vector`, `density_homogeneity`, and an updated `gaps()` function. A per-fold residualization procedure separated spatial features from density-only features to isolate orthogonal signal. Six classification experiments were run on ENCODE (n ≈ 744) and ChIP-Atlas (n ≈ 13,000) BED files covering assay type, cell type / cell type class, biosample name, histone mark, and ChIP target labels.

**Stage 1 (not started).** The CNN encoder-decoder with span-MLM pretraining was planned (see the 2026-04-10 Q1 filter-interpretability plan in the archived source) but never implemented.

## What we found

Residualized spatial features added linear-classifier signal above density-only features across all six experiments. The effect was largest on within-assay cell-type classification — the regime where density features cannot win by exploiting assay differences. Original result artifacts (`classifier_results_*.parquet`) live in the analysis repo's `sam` branch at `analysis/region-cnn/results/baseline/`.

A separate downstream contribution: the `gtars-genomicdist` Rust library with the new spatial-statistics functions was contributed upstream to `databio/gtars`. That infrastructure contribution is durable and independent of this project's retirement.

## Why we retired it

Two methodological problems surfaced during the universe-as-scope framing work in April 2026:

1. **Bottom-up region definition.** The spatial features were computed over regions constructed bottom-up from the BED file — peak clusters at various radii, density windows — that have no biological identity. These constructed regions are not known promoters, enhancers, CTCF insulators, or any other regulatory element. Aggregating spatial features over groupings of unknown biology produces numbers that are statistically valid but semantically empty: "mean inter-peak spacing across stitched groupings" does not translate into any biological statement.

2. **The "spatial beats density" finding is narrower than it looked.** If the features themselves measure bottom-up geometric proxies rather than biology, the finding reduces to "this specific bottom-up spatial proxy adds signal above density" — much narrower than "spatial structure in regulatory biology carries information beyond density."

The fix is to anchor spatial feature computation to **known regulatory-element intervals** (ENCODE cCREs, tissue-matched consensus peaks, or a validated universe) rather than genome-wide bottom-up groupings. That work continues in the parent `spatial-region-features/` workspace under the universe-as-scope framing.

## What to carry forward, what not to

**Carry forward:**
- The idea that hand-engineered spatial features can exceed a density baseline. This is still plausibly true under a better feature definition; the Stage 0 result is suggestive evidence, not conclusive.
- The general infrastructure shape (per-fold residualization, multi-assay comparison panels).

**Do not carry forward:**
- `peak_clusters` stitching at arbitrary radii as a feature.
- Any feature computed over genome-wide groupings that do not correspond to known regulatory elements.
- The "no universe needed" framing — the universe-as-scope work in this workspace directly contradicts it.
