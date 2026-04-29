---
date: 2026-04-28
status: complete
description: Methods reference for the region UMAP color encodings in genomic-regions. Per-toggle data sources, formulas, color schemes, range bounds, and interpretation caveats.
---

# Region UMAP color encodings — methods reference

The region UMAP in `genomic-regions/src/index.md` (Section 2) has a `Color regions by` radio with **9 options** (3 inherited from earlier viz, 6 added 2026-04-28). This doc says how each encoding is computed, where the data comes from, and what to read into it.

For all encodings, the underlying scatter is **`classedChr16`** — chr16 tokens with non-null `cclass`, augmented with derived fields. UMAP coords (`umap_x`, `umap_y`) come from `viz_chr16.parquet` (stage 08, run on chr16 R2V embeddings via `umap-learn`).

## Source layer summary

| Encoding | Underlying parquet | Pipeline stage |
|---|---|---|
| SCREEN class | `viz_chr16.parquet` | stage 08 |
| Genomic midpoint | `viz_chr16.parquet` | stage 08 (start/end inline) |
| File count (log) | `region_stats.parquet` | retired stage 12 v1 (kept as orphan) |
| Anchor score | `region_concept_axes.parquet` | `embedding_features.py` |
| Activity score | `region_concept_axes.parquet` | `embedding_features.py` |
| K562 / GM12878 / HepG2 specificity | `region_concept_axes.parquet` | `embedding_features.py` |
| Target evidence count (log) | `region_target_evidence_summary.parquet` | stage 06 |

All concept-axis values are derived from the **pretrained R2V model** (`databio/r2v-encode-hg38`, 100-dim embeddings per region) via cosine projection onto direction vectors; see `embedding_features.py` for the full code.

---

## 1. SCREEN class (`cclass`)

**Source.** `viz_chr16.parquet` column `cclass`. Inherited from the SCREEN registry of cCREs by intersecting each chr16 R2V universe region with overlapping cCREs and keeping the largest-overlap class assignment (stage 07 for the universe-wide tagging; stage 08 carries it into `viz_chr16`).

**Values.** Categorical, one of:
- `PLS` — promoter-like signature (n=1,770 on chr16)
- `pELS` — proximal enhancer-like (n=7,475)
- `dELS` — distal enhancer-like (n=25,055; 70% of chr16 universe)
- `CA-CTCF` — chromatin-accessible CTCF-bound (n=899)
- `CA-H3K4me3` — chromatin-accessible H3K4me3-marked (n=600)
- `unclassed` — no SCREEN overlap (filtered out of `classedChr16`)

**Color.** Explicit palette from `genomic-regions/src/components/dictionary.js`:

```
PLS         #ff0000  red
pELS        #ffa700  orange
dELS        #ffcd00  yellow
CA-CTCF     #00b0f0  blue
CA-H3K4me3  #ffaaaa  pink
```

**Caveat.** Hard categorical labels with no sub-structure. The dELS class is 70% of the universe and visually dominates as one yellow blob. Sub-class structure (e.g., lineage-specific dELS subpopulations) is invisible here — that's what the concept-axis encodings expose.

---

## 2. Genomic midpoint (`midpoint`)

**Source.** Derived inline in the viz from `viz_chr16.start` and `viz_chr16.end`:

```
midpoint = (Number(start) + Number(end)) / 2
```

**Values.** Continuous; range covers chr16 (~0–90M bp).

**Color.** Turbo scheme, label formatted as Mb (`(d / 1e6).toFixed(1)`).

**Interpretation.** This is a **negative diagnostic**. If R2V had memorized chromosomal position, you'd see a smooth turbo gradient across the UMAP. Instead the color is uniformly mixed — confirming the embedding learned function rather than position. Useful as a sanity check; not as a feature.

---

## 3. File count (log) (`n_files_total`)

**Source.** Stale `region_stats.parquet` from the retired stage 12 v1 (corpus stats with nearest-gene approach). Kept in the bundle as an orphan per Locked-in decision #1 in the pipeline plan; will be replaced once the new stage 12 PMI parquet is wired in.

**Definition (current stale version).** For each chr16 token, the count of BED files in the 16,799-file curated corpus whose `chr16_active_token_ids` list contains this token. Equivalent to the per-(token, `corpus_baseline` stratum) `n_files_active` value that stage 12 v3's `region_cooccurrence_pmi.parquet` carries — **the swap target when the orphan is retired**.

**Color.** Viridis sqrt scale (`vg.colorScale("sqrt")`). Sqrt compresses the long tail without going fully log; a token in 5,000 files is visually distinguishable from one in 100, but the gradient stays usable.

**Interpretation.** Tier D "popularity" proxy. Bright = housekeeping / broadly-accessible (CpG-island promoters, ubiquitous open chromatin). Dark = cell-type-specific or rarely-called. Doesn't tell you what kind of regulation; only how often the token shows up.

---

## 4. Anchor score (`anchor_score`)

**Source.** `region_concept_axes.parquet`, column `anchor_score`. Produced by `embedding_features.py`.

**Definition.** Cosine projection of each region's R2V embedding onto a learned **promoter-vs-enhancer direction vector** in 100-dim space:

```
anchor_axis = mean(R2V[t] for t in PLS-tokens) − mean(R2V[t] for t in dELS-tokens)
anchor_score(t) = cos(R2V[t], anchor_axis)
```

The "PLS centroid" and "dELS centroid" are the unweighted means of the chr16 R2V embeddings for tokens with `cclass == "PLS"` (n=1,770) and `cclass == "dELS"` (n=25,055). The axis vector is their difference.

**Range.** Bounded [-1, +1] (cosine). Empirical range ~[-0.7, +0.7] in practice.

**Color.** RdBu divergent. Domain `[-0.7, +0.7]` (symmetric). Blue = positive (PLS-like); red = negative (dELS-like). White = ~0 (no preference).

**Validation.** Per-class median anchor_score from the script's sanity check:

| cclass | median anchor_score |
|---|---|
| PLS | **+0.428** |
| pELS | +0.040 |
| CA-H3K4me3 | −0.189 |
| dELS | −0.212 |
| CA-CTCF | −0.249 |

Monotonically PLS > pELS > CA-* / dELS — clean. The embedding has internalized a continuous promoter-vs-distal-enhancer direction without supervision.

**Interpretation.** Direct visual evidence the embedding learned the promoter ↔ distal-enhancer gradient as a *single direction* in 100-dim space. The UMAP shows it as a horizontal gradient (right = blue/PLS-rich, left = red/dELS-rich).

---

## 5. Activity score (`activity_score`)

**Source.** `region_concept_axes.parquet`, column `activity_score`. Produced by `embedding_features.py`.

**Definition.** Cosine projection onto a **K562-active vs K562-repressed direction**:

```
case_set = top 25% of chr16 tokens by max(K562 H3K27ac, K562 H3K4me1, K562 H3K4me3)
control_set = top 25% by K562 H3K27me3
activity_axis = mean(R2V[t] for t in case_set) − mean(R2V[t] for t in control_set)
activity_score(t) = cos(R2V[t], activity_axis)
```

The signal-mark columns come from `viz_chr16.parquet`'s 18 bigwig-mean columns (stage 04). Top-25% thresholds are computed via numpy `percentile(75)`.

**Range.** Bounded [-1, +1]. Empirical range ~[-0.5, +0.5].

**Color.** PiYG divergent. Domain `[-0.5, +0.5]`. Green = active-like, magenta = repressed-like.

**Why K562 specifically?** The threshold sets `case_set` and `control_set` are K562-defined, so this axis is technically a K562-active-vs-repressed gradient, not a generic activity gradient. In practice it generalizes — most regions that are H3K27ac-active in K562 are also H3K27ac-active broadly — but the axis carries a slight K562 bias by construction.

**Interpretation.** Visible as a roughly **vertical** gradient in the UMAP (active-like at the bottom, repressed-like at the top), **orthogonal to the anchor axis's left-right gradient**. The two concept axes are picking up genuinely independent dimensions: anchor describes class identity (promoter vs enhancer), activity describes chromatin state (active vs repressed). A region can be promoter-like AND repressed (Polycomb-bivalent), or distal-enhancer-like AND active (typical pELS) — four quadrants in the UMAP, each with biological meaning.

---

## 6. Cell-line specificity scores (`K562_specificity_score`, `GM12878_specificity_score`, `HepG2_specificity_score`)

**Source.** `region_concept_axes.parquet`. One column per viz cell line.

**Definition (per cell line `C`, others = the two non-`C` cell lines):**

```
this_max(t) = max over {H3K27ac, H3K4me1, H3K4me3} of `C__{mark}__mean` for region t
others_max(t) = max over {H3K27ac, H3K4me1, H3K4me3} for the *other* cell lines (max-then-max)

case_set    = tokens where this_max ≥ percentile(75, this_max) AND others_max < median(others_max)
control_set = tokens where others_max ≥ percentile(75, others_max) AND this_max < median(this_max)

specificity_axis(C) = mean(R2V[t] for t in case_set) − mean(R2V[t] for t in control_set)
specificity_score_C(t) = cos(R2V[t], specificity_axis(C))
```

Case captures "specifically high in C, not high in others"; control captures "high in others, not in C." The difference of their embedding centroids is the lineage-specificity axis.

**Range.** Bounded [-1, +1]. Empirical range ~[-0.3, +0.3] (smaller than anchor / activity — most regulatory elements aren't sharply lineage-specific in their R2V signature).

**Color.** PuOr divergent. Domain `[-0.3, +0.3]`. Orange = `C`-specific; purple = active in non-`C` contexts.

**Empirical case/control sizes (chr16):**

| Cell line | case_set | control_set |
|---|---|---|
| K562 | 2,110 | 2,289 |
| GM12878 | 2,269 | 2,648 |
| HepG2 | 2,242 | 2,120 |

**Interpretation.** Even though the absolute magnitudes are small, **spatial coherence in the UMAP carries the signal**. The dELS class (yellow blob under SCREEN coloring) splits into different geographic islands depending on which lineage's specificity axis you project on:

- **K562**: sharp, localized orange islands (lower-left + right edge). Erythroid-specific dELS regions cluster tightly.
- **GM12878**: more diffuse spread across the upper wing; different regions than K562.
- **HepG2**: broadest spread, less focal than K562.

The K562 pattern being sharpest is partly a corpus-coverage artifact: ~3,200 K562 files vs ~900 GM12878 vs ~1,800 HepG2. More training signal for K562 = sharper learned lineage representation.

**The dictionary's argument made visible.** Three lineages, three separable directions in the embedding, even though no lineage label appeared in R2V's training. SCREEN's class label says "this is dELS"; the concept-axis projection says "this is a K562-specific dELS" — sub-class structure invisible to the categorical view.

---

## 7. Target evidence count (log) (`n_evidence_rows`)

**Source.** `region_target_evidence_summary.parquet` (from stage 06), column `n_evidence_rows`. Augmented inline in the viz with `+1` to permit log scaling (tokens with no evidence map to log(1) = 0, the lowest color).

**Definition.** Per chr16 token: count of `(region, target_gene, evidence_type, biosample/tissue)` rows in `region_target_evidence.parquet`. Stage 06 sources:
- **3D Chromatin contacts**: Hi-C / ChIA-PET / HiChIP / pCHi-C entries from ENCODE V4 cCRE-Gene Links
- **eQTL associations**: GTEx + similar variant→gene links

CRISPR perturbation entries are absent for chr16 (CRISPRi tile-screens have only covered chr3/10/12/19/X so far in V4).

**Range.** Empirical: 1 (no evidence) to ~200+ for hub regions. Median ~56 across the 86% of tokens that have any evidence.

**Color.** Viridis log scale. Bright = many evidence rows; dark = few.

**Interpretation.** A "well-studiedness" / "hub-like" indicator. Bright regions are often canonical regulatory elements with many measured 3D contacts and eQTL associations across many biosamples. Dark regions are under-studied or genuinely solitary in the corpus. This is a coverage/study-bias signal, not strictly a biology signal, but useful for navigating the dictionary toward regions with rich evidence.

**Caveat — what `target_gene` means.** For PLS focal regions, the "target genes" in the underlying evidence are mostly other genes whose promoters this PLS contacts in 3D space (TAD-mate genes), NOT genes regulated by this region. For dELS / pELS focal regions, target_gene is more often the regulated promoter. The count metric collapses both readings.

---

## Implementation references

- **Concept axes**: `genomic-dict/pipeline/scripts/embedding_features.py`. Standalone non-stage script; loads `databio/r2v-encode-hg38` from HF cache, slices to chr16's 35,934 tokens via `viz_chr16.token_id`, computes prototype centroids and direction vectors, projects, writes `region_class_prototypes.parquet` + `region_concept_axes.parquet`.
- **Target evidence**: `genomic-dict/pipeline/scripts/06_target_evidence.py`. Maps viz tokens to V4 EH38E accessions via coordinate match (94.5% exact + overlap fallback to 99.2% total), filters chr16 ENCODE V4 cCRE-Gene Links to mapped tokens, writes `region_target_evidence.parquet` + `_summary.parquet`.
- **Bigwig means** (used by activity / specificity axes): `genomic-dict/pipeline/scripts/04_extract_intrinsic.py`. Sources from per-(cell_line × mark) bigwigs at chr16 cCRE coordinates.
- **Color directives in viz**: `genomic-regions/src/index.md` Section 2, in the `regionUmap` cell — see the `colorChannel` and `colorDirectives` selectors for the full mapping.

## Future extensions

Color encodings to consider as follow-ups:
- **PMI degree** — count of strong PPMI partners in `corpus_baseline` (replaces `n_files_total` once stage 12 swap happens).
- **Module size** — size of the Leiden module the focal token belongs to in a chosen stratum. Surfaces "this region is in a tight module" vs "in a giant catchall."
- **NPMI dispersion** — variance of NPMI scores across strata for this token. High dispersion = context-dependent grammar; low = robust grammar.
- **Stratum-specific marginal** — once the stratum lens picker is built, color regions by their `n_files_active / n_files_in_stratum` ratio in the selected lens.
- **Per-cclass softmax** — encode the dominant cclass *centroid distance* directly (closest-class soft profile), bypassing the hard `cclass` label.
