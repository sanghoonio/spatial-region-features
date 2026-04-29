---
date: 2026-04-28
status: in-progress
description: Pipeline + viz plan for "A Grammar of Regulatory Genomics" — replacing the dropped stage 12/13 with curated-strata PMI cooccurrence and Leiden community detection, and rebuilding the Step 5 / per-region experience around a Fudoki-inspired dictionary-entry UX. Self-contained design doc for a fresh session.
---

# A Grammar of Regulatory Genomics — pipeline + viz plan

This plan is the methodological backbone for the next iteration of the **genomic-regions** visualization (Observable Framework site at `/Users/sam/Documents/Work/genomic-regions`, deployed as "A Dictionary of Regulatory Genomics"). It explains why the existing per-region relationship surfaces aren't enough for the title's claim, what new corpus-derived structures we need, and how to build them.

The plan covers two sides:
- **Pipeline (this repo, `spatial-region-features/genomic-dict/`)** — new precomputed parquets that bring the corpus's grammar into legible form.
- **Viz (sibling repo, `genomic-regions/`)** — UX changes that turn the per-region experience into something that earns the word "dictionary."

The viz changes are framed at the design level here; full code-level execution belongs to a separate session in that repo, after the data layer lands.

---

## TL;DR

1. **The title commits us to a real dictionary.** Not a UMAP with tooltips — a place where every regulatory region has a substantive entry: class, chromatin signature, multi-context senses, grammatical relations, usage citations. The current viz has a strong index (UMAPs) and a weak entry (deleted Step 5 cooccurrence panel). Fudoki is the UX reference for what a strong entry looks like.
2. **The previous Step 5 cooccurrence approach failed for two diagnosed reasons:**
   - **Popularity bias** — broadly-active regions co-occur with everything. Raw counts and Jaccard reward this.
   - **Modality aggregation** — pooling cooc across DNase + ATAC + ChIP + every histone mark mixes activate/repress/insulator signals into a single bag.
3. **The fix is two methodological upgrades to corpus cooccurrence**, not abandonment of it:
   - **PMI** (pointwise mutual information) instead of Jaccard or raw counts. Discounts popularity directly.
   - **Multi-level curated strata** (L1 corpus → L2 cell-line → L3 assay → L4 biology pan-cell → L5 cell × biology). Each stratum is a named biological question at a stated level of zoom. The dictionary entry reads a region across multiple levels, surfacing where they agree or diverge — that's the grammar doing real work. Broader strata aren't broken; they answer broader questions.
4. **A "regulatory sentence" is the new unit of grammar** — a small, ordered set of regions whose ordering encodes biological role (anchor → modifiers → boundary), built from corpus statistics, not from genomic position. Three formulations (A, B, C) of varying ambition; we pursue all three. Sentence renderings include both Fudoki-style colored token streams *and* UMAP focal-lens overlays (the embedding's geometric view of the same sentence).
5. **Pipeline footprint**: new stage 12 (curated-strata PMI cooccurrence), new stage 13 (Leiden modules), unchanged stage 14 (featured signal). Plus inline derived-from-embedding features — prototype distances, concept axes — that ride alongside without new pipeline stages. Manifest `target` is partially populated (verified: 100% TF ChIP, 91% Histone ChIP, ~0% generic ChIP); curated strata that need histone marks may want filename-heuristic backfill, but not as a blocker.
6. **Proof-of-concept**: the α-globin cluster at chr16:218k–238k must surface as a coherent module with HBA1 promoter as anchor and HS-40 enhancer as a top partner. If it doesn't, the corpus or the strata are inadequate and the framing needs to change. Plus a cheap **R2V faithfulness pre-flight** (compare embedding kNN to corpus-wide PPMI partners) as a sanity check before locking strata.

---

## What's already happened in this session (2026-04-28)

These deletes have shipped before this plan was written. Captured here so a fresh session sees the starting state:

- **Old stage 13 (`13_cooccurrence.py`)** — produced `region_cooccurrence.parquet` with top-30 partners per (token, context) where context ∈ {`all`, `K562`, `GM12878`, `HepG2`} and metric was Jaccard. **Deleted.** Outputs felt random because of the popularity bias and modality aggregation issues above. Config block, results dir, and dead `cooccurrence_*` knobs in stage 08 also removed.
- **Old stage 12 (`12_corpus_stats.py`)** — produced `region_stats.parquet` with corpus-wide + per-cell-line + per-assay-class activation counts per chr16 token. **Deleted.** The viz used exactly **one** of its sixteen columns (`n_files_total`, as a UMAP color toggle). The new stage 12 (PMI cooccurrence) computes per-token marginals as a byproduct, so that one usage gets reproduced naturally.
- **Renumbering convention agreed:** new stages slot into the freed numbers. New stage 12 = PMI cooccurrence; new stage 13 = Leiden modules; stage 14 (featured signal, already done, drives Step 1's continuous/peaks/tokens toggle) stays put.
- **Genomic-regions viz still references `region_stats.parquet`** in five places in `src/index.md`. The bundle copy at `genomic-regions/src/data/dictionary/region_stats.parquet` is now an orphan. **Decision (see Locked-in decisions, #1): keep the orphan in place** until new stage 12 lands; the viz keeps reading the stale file rather than break the "File count (log)" UMAP color toggle.

---

## Why the title matters: dictionary as design constraint

The viz is titled "A Dictionary of Regulatory Genomics." This isn't decoration — it's a methodological claim. A dictionary, taken seriously as a UX, demands:

| Dictionary affordance | What it means for regulatory regions |
|---|---|
| **Lookup** any word | Find any region by location, class, or similarity |
| **Read** the entry | A substantive description of what kind of element this is |
| **See the senses** | The same region behaves differently in different cell lines / contexts |
| **Find related entries** | Grammatical relations — what does this region pattern with |
| **Read usage examples** | Which BED files / experiments produce this region |
| **Browse** | Navigate from one entry to neighbors and keep going |

Translated to a Fudoki-style UX (https://github.com/iamcheyan/fudoki, the working Japanese-text-analysis app we're using as a UX reference): the document is the canvas; every token is colored by class; every token is clickable to open a persistent dictionary card; the card is multi-section and substantive.

**The current genomic-regions viz** has a strong *index* (the file UMAP, the region UMAP colored by SCREEN class) but a weak *entry*. Step 5 was supposed to be the entry; the deleted cooccurrence panel didn't earn that role.

**The work in this plan is to build the entry.** Specifically: turn each region into a navigable dictionary entry whose grammatical-relations section is corpus-grounded, biology-coherent, and substantive enough to actually look up.

---

## Failure analysis: why the previous Step 5 didn't work

The deleted stage 13 produced top-30 partners per token using Jaccard, with contexts `[all, K562, GM12878, HepG2]`. The output looked plausible at column-schema level and had ~143k rows of "partner lists." But on inspection the partners felt random — like reading a thesaurus that says "near 'house': building, structure, dwelling, edifice, residence, …" — technically correct, lazy.

Two diagnosed mechanisms:

### Mechanism 1: popularity bias in raw cooccurrence

A region that's broadly accessible (e.g., a generic open-chromatin spot) gets peak-called in many BED files for trivial reasons. Its top-Jaccard partners end up being other broadly-accessible regions — co-occurrence dominated by file-level co-membership noise rather than biological association.

PMI (pointwise mutual information) directly discounts popularity:

$$\text{PMI}(a, b) = \log \frac{P(a, b)}{P(a) \cdot P(b)}$$

Read: "how much more often than chance do a and b appear together?"
- `> 0` → associated more than chance (biological signal)
- `0` → independent (just popular)
- `< 0` → anti-associated (often filtered out as PPMI = max(PMI, 0))

**Concrete toy example.** Imagine 1000 NYC-text documents, with:
- "the": 950 docs (P = 0.95)
- "Brooklyn": 50 docs (P = 0.05)
- "Bridge": 30 docs (P = 0.03)
- (the, Brooklyn) cooc = 48 (almost every Brooklyn doc has "the")
- (Bridge, Brooklyn) cooc = 25

| Metric | (the, Brooklyn) | (Bridge, Brooklyn) | Verdict |
|---|---|---|---|
| Raw count | 48 | 25 | "the" looks more associated. **Wrong.** |
| Jaccard | 48/952 = 0.05 | 25/55 = 0.45 | Better, but still rewards moderately-popular partners |
| PMI | log(0.048 / (0.95·0.05)) ≈ 0.01 | log(0.025 / (0.03·0.05)) ≈ 2.81 | Correctly identifies "Bridge" as the real partner |

PMI is the standard NLP fix for the same problem. Word2Vec's negative-sampling objective is mathematically related to factorizing the PPMI matrix — so R2V (Word2Vec-derived, trained on bedbase) implicitly learned a smoothed version of PPMI. The kNN partners in `viz_chr16.knn_*` are the model's compressed/smoothed view; raw PPMI on the corpus is the unsmoothed evidence.

### Mechanism 2: modality aggregation — and the limits of how big a deal this is

Even after stratifying by cell line (the "K562" / "GM12878" / "HepG2" contexts), the previous approach pooled together every assay and every mark within each cell line. "K562" cooccurrence mixed:
- DNase peaks (open chromatin)
- ATAC peaks (open chromatin, similar to DNase)
- TF ChIP peaks (specific factor binding)
- H3K27ac peaks (active enhancer)
- H3K27me3 peaks (Polycomb-repressed)
- H3K9me3 peaks (heterochromatin)
- H3K4me1, H3K4me3, etc.

The naive worry is that mixing H3K27ac (active) with H3K27me3 (repressive) within "K562" produces incoherent partner lists, since these marks have opposite biological meaning. **In practice this concern is overstated.** A region active in K562 H3K27ac files and a region repressed in K562 H3K27me3 files mostly *don't appear in the same files at all*. Active-context files and repressed-context files have largely disjoint region complements, so cooccurrence within K562 doesn't mix them as much as a flat description would suggest. Each region's partner list within K562 tends to be naturally coherent — its file co-memberships reflect either its active or its repressed context, not both.

**So K562-aggregate cooccurrence isn't broken; it's just answering a different (broader) question.** It tells you "what does K562's overall regulatory landscape look like for this region?" rather than "what's this region's active-enhancer grammar?" Both are valid; the second is sharper and more interpretable for a dictionary entry.

**The real fix is multi-scale curated strata.** The strata library is research-design that names biological questions at different levels of zoom. The dictionary entry's "senses" section reads a region's grammar at multiple levels and shows where the answers agree vs differ.

The library spans levels (full level-by-level table appears in stage 12 below):

- **Corpus level**: `corpus_baseline` — pan-everything, broadest reference.
- **Cell-line broad**: `featured_lineage_K562`, `featured_lineage_GM12878`, `featured_lineage_HepG2` — "what's this cell line's overall regulatory landscape." The strata I previously called "broken." They're not broken; they answer a coarse but real question.
- **Assay broad**: `open_chromatin_pan_cell`, `tf_bound_pan_cell` — "what pairs in accessibility / TF-binding contexts pan-cell."
- **Biology pan-cell**: `active_enhancers_pan_cell`, `active_promoters_pan_cell`, `polycomb_repressed_pan_cell`, `heterochromatin_pan_cell`, `ctcf_boundaries` — "what's the grammar of this specific biological state across cell types."
- **Cell × biology**: `erythroid_active`, `lymphoid_active`, `hepatic_active` — "what's this lineage's specific-context grammar."

The principle: **each stratum is a named biological question at a stated level of zoom**, and the dictionary surfaces multiple levels for the same region. Different levels answer different questions; agreement across levels signals robustness; divergence across levels reveals context-dependence.

### Combined fix

The new stage 12 will compute **PMI** (popularity-discounted metric) within **multi-level curated strata** (each stratum a named biological question at a stated level of zoom). PMI handles popularity bias; the level-by-level strata structure makes the dictionary's grammar legible at multiple scales. Each level answers a different question; the entry exposes them side-by-side.

---

## The new framing: regulatory sentences

A "sentence" here is a small, ordered set of regions whose ordering encodes grammatical role (anchor → modifiers → boundary), constructed from corpus distributional evidence, **not from genomic position**. Genomic position is not the right ordering principle for regulatory grammar — promoter-enhancer pairs can be megabases apart, and reading regions in chromosome order is the boring genome-browser view this project is trying to escape.

For the α-globin cluster (chr16:218k–238k), the sentence we'd hope to construct from corpus statistics:

| Role | Region | Class | Why it's there |
|---|---|---|---|
| Head | HBA1 promoter | PLS | The thing the sentence is about |
| Co-head | HBA2 promoter | PLS | Parallel duplicated gene |
| Primary modifier | HS-40 enhancer | dELS | Canonical α-globin master enhancer |
| Secondary modifiers | HS-33, HS-48 | dELS | Auxiliary enhancers in MCS-R |
| Pre-modifier | HBQ1 promoter | PLS | Bidirectional partner |
| Boundaries | flanking CTCF sites | CA-CTCF | Domain delimiters |

The genomic order would put HBQ1 at one end and CTCF sites bracketing — a linear strip. The grammatical order is hierarchical: HBA1 anchor → HS-40 (strongest enhancer) → HS-33, HS-48 (other enhancers) → CTCF (boundary) → HBQ1 (parallel head). The grammatical reading reveals structure that linear reading hides.

We pursue **three formulations** of "sentence," each appropriate at different scales:

### A — Anchor-based sentences (per-region, on demand)

User clicks any region (the **anchor**); the system builds a sentence:
1. Pull anchor's top-K corpus partners from new stage 12's PMI cooccurrence (the relevant stratum, see below).
2. Order partners by SCREEN class hierarchy: `PLS → pELS → dELS → CA-H3K4me3 → CA-CTCF → unclassed`.
3. Within class, order by PMI score (strongest first).

Stratum selection: by default, pick the stratum where the anchor is most active (e.g., if HBA1 is most active in `active_promoters_pan_cell` BEDs, default to that stratum). Allow the user to switch lenses.

**Two parallel rendering modes for the same sentence:**

- **Fudoki stream** — anchor + partners as a horizontal flow of colored class-tagged tokens, in grammatical order. The "linguistic" view. Each token clickable → re-anchor → new sentence.
- **UMAP focal lens** — same anchor + partners highlighted on the existing region UMAP, with other tokens dimmed. Optional convex-hull overlay of the partner set. The "geographic" view of where this sentence lives in semantic space. Reveals shape information the linear stream can't (compact patch = locally coherent neighborhood; spread/split patch = semantically diverse partner set).

A — exists in two forms, based on **two evidence sources**:

- **A-PMI**: partners from stage 12's PPMI cooccurrence (direct corpus evidence; popularity-discounted).
- **A-kNN**: partners from `viz_chr16.knn_*` (R2V's smoothed similarity; already in bundle).

Both are corpus-grounded — R2V was trained on cooccurrence, so kNN is a smoothed view of the same evidence PMI exposes directly. Showing them side-by-side in the dictionary entry is a methodological feature: agreement → robust grammar; disagreement → diagnostic (model-smoothed transitive similarity vs raw cooccurrence evidence).

### B — Module-based sentences (precomputed, browsable)

Pre-cluster the cooccurrence graph into communities (Leiden algorithm). Each community = a sentence body, conceptually a "regulatory module."

For each module:
- Identify an **anchor** (highest within-module centrality, or the dominant PLS if any).
- Identify **role-tagged members** by class.
- Identify a **dominant signature** (which marks / cell lines characterize this module).
- Generate an **auto-label** ("H3K27ac-active K562 module containing HBA1 / HBA2 / HS-40").

Modules can be computed per stratum, giving different "languages" of grammar. The `active_enhancers_pan_cell` stratum produces one set of modules (active enhancer-promoter clusters); `ctcf_boundaries` produces another (loop boundaries); `polycomb_repressed_pan_cell` produces another (Polycomb domains). The viz lets the user switch which language they're viewing.

This is the **catalogue** of the dictionary — the alphabetical-by-topic browse view. Clicking a module opens its sentence; clicking a member of a sentence is option A.

### C — Path-traversal sentences (region-to-region paths)

Given two regions A and B, find a path connecting them through the corpus-grammar graph — a chain of regions where each consecutive pair has an explicit similarity link.

**Critical interpretation note:** A path is a **chain of pairwise similarities**, not a biological mechanism. C is a *guided-tour / browsing* feature, not a claim about mechanism, temporal sequence, or causal regulation. Reading "HBA1 → HBA2 → HS-40 → CTCF" tells you each consecutive pair shares strong PPMI cooccurrence; it does NOT claim HBA1 regulates HBA2 regulates HS-40. UI copy must convey this honestly.

**Where paths are computed:**

- **Primary: PMI graph paths.** Edges from stage 12 with PPMI above threshold; weighted shortest path with `weight = 1 / (PPMI + ε)`. Each edge has direct corpus evidence (these regions actually appear together in BED files more than chance). Most biologically grounded.
- **Fallback: 100-dim embedding space.** Greedy nearest-neighbor walk or A* with cosine distance. Useful when the PMI graph is sparse or the endpoints are in disconnected components. Each edge is "R2V thinks these regions are similar" — smoother but more interpretively distant.
- **Never: UMAP 2D coordinates.** UMAP distorts global distances; paths computed on UMAP layout are geometric noise, not signal. UMAP serves as a **canvas** for visualization only.

**Rendering:** path nodes overlaid on the region UMAP, drawn as a connected polyline through the 2D positions of each path node. The path was computed in PMI/embedding space; the UMAP gives it a spatial reading. Each edge labeled with its similarity score (PPMI value or cosine distance). Path nodes also rendered as a Fudoki stream below the UMAP overlay — same chain, two visual modes.

**What paths can reveal:**

- **Bridge regions.** Nodes that paths between distant endpoints often pass through are bridges between regulatory contexts — biologically interesting because they connect otherwise-separate parts of the landscape.
- **Modular structure.** Comparing paths across multiple endpoint pairs surfaces shared intermediates → candidate "hubs."
- **Lens-dependent connectivity.** Paths in the `active_promoters_pan_cell` stratum vs the `polycomb_repressed_pan_cell` stratum can be drastically different, and that contrast is itself the dictionary's grammar at work.

**What paths do NOT reveal:** causal cascades, temporal sequences, genomic adjacency, regulatory mechanism. Honest framing in copy is load-bearing.

This is the most exploratory of the three formulations. UX bottleneck is real (how to let the user pick endpoints, when to surface paths automatically, etc.). Build A and B first; layer C on top.

### Why all three are corpus-grounded (not the 18-column trap)

An earlier idea in this conversation was to use the 18 `{cell_line}__{mark}__mean` signal columns in `viz_chr16` as a region similarity surface (cosine on the 18-dim vector). **That was rejected** because it abandons the project's distributional-semantics thesis: it substitutes a thin handpicked feature panel for the corpus evidence. Worse, it would mostly re-discover SCREEN class labels (which are already derived from those marks).

A, B, and C all stay inside the corpus framing. The 18-dim signal columns can show up as part of the dictionary entry's *signature* section — describing a region's chromatin state — but they don't drive the relations, because relations are corpus-distributional. This is the load-bearing methodological commitment.

---

## Methodology: PMI + Leiden community detection

### PMI computation (stage 12)

For each `(token_a, token_b, stratum)`:

$$P(a) = \frac{|\text{files in stratum where a is active}|}{|\text{files in stratum}|}$$

$$P(a, b) = \frac{|\text{files in stratum where both a and b are active}|}{|\text{files in stratum}|}$$

$$\text{PMI}(a, b | \text{stratum}) = \log \frac{P(a, b)}{P(a) \cdot P(b)}$$

$$\text{PPMI} = \max(0, \text{PMI})$$

We store PPMI; the negative tail is dominated by noise from low-count pairs and rarely useful in NLP applications.

For each (token, stratum), keep top-K partners by PPMI (K=30). Also store the raw count `c(a, b)` and Jaccard for backwards compatibility.

**Statistical floor.** Tokens with `n_files_active < 5` in a stratum are excluded — too little evidence to estimate `P(a)` reliably. Sweep this threshold during validation.

### Leiden modules (stage 13)

Build a graph from stage 12:
- Nodes: chr16 tokens
- Edges: pairs with PPMI above threshold (e.g., > 1.0; sweep during validation)
- Edge weights: PPMI

Run Leiden (Traag et al. 2019) at multiple resolutions γ ∈ {0.5, 1.0, 1.5, 2.0}. Output a partition: each token assigned a `module_id` per (stratum, γ) combination.

Leiden optimizes modularity:

$$Q = \frac{1}{2m} \sum_{ij} \left[ A_{ij} - \frac{k_i k_j}{2m} \right] \delta(c_i, c_j)$$

where `A_ij` = edge weight, `k_i` = node degree, `m` = total edges. Maximizing Q finds groupings denser than random expectation given the degree sequence.

For each module: compute within-module **eigenvector centrality**; the highest-centrality node is the **anchor**. Label modules by their dominant SCREEN class and dominant active cell lines (auto-label until stage 06 / GENCODE proximity is wired up).

### Worked NLP example for intuition

Eight NYC-related tokens with PPMI edges (illustrative numbers):

```
Brooklyn ─── 2.8 ─── Bridge ─── 2.1 ─── Heights
   │                                       │
   2.5                                    1.9
   │                                       │
  Bronx ─── 2.7 ─── Yankees ─── 3.1 ─── Stadium

Manhattan ─── 2.4 ─── Empire ─── 3.0 ─── State ─── 2.8 ─── Building
```

Faint cross-cluster edges (PPMI < threshold) are dropped. Leiden output:
- Community A: {Brooklyn, Bridge, Heights}
- Community B: {Manhattan, Empire, State, Building}
- Community C: {Bronx, Yankees, Stadium}

Three "topics" emerge from cooccurrence statistics alone. No labels supplied. This is the same pattern PPMI + Leiden produces on Wikipedia text — topic communities without supervision. Translated to the regulatory case: regulatory modules without supervised tissue/element labels.

### Sanity-check dataset (optional)

Zachary's Karate Club graph (`networkx.karate_club_graph()`) — 34 nodes, ground-truth two-faction split. Leiden recovers it ~85–90% of runs. Useful for sanity-checking the algorithm and modularity sweep before pointing it at noisy biology data.

---

## Pipeline plan

### Repository state after this plan

```
genomic-dict/pipeline/scripts/
├── 00_inspect_metadata.py       (unchanged)
├── 01_curate_corpus.py          (`target` verified populated for TF ChIP / Histone ChIP; 4000 generic ChIPs untagged — accepted gap)
├── 02_prepare_universe.py       (unchanged)
├── 03_fetch_bigwigs.py          (unchanged)
├── 04_extract_intrinsic.py      (unchanged)
├── 05_extract_sequence.py       (unchanged; not run; gap)
├── 06_extract_extrinsic.py      (unchanged; not run; gap; recommended for module labels)
├── 07_load_pretrained.py        (unchanged)
├── 08_precompute_viz.py         (unchanged)
├── 09_prepare_file_viz.py       (unchanged)
├── 10_featured_narrative.py     (unchanged)
├── 11_tokenize_corpus_chr16.py  (unchanged)
├── 12_cooccurrence_pmi.py       NEW
├── 13_modules.py                NEW
└── 14_featured_signal.py        (unchanged)
```

### Manifest target coverage (verified 2026-04-28)

The corpus manifest at `data/corpus/manifest.parquet` (16,799 rows, 21 columns) has a `target` column. Coverage by assay:

| Assay | Target populated | Notes |
|---|---|---|
| TF ChIP-seq | 4000/4000 (100%) | every file has its TF (MYC, GATA1, etc.) |
| Histone ChIP-seq | 724/799 (91%) | most have the mark (H3K4me3=158, H3K27me3=76, H3K27ac=64, H3K4me1=60, H3K9me3 small) |
| ChIP-seq (generic) | 1/4000 (0.03%) | almost certainly mostly histone ChIP that wasn't categorized tightly by bedbase |
| ATAC-seq | 0/4000 (0%) | expected — no target concept |
| DNase-seq | 0/4000 (0%) | expected — no target concept |

The 4000 generic "ChIP-seq" files are the gap. Decision: **proceed without filename-heuristic backfill for v1**; if α-globin validation needs more files in the histone strata, add backfill as a step.

No re-curation of stage 01. Corpus balance preserved.

### Stage 12: curated-strata PMI cooccurrence

**Inputs:**
- `data/corpus/manifest.parquet` (with `target` for ChIP files)
- `data/annotations/tokenized_corpus_chr16.parquet`

**Strata as named biological questions, at five levels of zoom.** Each stratum is a curated filter rule answering a specific question. The library is **multi-level**: broader strata describe overall landscapes, narrower strata isolate specific biological contexts. The dictionary entry's "senses" section reads a region across multiple levels, surfacing where they agree vs diverge. Sizes are estimates; `target`-dependent strata may grow if filename-heuristic backfill is added later.

| Level | Stratum | Biological question | Filter | Est. files |
|---|---|---|---|---|
| **L1** corpus | `corpus_baseline` | What pairs corpus-wide (broadest reference)? | (none) | 16,799 |
| **L2** cell-line broad | `featured_lineage_K562` | What's K562's overall regulatory landscape? | cell_line == K562 | (full subset) |
| **L2** cell-line broad | `featured_lineage_GM12878` | … GM12878's? | cell_line == GM12878 | (full subset) |
| **L2** cell-line broad | `featured_lineage_HepG2` | … HepG2's? | cell_line == HepG2 | (full subset) |
| **L3** assay broad | `open_chromatin_pan_cell` | What pairs in accessible chromatin? | assay ∈ {DNase-seq, ATAC-seq} | 8000 |
| **L3** assay broad | `tf_bound_pan_cell` | What pairs in any TF-bound region? | assay == TF ChIP-seq | 4000 |
| **L4** biology pan-cell | `active_enhancers_pan_cell` | What pairs in active-enhancer chromatin? | target ∈ {H3K27ac, H3K4me1} | ~120 |
| **L4** biology pan-cell | `active_promoters_pan_cell` | What pairs in active-promoter chromatin? | target == H3K4me3 | ~158 |
| **L4** biology pan-cell | `polycomb_repressed_pan_cell` | What pairs in Polycomb-repressed regions? | target == H3K27me3 | ~76 |
| **L4** biology pan-cell | `heterochromatin_pan_cell` | What pairs in constitutive heterochromatin? | target == H3K9me3 | small |
| **L4** biology pan-cell | `ctcf_boundaries` | What pairs in CTCF-bound elements? | target == CTCF | ~111 |
| **L5** cell × biology | `erythroid_active` | What's the grammar of erythroid-active regions? | cell_line == K562 AND target ∈ {H3K27ac, H3K4me1, H3K4me3} | ~25–40 |
| **L5** cell × biology | `lymphoid_active` | … of lymphoid-active regions? | cell_line == GM12878 AND target ∈ active marks | ~25–40 |
| **L5** cell × biology | `hepatic_active` | … of hepatic-active regions? | cell_line == HepG2 AND target ∈ active marks | ~25–40 |

= **14 starter strata across 5 levels.** Per-(mark × cell-line) cross-cuts beyond active-marks (e.g., `K562_H3K27me3` specifically) and per-TF strata (e.g., `GATA1_bound`) are deferred. The strata library is research-design that you (and a future paper) own — these are starter values; the final list can grow or shrink based on what's biologically informative.

**Why levels matter for the dictionary entry.** A reader looking up a region wants to read its grammar at multiple zooms — the L2 K562 grammar tells one story; the L4 active-enhancer grammar another; the L5 erythroid-active grammar a third. Where they agree, the grammar is robust across zoom. Where they diverge, the entry has surfaced context-dependent grammar — this region behaves differently in different biological frames. That's the dictionary doing real work.

**Output: `data/precomputed/region_cooccurrence_pmi.parquet`**

| col | type | meaning |
|---|---|---|
| `token_id` | Int64 | focal region |
| `stratum` | String | named curated stratum (e.g., `active_enhancers_pan_cell`) |
| `n_files_in_stratum` | Int64 | denominator for P(a) — files matching this stratum's filter |
| `n_files_active` | Int64 | files in stratum where token is active |
| `partner_token_ids` | List[Int64] | top-K partners (K=30) |
| `weights_ppmi` | List[Float32] | PPMI score, parallel |
| `weights_jaccard` | List[Float32] | Jaccard, parallel |
| `counts` | List[Int64] | raw cooccurrence count, parallel |

Approx size: 35,934 tokens × 14 strata × 30 partners ≈ 15M rows max; expect 25–50 MB after parquet compression. (Rows where `n_files_active < 5` are dropped from output. Rare-mark strata like `heterochromatin_pan_cell` will contribute few rows.)

**Per-token marginals as a byproduct.** Computing PMI requires `n_files_active` per (token, stratum), which is exactly what the deleted stage 12 produced. We don't need a separate `region_stats` parquet anymore — those numbers live in this output, queryable via `WHERE token_id = ? AND stratum = 'all'`.

**Algorithm sketch:**
1. Inner-join manifest with tokenized corpus to get `(file_id, cell_line, target, [active token ids])`.
2. For each stratum:
   - Filter to files matching stratum criteria.
   - Build sparse `(file × token)` activation matrix `X` (CSR).
   - Compute marginals: `n_a = X.sum(axis=0)` (per-token file counts), `N = X.shape[0]` (file count in stratum).
   - Compute joint counts: `C = X.T @ X` (sparse matmul; diagonal is `n_a`; off-diagonal is `c(a,b)`).
   - Compute PPMI per nonzero off-diagonal entry: `PPMI(a,b) = max(0, log(c(a,b) * N / (n_a * n_b)))`.
   - For each row, keep top-K by PPMI (with raw count > 0).
3. Concat all strata, write parquet.

Compute estimate: each stratum is one sparse matmul. The `all` stratum is 16,799 × 35,934 ~= 600M cells max, but very sparse (~3000 active per row on average) → ~50M nonzero. Sparse matmul in scipy is feasible; ~10–30 seconds per stratum. Total ~5–10 minutes. Memory: peak ~few GB if naive; with chunking, well under that. Runs locally; no SLURM needed.

**Config block (replaces deleted stage 13 block):**

```yaml
12_cooccurrence_pmi:
  top_k_partners: 30
  ppmi_threshold: 0.0          # drop pairs with PPMI <= this; 0 = keep all positive
  min_files_active: 5          # statistical floor per (token, stratum)
  strata:
    corpus_baseline:
      filters: {}
      description: "Corpus-wide baseline; no filter"
    active_enhancers_pan_cell:
      filters:
        target: ["H3K27ac", "H3K4me1"]
      description: "Active enhancer signatures across cell types"
    active_promoters_pan_cell:
      filters:
        target: ["H3K4me3"]
      description: "Active promoter mark across cell types"
    polycomb_repressed_pan_cell:
      filters:
        target: ["H3K27me3"]
      description: "Polycomb-repressed regions across cell types"
    heterochromatin_pan_cell:
      filters:
        target: ["H3K9me3"]
      description: "Constitutive heterochromatin across cell types"
    ctcf_boundaries:
      filters:
        target: ["CTCF"]
      description: "CTCF-bound boundary elements"
    erythroid_active:
      filters:
        cell_line: ["K562"]
        target: ["H3K27ac", "H3K4me1", "H3K4me3"]
      description: "Active chromatin in erythroid-lineage cells"
    lymphoid_active:
      filters:
        cell_line: ["GM12878"]
        target: ["H3K27ac", "H3K4me1", "H3K4me3"]
      description: "Active chromatin in lymphoblastoid cells"
    hepatic_active:
      filters:
        cell_line: ["HepG2"]
        target: ["H3K27ac", "H3K4me1", "H3K4me3"]
      description: "Active chromatin in hepatic cells"
    open_chromatin_pan_cell:
      filters:
        assay: ["DNase-seq", "ATAC-seq"]
      description: "Accessible chromatin across cell types"
    tf_bound_pan_cell:
      filters:
        assay: ["TF ChIP-seq"]
      description: "TF-bound regions across all surveyed TFs"
    featured_lineage_K562:
      filters:
        cell_line: ["K562"]
      description: "L2 — Full K562 slice (overall regulatory landscape, all assays)"
    featured_lineage_GM12878:
      filters:
        cell_line: ["GM12878"]
      description: "L2 — Full GM12878 slice (overall regulatory landscape, all assays)"
    featured_lineage_HepG2:
      filters:
        cell_line: ["HepG2"]
      description: "L2 — Full HepG2 slice (overall regulatory landscape, all assays)"
```

Description fields prefixed with their level (L1–L5) for legibility.

**Filter rule semantics.** Each filter key is a manifest column; each value is a list of acceptable values (logical AND across keys, OR within a key). An empty `filters: {}` means "no filter" (full corpus). Adding a stratum is a one-line config edit.

### Stage 13: regulatory modules (Leiden)

**Inputs:** `data/precomputed/region_cooccurrence_pmi.parquet`

**Approach:**
1. For each stratum (and each γ value), build an undirected graph: nodes = tokens with at least one PPMI > threshold edge in this stratum; edges = symmetrized top-K partner pairs with PPMI as weight.
2. Run Leiden via `igraph` or `leidenalg`. Set `partition_type=ModularityVertexPartition` (or RBConfigurationVertexPartition with γ).
3. For each `(stratum, γ, module_id)`:
   - Compute eigenvector centrality on the within-module subgraph; flag the highest as `is_anchor`.
   - Compute dominant SCREEN class, dominant active cell lines, dominant target marks (from manifest of files where module members co-activate).

**Output: `data/precomputed/region_modules.parquet`**

| col | type | meaning |
|---|---|---|
| `token_id` | Int64 | foreign key to viz_chr16 |
| `stratum` | String | which language of grammar |
| `gamma` | Float32 | resolution parameter |
| `module_id` | Int64 | community assignment (per stratum, γ) |
| `within_module_centrality` | Float32 | eigenvector centrality |
| `is_anchor` | Bool | top-centrality node in module |

**Output: `data/precomputed/module_summary.parquet`**

| col | type | meaning |
|---|---|---|
| `stratum` | String | |
| `gamma` | Float32 | |
| `module_id` | Int64 | |
| `n_tokens` | Int64 | size |
| `anchor_token_id` | Int64 | |
| `dominant_class` | String | most common SCREEN class |
| `dominant_cell_lines` | List[String] | from member coactivity |
| `dominant_target_marks` | List[String] | from member coactivity |
| `auto_label` | String | e.g., "H3K27ac-active K562 module (anchor: HBA1 promoter)" |

**Compute estimate:** Leiden on a sparse graph with ~35k nodes and ~1M edges runs in seconds via `igraph`. For v1 (γ = 1.0 hardcoded per Locked-in decision #3): 14 strata × 1 γ = 14 runs, well under a minute total. The full γ sweep `[0.5, 1.0, 1.5, 2.0]` is available as a dev-side diagnostic if module quality is suspect; not user-facing.

**Config block:**

```yaml
13_modules:
  ppmi_threshold: 1.0          # edge filter for module graph
  resolutions: [1.0]           # v1 hardcode; expand to [0.5, 1.0, 1.5, 2.0] for offline diagnostic sweep
  centrality_method: "eigenvector"
  random_seed: 42
```

### Embedding-derived features (no new stages; inline derivation)

The R2V embeddings already in the bundle (`pretrained_universe.parquet` and `viz_chr16.knn_*`) can produce additional dictionary-entry surfaces *without* new pipeline stages. These are computed inline either as part of stage 12/13 or directly in the viz layer.

**Prototype distances (per-cclass).** For each SCREEN class, compute the centroid embedding (mean of token embeddings within that class). For every chr16 token, store its cosine distance to each centroid. Output: a small `region_class_prototypes.parquet` (35,934 rows × {distance per cclass}, ~few hundred KB). UI: a soft-classification bar in the dictionary entry — e.g., "84% PLS, 11% pELS, 3% dELS" — revealing sub-class structure that hard cclass labels hide.

**Concept axes (direction vectors in embedding space).** Define a small set of meaningful directions:
- `activity_axis` = mean(tokens with high K562 H3K27ac signal) − mean(tokens with high K562 H3K27me3 signal)
- `K562_specificity_axis` = mean(K562-active tokens) − mean(non-K562-active tokens)
- `GM12878_specificity_axis`, `HepG2_specificity_axis` — analogous
- `anchor_axis` = mean(PLS centroids) − mean(dELS centroids) ("promoter-like vs distal-enhancer-like")

Project every chr16 token onto each axis → interpretable scalar score per axis. Output: a small `region_concept_axes.parquet` (35,934 rows × {score per axis}, ~few hundred KB). UI: a radar chart or score table in the entry — "(activity: +0.7, K562-specific: +0.4, anchor: +0.6)."

This direction-vector technique is standard in NLP (gender / sentiment / formality axes from Word2Vec) and is genuinely candidate-novel for regulatory genomics. The axes are research-design artifacts; the starter set above is illustrative, the final list is something to commit to deliberately.

**Cluster validation against PMI modules.** Run k-means or HDBSCAN on the 100-dim chr16 token embeddings → embedding-derived clusters. Compare to stage 13's Leiden PMI modules. Where they agree → robust grammar. Where they disagree → diagnostic (R2V's smoothing pulled in transitively-similar, OR corpus signal R2V didn't compress). Output: agreement metric per module (Jaccard between embedding cluster and PMI module memberships) added to `module_summary.parquet`.

### Out-of-corpus extension (separate spike, not in this plan's main scope)

**Use case:** user pastes a non-coding GWAS hit at any hg38 location → R2V embeds the new region → look up nearest universe tokens → read their dictionary entries. Makes the dictionary actually usable for downstream variant interpretation rather than self-contained demo. This is the project's strongest external pitch.

**Implementation paths to investigate (1–2 day spike):**
- **Client-side R2V via ONNX or transformers.js.** Browser-only inference; preserves the static-site GitHub Pages model. Preferred if model size permits.
- **Server-side R2V endpoint.** Easier but breaks the static-site assumption.
- **Pre-compute genome-wide embedding cache.** R2V every 200bp window in hg38 → ~15M tokens. Bigger bundle but eliminates inference at query time.

Decision deferred. Worth a separate spike before scope-extending this plan.

### Optional supporting work

**Stage 06 (gene proximity).** Would let module summaries auto-label by nearest gene ("α-globin cluster module" instead of "module 47"). The script is written but not run. ~30 min compute. Locked-in decision: prompt user to run when label-quality matters; not blocking.

**Stage 05 (motifs / ChromHMM).** Would let dictionary entries say things like "this region is GATA1-bound" or "this region is in ChromHMM active-promoter state in K562." Bigger compute (FIMO motif scan ~1–2 hours). Locked-in decision: same as 06 — prompt to run when needed for entry richness.

**Featured-region curated sentences.** A small `featured_sentences.parquet` (~80 rows) with hand-curated grammatical-role tags for the 4 featured intervals. Drives the demo path deterministically without depending on whatever the Leiden run picks. Low priority; only worth doing after auto-output is validated.

### Bundle changes

The genomic-regions data bundle currently has 8 parquets (`viz_files`, `viz_chr16`, `region_stats`, `tokenized_corpus_chr16`, `featured_files`, `featured_intervals`, `featured_tracks`, `featured_signal`). After this plan:

| Parquet | Status |
|---|---|
| `viz_files` | unchanged |
| `viz_chr16` | unchanged |
| `region_stats` | **kept as orphan for now**; rewire viz to PMI marginals after stage 12 lands, then drop |
| `tokenized_corpus_chr16` | unchanged |
| `featured_files` | unchanged |
| `featured_intervals` | unchanged |
| `featured_tracks` | unchanged |
| `featured_signal` | unchanged |
| `region_cooccurrence_pmi` | **new (from stage 12)** |
| `region_modules` | **new (from stage 13)** |
| `module_summary` | **new (from stage 13)** |
| `region_class_prototypes` | **new (inline derivation from R2V embeddings)** |
| `region_concept_axes` | **new (inline derivation from R2V embeddings)** |

Net bundle size change: keep ~700 KB orphan (region_stats), add ~30–50 MB (PMI cooc) + ~5 MB (modules + summary) + ~1 MB (embedding-derived). New bundle ≈ 150–170 MB, up from ~117 MB. Acceptable for a static-site deploy.

---

## Viz plan (genomic-regions repo)

This section is design-level. Implementation belongs to a separate session in `/Users/sam/Documents/Work/genomic-regions` after the data layer lands.

### The dictionary card (Fudoki-inspired)

The new centerpiece is a **persistent dictionary-card sidebar** that opens on click of any region (in any panel). The card replaces the deleted Step 5 cooccurrence panel and reorganizes the per-region experience into a single coherent entry.

**Card sections** (top to bottom):

1. **Headword.** Genomic location, token ID. e.g., "chr16:234,609–234,959 (token 255786)".
2. **Class.** SCREEN cclass + a one-line description ("PLS — promoter-like signature").
3. **Chromatin signature.** Narrative summary of the 18 signal-mean columns: "H3K4me3 strong (K562 17.7, HepG2 7.0, GM12878 3.5), H3K27ac present, H3K27me3 absent — classic active-promoter pattern." Optionally a small heatmap or bar chart.
4. **Senses across contexts.** Per-cell-line activity summary derived from the PMI cooc parquet's per-cell-line strata: "Active in K562 (8% of K562 files), quieter in GM12878 (2%)."
5. **Grammatical relations** — two parallel surfaces, both primary, shown side-by-side:
   - **Functional similarity (R2V kNN).** Class-coherent "see also" partners across chr16. Pre-flight (2026-04-28) showed kNN delivers strong class-coherent functional partners (80% class-match for dELS, 68% overall) even where corpus PPMI on broad strata fails. This is *not* a fallback — it's a primary surface answering "what other regions are functionally similar to this one in R2V's view."
   - **Corpus grammar (curated PMI).** Top PPMI partners from stage 12, in the most relevant curated stratum (e.g., `active_promoters_pan_cell` for promoters, `active_enhancers_pan_cell` for enhancers, `erythroid_active` if the user picks that lens). Sentence form: anchor + ordered modifiers (option A formulation).
   - The two surfaces answer **different questions** (pre-flight Jaccard between them = 0.038, near-orthogonal). Disagreement is information; the dictionary exposes both with provenance and lets the reader see where they agree vs differ.
6. **Module.** Which Leiden module(s) this region belongs to (from stage 13, per stratum), with a link to each module's full sentence (option B).
7. **Usage examples.** Top files where this region is most distinctive — derived by joining `tokenized_corpus_chr16` with `viz_files`.

The card is the **entry**. Browsing the dictionary is following links between entries.

### Step structure

Existing 5-step narrative reframes around the entry concept:

| Step | Existing | After this plan |
|---|---|---|
| 1 — peaks per file | continuous/peaks/tokens toggle (stage 14) | unchanged |
| 2 — universe + tokenization | ditto | unchanged |
| 3 — region UMAP | colored by class | unchanged for v1 (orphan `n_files_total` color toggle reads stale `region_stats.parquet`); after stage 12 lands, rewire to read `n_files_active` from `region_cooccurrence_pmi.parquet` and add a "PMI degree" toggle option |
| 4 — file UMAP | unchanged | unchanged |
| 5 — relations | (deleted) | **the card.** Click any region (in any UMAP, in any featured interval) → card opens. |

Plus optional new browsing surfaces:

- **Module catalogue** (option B): a paginated list of modules per stratum, ordered by size or by anchor's class. Each module preview shows anchor + member count + auto-label. Click → open module's sentence in the card.
- **Path explorer** (option C): pick two regions; visualize the path of pairwise similarities connecting them (PMI graph primarily, embedding fallback). Path nodes overlaid on the region UMAP as a connected polyline through 2D positions. UI copy frames it as guided browsing, not biological mechanism. Defer until A and B are landed.

**Anchor-sentence rendering modes.** Both for option A (per-region sentences) and option B (module sentences), the dictionary card renders sentences in two parallel forms:

- **Fudoki stream** — horizontal flow of class-colored, role-tagged tokens.
- **UMAP focal lens** — same tokens highlighted on the region UMAP with non-sentence tokens dimmed; optional convex-hull overlay of the partner set. Reveals geographic shape information (compact patch vs spread) the linear stream can't.

### Existing region_stats dependency

**Decision (locked in 2026-04-28):** keep the orphan `region_stats.parquet` in the bundle for now. The viz keeps working off the stale file until new stage 12 lands. At that point, rewire the "File count (log)" UMAP color toggle to read `n_files_active` from `region_cooccurrence_pmi.parquet` (per-token, `corpus_baseline` stratum). At that swap point, also consider replacing the metric with **PMI degree** (number of strong PPMI partners) — a "how grammatically connected is this region" toggle that's more native to the new framing than raw popularity.

---

## Validation tests

### Test 0 (pre-flight): R2V faithfulness on our corpus

**Before locking the curated strata library**, run a cheap pre-flight: compute corpus-wide PPMI on the `corpus_baseline` stratum only and compare to `viz_chr16.knn_*` partners. For each chr16 token, compute the Jaccard between its R2V top-30 kNN and corpus-wide PPMI top-30 partners. Distribution of Jaccards across all 35,934 tokens tells you:

- **High agreement (mean Jaccard > 0.5)** → R2V faithfully captured corpus structure. The kNN already in the bundle does most of the work; per-stratum cooc is value-add narrative richness.
- **Low agreement (mean Jaccard < 0.3)** → R2V and our curated corpus diverge significantly (probably because R2V was pretrained on a broader bedbase distribution, not our 16,799-file balanced sample). Per-stratum PMI is essential, not optional.

Either result is informative. ~10 minutes of compute. Run before implementing the full curated-strata stage 12.

### Test 1: stratification value-add

Once curated strata land, for each token compare R2V kNN partners to each stratum's PPMI partners. **Strata should diverge from R2V** — that's the point. If all strata's PPMI agree closely with R2V kNN, the strata aren't doing biology-specific work.

### Test 2 (proof of concept): α-globin cluster

The α-globin gene cluster at chr16:218,000–238,000 is the primary validation target. Expected behavior at each stage:

### Stage 12 (PMI cooc) — `active_promoters_pan_cell` and `erythroid_active` strata

For the HBA1 promoter token (255786):
- Top-K PPMI partners should include alpha-globin cluster members the kNN missed: at least one of HBA2 / HS-40 / HS-33 / HS-48 / HBQ1 promoter, ideally several. (kNN finds 255784, an AG sister PLS, as #2 — but misses HS-40 and the local cluster pELS regions; curated PMI must close that gap.)
- PPMI values for these should be substantially above 1.0 (strong association).
- The `erythroid_active` stratum should specifically reveal cluster local biology even if the file count is small (~25–40), because the question is sharp: "in erythroid active-mark contexts, what does HBA1 partner with?"
- **The decisive pass criterion**: do the curated strata recover at least one alpha-globin cluster partner that R2V kNN missed AND that `corpus_baseline` PMI missed (which had zero AG cluster members in HBA1's top-10)? If yes, the curated-strata work is unambiguously additive. If no, we need filename-heuristic backfill or a different framing.

### Stage 13 (modules) — `active_promoters_pan_cell` or `erythroid_active`, γ = 1.0

A module containing HBA1, HBA2, HS-40, HS-33, HS-48 should appear:
- Member count ~5–15 tokens.
- `dominant_class`: PLS or pELS.
- `dominant_cell_lines`: K562 prominently (erythroid lineage).
- `anchor`: ideally one of the HBA promoter tokens.

### What pass / fail looks like

**Pass:** the α-globin cluster recovers as a coherent module in `active_enhancers_pan_cell` or `erythroid_active`. HS-40 ranks among HBA1's top-3 PPMI partners. The auto-label is sensible.

**Partial pass:** HBA1's top partners are biology-coherent but the module assignment is messy (e.g., merged with an unrelated cluster). Tells us: PPMI is good; module resolution needs tuning.

**Fail:** HBA1's top partners look random or are dominated by housekeeping promoters genome-wide. Tells us: corpus or strata are inadequate. Possible causes:
- Histone strata are underpowered because the 4000 untagged generic ChIPs don't contribute (the known target gap). Mitigation: filename-heuristic backfill per Decision #5.
- Statistical floor too low → noise from rarely-active tokens dominating.
- The 16,799-file corpus might not have enough erythroid representation; the original 84,698-file pool might be needed.

---

## Failure modes worth flagging up front

### Hub dominance

PMI helps but doesn't fully eliminate broadly-active hub regions. They show up as high-degree nodes that bridge communities and can pull modules together that should be separate. Mitigations:
- Statistical floor (`min_files_active >= 5`) drops the noisy long tail.
- PPMI threshold (`> 1.0` for graph edges) drops weak associations.
- If hubs persist, can fall back to `weighted PMI` (`PPMI × log(c)`) or PMI-with-context-discount (PMI minus the median PMI of the focal token's edges).

### Resolution parameter γ

Leiden's γ is a knob, not a derived value. **For v1, hardcoded at γ = 1.0** (Locked-in decision #3); not exposed in the UI. Run the full sweep `[0.5, 1.0, 1.5, 2.0]` offline as a diagnostic if module quality looks off — modules whose membership shifts dramatically across γ are likely artifacts; modules stable across γ are robust. Revisit user-facing γ control only if v1 reveals real ambiguity that requires user judgment.

### Stratification choices

The 14 curated strata across 5 levels are a starting library, not a final list. Things that may change:
- **Add lineage-specific cross-cuts** if cross-tissue strata are too aggregated (e.g., split `active_enhancers_pan_cell` into per-cell-line versions).
- **Add per-TF strata** for prominent TFs (GATA1, CTCF, MYC) if TF-specific grammar is informative.
- **Drop underpowered strata** with fewer than ~30 files where most tokens fail the `min_files_active` floor.
- **Aggregate sparse strata** (e.g., "active marks pan-cell" combining H3K27ac + H3K4me1 + H3K4me3) if mark-specific strata are too sparse.

The validation step (α-globin pass/fail + R2V faithfulness pre-flight) should drive these choices. The strata library is research-design — keep it editable, treat each stratum as a hypothesis under test.

### Manifest gaps

`target` is missing for ATAC/DNase (expected, no target concept) and for the 4000 generic ChIP-seq files (problematic — they're likely mostly histone ChIP from labs that didn't tag them tightly). Mark-stratified strata (`active_enhancers_pan_cell`, `active_promoters_pan_cell`, `polycomb_repressed_pan_cell`) miss this 24% of the corpus.

Mitigations if validation reveals underpower:
- **Filename-heuristic backfill** — regex-mine `name` and `description` columns of the 4000 generic ChIP rows for mark strings (e.g., `H3K27ac`, `H3K4me3`). Likely recovers 50%+ of them. Cheap engineering, no re-curation.
- **Description-based fuzzy classification** — use the `description` field's text to infer marks via simple keyword matching.

Decision: defer. α-globin validation tells us whether the 60–160-file mark strata are viable as-is.

### Embedding kNN vs corpus PPMI mismatch

The R2V kNN partners (already in `viz_chr16.knn_*`) and the new PPMI partners may disagree. This is **not a bug**; the disagreement is itself a signal. Where they agree → robust grammar. Where R2V predicts a partner that corpus PPMI doesn't support → R2V is generalizing past direct evidence. Where corpus PPMI shows a partner that R2V missed → the corpus has signal R2V didn't compress. The viz should expose this contrast in the dictionary card.

---

## Order of operations

1. **Pre-flight: R2V faithfulness test.** Compute `corpus_baseline` PPMI only, compare to `viz_chr16.knn_*`. ~10 min compute. Tells us whether per-stratum PMI is value-add or essential before we commit to the full curated strata library.
2. **Implement stage 12 (curated-strata PMI cooccurrence).** Spot-check the `corpus_baseline` and `active_promoters_pan_cell` strata first; verify HBA1 partners look biological in `erythroid_active`. Iterate threshold / floor settings if not.
3. **Implement stage 13 (Leiden modules).** Run for all (stratum, γ) pairs at γ = 1.0 baseline. Confirm α-globin module recovers in at least one stratum.
4. **Compute embedding-derived features** (prototype distances, concept axes) inline. Small, additive — no separate stage needed.
5. **(When prompted by user) Run stage 06 (gene proximity)** for module auto-labels. Run stage 05 (motifs / ChromHMM) when entry-richness wants it. Both gated on user prompt; not blocking.
6. **Update bundle**: copy new parquets to `genomic-regions/src/data/dictionary/`. **Leave the orphan `region_stats.parquet` in place for now**; revisit when stage 12 lands and the viz can be rewired to pull `n_files_active` from `region_cooccurrence_pmi.parquet`.
7. **Hand off to genomic-regions session for viz changes.** That session implements the dictionary card UX on top of the new bundle. Out of scope for this plan beyond design intent.

---

## Locked-in decisions

(Open questions resolved with the user 2026-04-28.)

1. **Orphan `region_stats.parquet` in the viz** — **keep for now.** Don't strip the dependency in `src/index.md` until new stage 12 lands and the viz can rewire `n_files_total` to `region_cooccurrence_pmi.parquet`'s per-(token, `corpus_baseline`) marginals. Risk accepted: the orphan file in the bundle is stale but the viz keeps working off it.
2. **Cross-cut strata (mark × cell-line)** — **defer.** Validate the 12-stratum curated baseline first. If α-globin validation needs sharper biology, add cross-cuts as a follow-up.
3. **Resolution γ exposure to user** — **hardcode at γ = 1.0** for v1. Don't expose γ in the UI. Revisit if module quality is sensitive to γ in ways that demand user control.
4. **Stages 05 / 06 wire-up** — **prompt user to run when needed for label quality.** Proceed without for v1; auto-labels use dominant class + dominant cell-line + dominant target mark from member coactivity. Stage 06 enables nearest-gene labels; stage 05 enables TF/motif richness in entries. Both gated on explicit user prompt during implementation.
5. **Stage 01 re-run for `target`** — **no.** Manifest already has target populated for TF ChIP (100%) and Histone ChIP (91%); the gap is the 4000 generic ChIP-seq files (~24% of corpus). Corpus balance is preserved as-is. Filename-heuristic backfill is reserved as a fallback IF α-globin validation reveals histone strata are underpowered; otherwise, accept the smaller (60–160 file) per-mark cohort sizes.
6. **Strata library finalization** — the 14 starter strata across 5 levels listed above are tentative. The user owns the final strata list; revisit after the pre-flight test and α-globin validation. Levels (L1–L5) are the structural commitment; specific strata within each level can be edited.

---

## Appendix: bundle inventory after this plan lands

For reference, here's what the genomic-regions data bundle looks like after stage 12 / 13 ship.

| Parquet | Source stage | Size (est) | Used for |
|---|---|---|---|
| `viz_files.parquet` | 09 | 580 KB | File UMAP (Step 4) |
| `viz_chr16.parquet` | 08 | 28 MB | Region UMAP (Step 3), per-region signature, embedding kNN |
| `tokenized_corpus_chr16.parquet` | 11 | 62 MB | Click-any-file lookup (Step 4), live cooc queries |
| `featured_files.parquet` | 10 | 108 KB | Featured 17-file scaffolding |
| `featured_intervals.parquet` | 10 | 8 KB | The 4 narrative hotspots |
| `featured_tracks.parquet` | 10 | 8 KB | Step 1 peak data |
| `featured_signal.parquet` | 14 | (small) | Step 1 continuous/peaks/tokens toggle |
| `region_cooccurrence_pmi.parquet` | **12 NEW** | ~30–50 MB | Sentence A (anchor) per stratum, PMI graph for B and C |
| `region_modules.parquet` | **13 NEW** | ~3 MB | Sentence B (modules), per-token module id × stratum × γ |
| `module_summary.parquet` | **13 NEW** | ~200 KB | Sentence B (module catalogue / labels) |
| `region_class_prototypes.parquet` | **embedding-derived (inline)** | ~few hundred KB | Soft-classification distances per cclass for the entry card |
| `region_concept_axes.parquet` | **embedding-derived (inline)** | ~few hundred KB | Concept-axis scores per region for the entry card |
| `region_stats.parquet` | **kept as orphan** | 700 KB | Legacy `n_files_total` for the UMAP color toggle; rewire to `region_cooccurrence_pmi` after stage 12 lands |
| ~~`region_cooccurrence.parquet`~~ | (13 old, dropped) | — | Replaced by PMI version |

Net change: +~30–50 MB for new corpus-derived layers, +~1 MB for embedding-derived layers, with sharper biology and a clear narrative role for every parquet.

---

## Implementation log

### 2026-04-28 — Step 1 (pre-flight R2V faithfulness test) complete

**Script:** `genomic-dict/pipeline/scripts/preflight_r2v_faithfulness.py`
**Output:** `genomic-dict/results/preflight_r2v_faithfulness/summary.json`
**Compute:** ~6 minutes via chunked sparse matmul (200-token slabs); first attempt OOM'd because full `X.T @ X` produced 1.29B nonzero entries (~36k partners per token; chr16 cooccurrence is denser than expected) — rewritten to never materialize full matrix.

**Setup:**
- 16,755 files × 35,934 chr16 tokens
- 36M (file, token) activations
- Statistical floor: `min_files_active >= 5` (passed by all 35,934 tokens)
- Top-K: 30 partners per method

**Headline result: mean Jaccard = 0.038** (between R2V kNN top-30 and `corpus_baseline` PPMI top-30, per token).

| Jaccard bin | Count | % |
|---|---|---|
| [0.0, 0.1) | 33,084 | 92.1% |
| [0.1, 0.2) | 2,383 | 6.6% |
| [0.2, 0.3) | 393 | 1.1% |
| [0.3, 0.4) | 71 | 0.2% |
| ≥ 0.4 | 3 | 0.008% |

- Median = 0.017
- 26% of tokens have **zero overlap**
- 99.8% of tokens are below the 0.3 "low agreement" threshold

**Verdict (per plan thresholds): per-stratum PMI is essential, not optional.**

### Spot-check (HBA1 promoter, token 255786)

**Data integrity:** 35,704 distinct kNN token IDs all map cleanly to `viz_chr16.token_id`. No off-by-one, no index-vs-id confusion. Comparison is methodologically sound.

**Critical surprise:** HBA1 is active in **5,651 / 16,755 files = 33.7%** of the corpus. Not an erythroid-specific footprint — HBA1's promoter is broadly accessible across cell types because it's a CpG-island promoter, nucleosome-depleted in many contexts, peak-called in DNase / ATAC / generic ChIP files genome-wide. The transcriptional output is erythroid-specific; the chromatin accessibility footprint is not. This explains why corpus_baseline PMI fails to surface erythroid biology for HBA1: the marginal P(HBA1) is too high; PMI partners are dominated by other broadly-accessible regions, not erythroid co-actors.

**HBA1 R2V kNN top-10:** 9 PLS + 1 pELS, scattered across chr16. **Rank #2 is token 255784 (chr16:229336–229559 PLS) — a sister alpha-globin cluster member.** R2V did find local cluster biology for at least one partner. The other 9 are functionally similar PLSs elsewhere on chr16.

**HBA1 corpus_baseline PMI top-10:** mixed cclasses (CA-H3K4me3, pELS, dELS, PLS), all PMI ~1.0 (barely above chance), scattered across chr16, **zero alpha-globin cluster members**.

**Overlap between HBA1's two top-10 lists: 0.**

### Aggregate R2V kNN behavior across all chr16 tokens

- **Cosine distances:** rank-1 median 0.153, rank-30 median 0.240. Partners are tight in 100-dim space.
- **Locality:** mean 6.8% of top-30 are within 50kb of focal; 36% of regions have zero local partners. R2V kNN is **spatially blind** by design — finds functional cousins, not chromosomal neighbors.
- **Class coherence:**
  - dELS: 80% of partners share class
  - pELS: 44%
  - PLS: 39%
  - CA-CTCF: 23%
  - CA-H3K4me3: 11%
  - Overall: 68%

dELS recovery is excellent. PLS and pELS get mixed promoter/enhancer partners (biologically reasonable). CA-* classes have low coherence — a minority-class effect (only 600–900 same-class candidates available; kNN reaches into adjacent dELS/pELS pools).

### Within-cluster recovery on alpha-globin (14 tokens)

6 of 14 alpha-globin tokens find at least one other cluster member in their top-15 kNN: the dELS pair (255780/823112), the pELS pair (823113/255781), and HBA1 (255786 → 255784). The other 8 get same-class partners scattered across chr16, no cluster locality. **kNN is inconsistent at within-cluster recovery** — it finds same-type regions much better than the actual cluster siblings.

### Decisions / interpretation updates

1. **The 0.038 Jaccard is real** — not a data corruption artifact. R2V kNN and corpus PPMI are answering **different questions**, not disagreeing about the same one. The dictionary should expose both as primary surfaces.

2. **R2V kNN repositioned as primary, not "embedding aside."** Card section 5 (grammatical relations) revised to put kNN and curated-PMI side-by-side as parallel surfaces. The functional-similarity reading (kNN) and the corpus-grammar reading (curated PMI) answer different questions; both are needed.

3. **`corpus_baseline` PMI is not useful for popular promoters** like HBA1 — popularity-of-shared-broad-coverage-files dominates. This stratum is now de facto **diagnostic-only**, not user-facing for entries on broadly-accessible regions. The dictionary entry should default to the most-relevant curated stratum, not `corpus_baseline`.

4. **α-globin pass criterion sharpened.** No longer "top partners look biological." Now: "do the curated strata recover at least one alpha-globin cluster partner that R2V kNN missed AND that `corpus_baseline` PMI missed?" That's a precise testable claim about whether stratification is unambiguously additive.

5. **The kNN's class coherence is itself a finding to expose in the dictionary.** A region's entry can show "kNN partners are 80% same-class" as a confidence signal — high coherence = R2V's view is consistent; low coherence = kNN reaching into adjacent classes (informative either way).

### Plan sections touched in this log update

- TL;DR (point 6) referenced the pre-flight; result now logged here
- Card section 5 (Grammatical relations) — kNN repositioned as primary parallel surface, not "embedding neighbors aside"
- Card section 7 (Embedding neighbors) — removed (folded into section 5)
- α-globin Stage 12 pass criterion — rewritten with the kNN-vs-PMI complementarity test
- Status: `draft` → `in-progress`

### Next step

Step 2 of order of operations: **implement stage 12 (curated-strata PMI cooccurrence)**. Reuses the chunked sparse matmul from preflight; adds per-stratum file filtering and writes the output parquet. Spot-check on `active_promoters_pan_cell` for HBA1 should be the first sanity check.

---

### 2026-04-28 — Step 2 (stage 12 PMI cooccurrence) complete after three iterations

**Script:** `genomic-dict/pipeline/scripts/12_cooccurrence_pmi.py`
**Output:** `genomic-dict/data/precomputed/region_cooccurrence_pmi.parquet` (102 MB, 381,185 rows)
**Strata expanded from 14 → 18** (added L6 contrast family).
**Ranking metric changed from PPMI → NPMI** after diagnosing tied-max pathology.
**Total compute across iterations:** ~45 minutes (3× ~15 min runs).

#### v1 (PPMI, 14 strata) — exposed two unanticipated problems

First run with the originally-planned 14 strata + PPMI ranking:
- 311,369 rows, 53.4 MB
- AG cluster pairwise cohesion was 7/182 in `corpus_baseline`, **0/143** in every narrow mark-coherent stratum
- HBA1's top PMI partners showed PPMI ≈ 1.0 (barely above chance) in `corpus_baseline`, all ≈ 0.0 in `active_promoters_pan_cell`/`erythroid_active`/etc.
- **Pass criterion FAILED**: corpus_baseline found 1 AG sibling for HBA1 (255785), narrow strata found zero

Two diagnoses emerged:

**Problem A: PMI saturation for hub tokens.** HBA1 is active in 158/158 = 100% of `active_promoters_pan_cell` files and 23/24 = 96% of `erythroid_active`. When P(a) → 1, PMI mathematically collapses to 0 for every partner — narrow biology-coherent strata force their canonical-element hub tokens to saturate by construction. Stratification doesn't fix the popularity-bias problem for hub tokens; it makes it worse.

**Problem B: HBA1 misframing.** Spot-check revealed HBA1's promoter is broadly H3K4me3-marked across cell types (CpG island), not erythroid-specific in its chromatin footprint. The transcriptional output is erythroid-specific; the chromatin accessibility / mark footprint is not. So no narrow stratum can isolate "erythroid-specific HBA1 grammar" because HBA1's *chromatin* isn't erythroid-specific.

#### v2 (PPMI, 18 strata with L6 case-control contrasts) — exposed the tied-max ceiling

Added 4 case-control contrast strata as a new L6 layer:
- `active_vs_repressive_pan_cell` — H3K4me3+H3K27ac+H3K4me1 (case) ∪ H3K27me3+H3K9me3 (control); 405 files
- `erythroid_vs_other_repressive` — K562 active marks ∪ GM12878+HepG2 repressive marks; 44 files
- `lymphoid_vs_other_repressive`, `hepatic_vs_other_repressive` — symmetric

This *did* fix HBA1's marginal — she sits at 57% in `active_vs_repressive_pan_cell` (Goldilocks band). But cluster cohesion was *still* 0/143 in the L6 strata.

Diagnostic showed: HBA1 has **529 partners tied at PPMI = 0.5615** (the analytical maximum log(N/n_a) achievable for any partner whose entire active set is contained in HBA1's). All AG cluster members ranked at PPMI 0.41–0.56 — *just below* the tied ceiling — because they have one or two non-co-active files (imperfect containment). `argpartition` picks 30 from the 529-tied set arbitrarily, and AG members never make the cut.

This is the **"perfect containment ceiling"** problem: raw PMI rewards "fully contained" partners over "high statistical mass with imperfect containment," even when the relationship is biologically weaker. Tie-breaking by joint count alone doesn't help — the tied set has lots of high-joint partners.

#### v3 (NPMI, 18 strata) — pass criterion met

Switched primary ranking metric to **NPMI** (normalized PMI):

$$\text{NPMI}(a, b) = \frac{\text{PMI}(a, b)}{-\log P(a, b)}$$

NPMI normalizes PMI by joint probability, so partners with higher empirical mass (higher joint counts) outrank low-mass perfectly-contained ones. Bounded [-1, 1]; NPMI = 1 only for perfect mutual dependence with full joint coverage.

For HBA1 in `active_vs_repressive_pan_cell`:
- PPMI: 529 tied at 0.5615; AG sibling 255785 at PPMI rank 530
- NPMI: **1 partner at max (0.9923) — token 255785** (the AG sibling). AG sibling at rank 1.

**Pass criterion now satisfied.** Across all strata's NPMI top-30, HBA1's AG cluster recovery surfaces 4 of 13 siblings: {255783, 255784, 255785, 255787}. R2V kNN found only 255784. Curated NPMI added **3 new AG partners**: {255783, 255785, 255787}.

Pairwise AG-AG cohesion improved 2.5× under NPMI:

| Stratum | NPMI v3 | (PPMI v1) |
|---|---|---|
| corpus_baseline | **18/182 (10%)** | 7/182 |
| featured_lineage_HepG2 | **18/182 (10%)** | (similar) |
| tf_bound_pan_cell | **18/182 (10%)** | (similar) |
| featured_lineage_K562 | 16/182 (9%) | 3/182 |
| active_vs_repressive_pan_cell | 5/143 (3%) | 0/143 |

`tf_bound_pan_cell` was the most informative single stratum for HBA1 — it surfaced 4 AG members (ranks 1, 2, 14, 21). TF-binding cohabitation around α-globin promoters is a strong corpus signal that the broad TF stratum captures cleanly.

#### Where NPMI still fails: the "small-stratum saturation" gotcha

The L6 *lineage-contrast* strata (40–44 files each) gave HBA1 her Goldilocks marginal but *too many partners* hit NPMI = 1.0 in such small samples — perfect alignment is statistically common when N is small. So `erythroid_vs_other_repressive`, `lymphoid_vs_other_repressive`, `hepatic_vs_other_repressive` produced 0/143 cohesion despite having the right marginal structure.

The **`active_vs_repressive_pan_cell`** stratum (405 files) is the only L6 contrast that worked — it's large enough for NPMI to discriminate among partners. The tinier L6 lineage variants are not useful as primary surfaces; they're left in as diagnostic / future-work options.

#### The Goldilocks rule for stratum selection

Empirically validated:
- **Stratum size ≥ 400 files AND focal marginal in 5–60% band** → NPMI is discriminative; partners rank with meaningful spread
- **Stratum size < 50 files** → NPMI saturates at 1.0 for too many partners regardless of marginal; not useful as primary surface
- **Focal marginal > 80%** → saturation; even NPMI collapses (no variance left to normalize)

The dictionary card's per-token stratum-picker should respect this rule: prefer strata that satisfy both criteria for the focal token. For hub tokens like HBA1, the workable strata are: `corpus_baseline`, `featured_lineage_K562/GM12878/HepG2`, `tf_bound_pan_cell`, `active_vs_repressive_pan_cell`, `active_enhancers_pan_cell`. The narrow per-cell-line × mark strata and the L6 lineage contrasts are *not* useful for hub focals; they may still be useful for non-hub focals where the statistical mass works out.

#### Schema changes from v1 to v3

The output parquet now has **5 weight columns**:

| col | type | meaning |
|---|---|---|
| `weights_npmi` | List[Float32] | **Primary ranking metric.** Bounded [-1, 1]; rewards high statistical mass + association |
| `weights_ppmi` | List[Float32] | Kept for diagnostic / NLP-comparability; saturates for hub tokens |
| `weights_jaccard` | List[Float32] | Kept for cross-method comparison |
| `counts` | List[Int64] | Raw cooccurrence count |
| `partner_token_ids` | List[Int64] | Top-K by NPMI desc, joint count desc as tie-break |

Sort order is `(NPMI desc, joint count desc)` — partners are ordered for consumption.

#### Decisions / interpretation updates

1. **NPMI is the primary corpus-grammar surface.** PPMI saturates pathologically; NPMI does not. The plan's text and the dictionary entry should refer to NPMI as the ranking metric, with PPMI relegated to "diagnostic only."

2. **Plan's "popularity bias" diagnosis (Mechanism 1) was incomplete.** The full picture has three failure modes:
   - **Popularity bias** (raw counts / Jaccard reward broadly-active partners) → fixed by PMI/NPMI
   - **Modality aggregation** (mixed strata blur biology) → fixed by curated strata
   - **Tied-max / perfect-containment ceiling** (raw PMI ties for hub tokens) → fixed by NPMI

3. **The L6 contrast family is a partial success.** `active_vs_repressive_pan_cell` works well; the lineage variants are too small to be useful as primary surfaces. Keep them in the library as future-work options but don't rely on them for the dictionary card.

4. **`corpus_baseline` is significantly more useful than I positioned it.** With NPMI, it's the most consistent AG-cluster recovery surface and should be a default fallback for hub focal tokens. Earlier "diagnostic-only" framing was wrong.

5. **`tf_bound_pan_cell` is unexpectedly strong** — surfaced 4 of 13 AG partners for HBA1, including ranks 1, 2, 14, 21. The pan-cell TF binding signal captures regulatory cohabitation that mark-coherent strata can't, because TFs cluster around the same elements regardless of cell line. The dictionary should treat this stratum as a primary surface, not just a comparison.

#### Next step

Step 3 of order of operations: **implement stage 13 (Leiden modules)**. Build NPMI-weighted graph from stage 12, run Leiden community detection, validate that the AG cluster recovers as a coherent module in at least one stratum (most likely `corpus_baseline`, `featured_lineage_HepG2`, or `tf_bound_pan_cell` given their 18/182 cohesion).

---

### 2026-04-28 — Step 3 (stage 13 Leiden modules) complete

**Script:** `genomic-dict/pipeline/scripts/13_modules.py`
**Outputs:**
- `genomic-dict/data/precomputed/region_modules.parquet` (2.2 MB, 374,650 rows)
- `genomic-dict/data/precomputed/module_summary.parquet` (0.59 MB, 599 rows)

**Compute:** ~30 seconds total (across 18 strata × 1 γ). Per-stratum Leiden runs are sub-second to 4 seconds.
**Dependencies added:** `igraph 1.0.0`, `leidenalg 0.11.0` (via `uv add`).

#### AG cluster module recovery — pass criterion fully met

| Stratum | AG members in stratum | Largest-module cohesion | Module size | Dominant class |
|---|---|---|---|---|
| **active_vs_repressive_pan_cell** | 11/14 | **11/11 (100%)** | 3,664 | pELS-dominant |
| active_enhancers_pan_cell | 11/14 | 9/11 (82%) | 2,433 | pELS-dominant |
| featured_lineage_HepG2 | 14/14 | 11/14 (79%) | 5,550 | pELS-dominant |
| lymphoid_active | 10/14 | 7/10 (70%) | — | — |
| corpus_baseline | 14/14 | 9/14 (64%) | — | — |
| featured_lineage_GM12878 | 14/14 | 8/14 (57%) | — | — |
| hepatic_active | 11/14 | 6/11 (55%) | — | — |
| hepatic_vs_other_repressive | 11/14 | 6/11 (55%) | — | — |
| featured_lineage_K562 | 14/14 | 7/14 (50%) | — | — |
| lymphoid_vs_other_repressive | 10/14 | 5/10 (50%) | — | — |
| tf_bound_pan_cell | 13/14 | 6/13 (46%) | — | — |
| erythroid_active | 5/14 | 2/5 (40%) | — | — |
| open_chromatin_pan_cell | 14/14 | 5/14 (36%) | — | — |
| erythroid_vs_other_repressive | 11/14 | 4/11 (36%) | — | — |
| ctcf_boundaries | 12/14 | 4/12 (33%) | — | — |
| active_promoters_pan_cell | 3/14 | 1/3 (33%) | — | — |

**Headline:** 100% AG cohesion in `active_vs_repressive_pan_cell`, 82% in `active_enhancers_pan_cell`, 79% in `featured_lineage_HepG2`. The cluster does cohere into a single Leiden community in multiple complementary strata.

#### Resolution caveat — modules are coarse at γ=1.0

The AG cluster's module in `active_vs_repressive_pan_cell` has **3,664 total members** dominated by pELS (66%) plus PLS (25%). That's not "the α-globin module specifically" — it's "active enhancer-promoter regulatory regions across chr16 that aren't deeply repressed." AG cluster members are *inside* this larger module rather than alone in their own.

Module count per stratum at γ=1.0 reveals the resolution issue:

| Stratum | Modules | Largest |
|---|---|---|
| featured_lineage_HepG2 | 7 | 10,655 |
| featured_lineage_GM12878 | 8 | 12,077 |
| featured_lineage_K562 | 9 | 9,878 |
| corpus_baseline | 13 | 4,603 |
| active_vs_repressive_pan_cell | 17 | 3,664 |
| open_chromatin_pan_cell | 18 | 6,970 |
| active_enhancers_pan_cell | 20 | 3,885 |
| active_promoters_pan_cell | 21 | 1,394 |
| **tf_bound_pan_cell** | **266** | **6,539** |

The lineage broads at γ=1.0 produce 7–9 giant modules — too coarse for biology-specific reading. `tf_bound_pan_cell` is the outlier with 266 modules, suggesting TF binding patterns naturally form fine-grained communities (different TFs cluster around different region sets).

For finer AG-only resolution, γ would need to be 1.5–2.0+. The plan's locked-in γ=1.0 stands for v1, with γ-sweep as offline diagnostic if needed.

#### What this means for the dictionary

- The AG cluster recovery validates the framework end-to-end: corpus → PMI → NPMI → Leiden → coherent biology-grounded module.
- The module's auto-label ("pELS-dominant ... 3,664 regions, anchor 276084") is honest but coarse. Better labels need either:
  (a) higher γ to get tighter modules,
  (b) stage 06 gene-proximity wire-up so labels become "α-globin region" rather than "pELS-dominant cluster," or
  (c) hand-curated featured-region sentences (the deferred `featured_sentences.parquet` work).
- Multiple strata give different cuts at the same biology. The dictionary card should expose 2–3 stratum-specific module memberships per region, not a single canonical module — different strata reveal complementary structure (active_vs_repressive vs tf_bound vs featured_lineage_HepG2 each give a different lens on AG).

#### Decisions / interpretation updates

1. **Pass criterion for stage 13: met.** AG cluster coheres as a Leiden community at γ=1.0 in multiple strata, with 100% cohesion in `active_vs_repressive_pan_cell`. Pipeline is doing biology-coherent work.

2. **γ=1.0 produces coarse modules in some strata.** Honest framing: the Leiden modules are "regulatory neighborhoods" rather than "specific α-globin clusters." For paper-grade tight clusters, future γ-sweep is needed.

3. **Different strata produce different module structures** — 7 modules in lineage broads vs 266 in tf_bound. The dictionary should expose multi-stratum module membership so users can see complementary clusterings, not collapse to one.

4. **Module summaries are auto-labeled by dominant class + anchor.** Adequate for v1 but not user-friendly. Stage 06 (gene proximity) would dramatically improve labels — turning "pELS-dominant cluster" into "α-globin region" — and is the natural next polish step.

5. **`tf_bound_pan_cell` produces an order-of-magnitude more modules** than the other strata at the same γ. Fine-grained TF-binding community structure is a real biological signal worth surfacing in the dictionary as its own lens.

#### Next step

Step 4 of order of operations: **compute embedding-derived features** (prototype distances per cclass, concept axes). These ride alongside stages 12/13 without new pipeline stages; produces small parquets that supplement the dictionary entry's class section and add interpretable scalar scores per region. After step 4, prompt user to run stage 06 (gene proximity) for module label quality before bundle handoff.

---

### 2026-04-28 — Step 4 (embedding-derived features) complete

**Script:** `genomic-dict/pipeline/scripts/embedding_features.py` (non-stage utility — sits alongside `preflight_r2v_faithfulness.py`).
**Outputs:**
- `genomic-dict/data/precomputed/region_class_prototypes.parquet` (678 KB, 35,934 rows × 6 cols)
- `genomic-dict/data/precomputed/region_concept_axes.parquet` (720 KB, 35,934 rows × 6 cols)

**Compute:** ~3 seconds total. Loads the cached R2V model (`databio/r2v-encode-hg38`, 1,063,880 × 100 embedding matrix), slices to chr16's 35,934 tokens.

**Approach:** rather than running stage 02 + 07 just to feed this step, wrote a small standalone script that loads the R2V model directly via `geniml.region2vec.main.Region2VecExModel` and slices the embedding matrix to chr16 tokens by `viz_chr16.token_id`. Avoids the 498 MB intermediate parquet that stage 07 produces (full pretrained universe). The features themselves stay small.

#### Output 1: class prototype distances

Per-token cosine distance to each SCREEN class centroid. Example for HBA1 (token 255786, true class PLS):

| class | distance |
|---|---|
| **PLS** | **0.2995** |
| pELS | 1.0069 |
| CA-CTCF | 1.0555 |
| CA-H3K4me3 | 1.0731 |
| dELS | 1.2349 |

PLS distance is ~3× smaller than the next-closest class. The embedding strongly identifies HBA1 as a PLS, consistent with cclass label.

#### Output 2: concept-axis projections

Five interpretable scalar scores per region:

- `activity_score` — projection onto K562-active − K562-repressed direction
- `K562_specificity_score`, `GM12878_specificity_score`, `HepG2_specificity_score` — per-cell-line specificity directions (case = top H3K27ac/H3K4me1/H3K4me3 only in this cell, control = top in others)
- `anchor_score` — projection onto PLS-centroid − dELS-centroid direction (promoter-like vs distal-enhancer-like)

**Sanity checks pass cleanly:**

| cclass | n | anchor_score median | activity_score median |
|---|---|---|---|
| **PLS** | 1,770 | **+0.428** | **+0.248** |
| pELS | 7,475 | +0.040 | +0.165 |
| CA-H3K4me3 | 600 | −0.189 | −0.310 |
| dELS | 25,055 | −0.212 | −0.125 |
| CA-CTCF | 899 | −0.249 | −0.356 |

`anchor_score` runs **monotonically PLS > pELS > CA-H3K4me3 > dELS > CA-CTCF** — the major-class promoter-vs-distal gradient is clean. The minor classes (CA-CTCF, CA-H3K4me3) cluster below dELS, reflecting their distinctness from the promoter-enhancer continuum.

`activity_score` separates active classes (PLS, pELS positive) from quieter / insulator classes (CA-* negative). Distributional structure tracks biology.

**HBA1's profile:** `anchor_score = +0.67` (more promoter-like than typical PLS at +0.43), `activity_score = +0.26` (active), `K562_specificity = +0.07` (mild erythroid bias, consistent with HBA1 being broadly active but slightly more so in K562). The cell-line specificity scores are subtle in absolute value because HBA1 is broadly active across cell types — the signal is real but small, which is the right reading.

#### Decisions / interpretation updates

1. **Embedding-derived features added to bundle inventory.** Two new small parquets (`region_class_prototypes.parquet`, `region_concept_axes.parquet`); ~1.4 MB combined.

2. **`anchor_score` cleanly orders the cclasses.** This is direct evidence the embedding has learned the promoter-vs-distal-enhancer gradient as a continuous direction. The dictionary card can render anchor_score as "this region is X% promoter-like" with high confidence in the interpretation.

3. **Per-cell-line specificity scores are subtle.** Absolute magnitudes are small (~±0.05 to ±0.15) because R2V's training corpus didn't sharply separate K562/GM12878/HepG2 at the embedding level for most tokens. The viz should render these as "slight bias toward X" rather than "X-specific" — an honest reading of weak signal.

4. **Soft classification beats hard cclass labels for entries.** A region whose `distance_PLS = 0.3` and `distance_pELS = 0.4` is "85% PLS, 13% pELS" in the embedding's view, not just "PLS." For regions near class boundaries this nuance matters — the dictionary should expose the soft profile.

5. **Stage 02 / 07 not needed for this step.** The standalone script approach avoided generating the 498 MB `pretrained_universe.parquet`. Future work that needs the full universe (cross-chrom analyses, etc.) can run stages 02 + 07 then; for chr16-only the small standalone is sufficient.

#### Next step

Step 5 of order of operations: **prompt user about stage 06 (gene proximity)** before bundle copy. Stage 06 would dramatically improve module auto-labels — turning "pELS-dominant cluster" into "α-globin region module" — but it's gated on user prompt per Locked-in decision #4.

After stage 06 (or skip), step 6: **copy the new parquets to `genomic-regions/src/data/dictionary/`** for the viz handoff.

---

### 2026-04-28 — Step 5 (stage 06 reframed: ENCODE cCRE V4 → target gene evidence) complete

**Decision change:** the original stage 06 design (nearest-gene + GTEx eQTLs + GWAS Catalog) was retired. User flagged that "nearest gene" is methodologically unreliable — exemplified by FTO's chr16 enhancer (acting on IRX3 ~1 Mb away, not its nearest gene). Replaced with **ENCODE Registry V4 cCRE-Gene Links** as the Tier A regulatory-target evidence source.

**Old script:** `genomic-dict/pipeline/scripts/06_extract_extrinsic.py` → moved to `retired/old_stage_06_nearest_gene/` (kept for reference; no longer in pipeline).
**New script:** `genomic-dict/pipeline/scripts/06_target_evidence.py`.

**Data source:** SCREEN V4 cCRE-Gene Links (https://screen.wenglab.org/downloads), specifically the file `cCRE-Gene Links (3D Chromatin, CRISPR, eQTLs).zip` (463 MB). User downloaded; chr16 entries pre-extracted to `data/annotations/encode_gene_links/`.

**Outputs:**
- `data/precomputed/region_target_evidence.parquet` (10.9 MB, 1,736,956 rows)
- `data/precomputed/region_target_evidence_summary.parquet` (0.3 MB, 30,970 tokens)
- Master cCRE V4 list at `data/universe/screen_v4_2024-07.bed.gz` (31 MB gzipped)

#### Coordinate-based viz↔V4 mapping

Our viz_chr16's `accession_hex` column was joined from an earlier SCREEN release (1.97M cCREs, EH38D-form accessions). The current V4 release (2.35M cCREs) re-issued accessions: **0% accession-based overlap**. Coordinate-based matching solved this:

- **94.5%** of viz tokens match a V4 cCRE on exact (chrom, start, end)
- **4.7%** match via highest-overlap-fraction fallback (boundaries shifted ≤ a few bp between releases)
- **99.2% total mapped** (35,648 / 35,934 tokens). 286 tokens (0.8%) have no V4 cCRE match — likely deprecated cCREs.

#### Evidence-type breakdown

The V4 zip contained three evidence types; chr16 yield:

| Evidence type | Raw rows | Chr16 filtered | Matched to viz tokens |
|---|---|---|---|
| 3D Chromatin (Hi-C / ChIA-PET / HiChIP / pCHi-C) | (5.9 GB raw) | 1,754,593 | 1,233,247 (70.3%) |
| eQTL (GTEx + similar) | (2.3 GB raw) | 949,015 | 503,709 (53.1%) |
| **CRISPR perturbation** | (414 KB raw) | **0 entries** | — |

CRISPRi tile-screens have only been run on **chr3, chr10, chr12, chr19, chrX** so far (3,018 total rows in the file). chr16 has zero CRISPR-perturbation evidence — so this stage skips CRISPR for chr16-scoped output. **Worth noting in dictionary copy:** "no direct CRISPRi-validated targets available for chr16 regions yet — evidence comes from 3D contacts and eQTL associations."

#### Coverage: 86.2% of tokens have at least one Tier A evidence

- 30,970 / 35,934 chr16 tokens (86.2%) have ≥1 evidence row
- Median tokens-per-evidence ratio ~56 evidence rows per token (mean: many higher)
- The 14% of tokens without evidence are mostly distal / less-studied regulatory regions

#### Methodological nuance: what "target_gene" means depends on focal cCRE class

3D Chromatin evidence reports `cCRE → gene` where the gene is the one whose **promoter the cCRE physically contacts** in 3D space. This means:

- **Promoter cCREs (PLS) report TAD-mate genes as targets.** HBA1's promoter (token 255786, PLS) lists AXIN1 (57 datasets), NPRL3 (6), CAPN15 (2) — these are the *other* genes whose promoters HBA1's promoter contacts in 3D. HBA1 itself doesn't appear because *it is the source cCRE*.
- **Enhancer cCREs (pELS / dELS) report regulated promoters as targets.** A distal enhancer in the AG cluster correctly lists HBA1 / HBA2 / HBQ1 as targets when 3D contacts go from the enhancer to those promoters.

The dictionary card has to render this read correctly: for a PLS focal region, the evidence section says "this promoter co-contacts (in 3D space) genes X, Y, Z across N tissue contexts" — i.e., its TAD-mate gene network. For a dELS / pELS focal, same data reads as "this enhancer's predicted targets are genes X, Y, Z."

#### α-globin cluster validation: full pass

**All 14 AG cluster tokens have target evidence.** Cross-cluster co-target analysis:

- **AXIN1** is targeted by **14/14 AG tokens** — every AG token contacts AXIN1's promoter in 3D. AXIN1 sits at chr16:~340-380k, in the same TAD as α-globin. The data recovers this TAD without any TAD calling.
- 9–12 of 14 AG tokens also share NPRL3, ARHGDIG, SNRNP25, MRPL28, RGS11, POLR3K, CAPN15, DECR2, LINC00235 as targets. **TAD-mate gene network identified.**
- HBQ1 (within AG window) is correctly predicted as a target by 11/14 AG tokens (distal-enhancer cCREs in the cluster contact HBQ1's promoter).
- Reverse view: **HBA1 is predicted as a target by 9 distinct chr16 tokens across 33 evidence rows** — likely some from the LCR area (chr16:165–180k, outside our AG window). Worth investigating which tokens those are.

Biosample / tissue diversity is broad: top contributors include WTC11 (iPSC), MCF-10A/MCF-7 (breast), thyroid, T cells, HepG2, kidney, brain, liver. The 3D contacts recovered are tissue-pleiotropic, with HepG2 specifically appearing among the top 5.

#### Decisions / interpretation updates

1. **Stage 06 reframed completely.** Nearest-gene approach retired; replaced with cCRE V4 Tier A target evidence. The new framing is "predicted regulatory targets from direct measurement" rather than "proximity context."

2. **Bundle inventory updated.** New parquets:
   - `region_target_evidence.parquet` (10.9 MB) — long-format per-token-per-target-per-evidence
   - `region_target_evidence_summary.parquet` (0.3 MB) — per-token aggregate
   - Master cCRE V4 BED (31 MB gzipped) — universe reference, also useful for future stages

3. **CRISPR evidence is structurally absent for chr16.** Not a bug; the experiments haven't been run. The dictionary should disclose this honestly per region: "Tier A evidence available: 3D-chromatin, eQTL. CRISPR perturbation: none on chr16 yet."

4. **The promoter-vs-enhancer reading of target_gene needs UI discipline.** A user looking at HBA1's entry should NOT think AXIN1 "is HBA1's regulatory target" — but should see "HBA1's promoter physically contacts AXIN1's promoter (and 9 other TAD-mate genes) in 3D space, supported by 57 datasets across 53 cell-line contexts." Clear copy here is load-bearing.

5. **The TAD-mate gene network IS the dictionary's "co-regulated unit" surface.** This is what stage 13's coarse Leiden modules were trying to capture but couldn't sharpen. The 3D chromatin evidence directly recovers it. AG cluster's TAD-mate genes (AXIN1, NPRL3, etc.) become a real dictionary entry for the α-globin region.

#### Cleanup

Temp files removed:
- `/tmp/encode_gene_links/` (zip extracts)
- `/tmp/chr16_eh38e.txt`, `/tmp/crispr_chr16.tsv`
- `data/universe/screen_v4_2024-07.bed` uncompressed (kept the .bed.gz)

#### Next step

Step 6 of order of operations: **bundle copy to `/Users/sam/Documents/Work/genomic-regions/src/data/dictionary/`**. New parquets to ship:
- `region_cooccurrence_pmi.parquet` (102 MB)
- `region_modules.parquet` (2.2 MB)
- `module_summary.parquet` (0.6 MB)
- `region_class_prototypes.parquet` (0.7 MB)
- `region_concept_axes.parquet` (0.7 MB)
- `region_target_evidence.parquet` (10.9 MB)
- `region_target_evidence_summary.parquet` (0.3 MB)

Total new bundle additions: ~117 MB. Combined bundle ~234 MB (was 117 MB, growing to capture the full dictionary infrastructure).

After bundle copy, step 7: hand off to a fresh genomic-regions session for the dictionary card UX work on top of the new bundle.

---

### 2026-04-28 — Step 6+7 (bundle copy + viz proof-of-concept) shipped

**Bundle copy.** The 7 new parquets copied to `/Users/sam/Documents/Work/genomic-regions/src/data/dictionary/`. Bundle now 203 MB total. Orphan `region_stats.parquet` left untouched per Locked-in decision #1.

**Dictionary card proof-of-concept (Section 4 in `genomic-regions/src/index.md`).** Persistent card that re-renders on `pickedTokenId` change; six sections built end-to-end:

1. Headword (location + class + token id + length)
2. Class soft profile (R2V cosine to per-cclass centroids, normalized to %)
3. Concept-axis projections (5 axes, divergent bars centered on 0)
4. R2V kNN top-10 partners (Fudoki stream, class-colored chips)
5. Corpus PMI top-10 partners in `corpus_baseline` lens (NPMI-ranked, with marginal context line)
6. Tier A target evidence (3D-chromatin + eQTL counts; top genes; PLS-vs-enhancer reading caveat inline)

**Validation against the AG cluster default.** First-loaded token (255780, dELS at chr16:218,549) renders all sections cleanly:
- Soft profile: dELS dominant (28%), then pELS (24%), then CA-CTCF (18%) — biology-consistent
- Anchor score: −0.215 (correctly distal-enhancer-like)
- kNN top-1: chr16:218855-219186 dELS (the adjacent AG cluster sibling, cos=0.135)
- NPMI top-1: same AG sibling at NPMI=0.788, cooc=254 (kNN and NPMI agree on the top partner)
- Target evidence: 6 rows / 6 contexts; top genes FAM234A, AXIN1, SNRNP25 — TAD-mate gene network

**Region UMAP color toggles expanded from 3 → 9.** The four new corpus-derived encodings (anchor, activity, K562/GM12878/HepG2 specificity, evidence count) all render without errors. Each toggle has a 1–2 sentence dynamic caption explaining what story that view surfaces.

**Empirical findings worth noting from the visual exploration:**
- **Anchor & activity carve orthogonal gradients on the UMAP.** Anchor runs horizontal (PLS-rich right, dELS-rich left); activity runs vertical (active down, repressed up). Four quadrants, each biology-meaningful.
- **The dELS class (70% of chr16 tokens, one yellow blob under SCREEN coloring) splits into geographically separable sub-populations by lineage** when projected onto K562 / GM12878 / HepG2 specificity axes. K562 shows the sharpest pattern (artifact: K562 has more files in the corpus); GM12878 spreads more diffusely; HepG2 broadest. **Three different lineages, three different islands within the same dELS mass.** This is direct evidence that R2V learned cell-type-specific enhancer programs as separable directions without seeing any lineage labels in training — a substantive validation of the dictionary's distributional-semantics argument.

**Methods reference for color encodings.** `plans/2026-04-28-color-encoding-methods.md`. Documents per-toggle data sources, formulas, color schemes, validation results, and known caveats. Standalone reference doc; updated when new encodings are added to the toggle.

#### Test scripts

- `genomic-regions/scripts/check-card.mjs` — page-load smoke test; verifies all 6 card sections render and no runtime errors.
- `genomic-regions/scripts/check-card-zoom.mjs` — focused screenshot of Section 4.
- `genomic-regions/scripts/check-color-toggles.mjs` — cycles through all toggle options and screenshots each.
- `genomic-regions/scripts/check-lineages.mjs` — focused comparison of K562 / GM12878 / HepG2 specificity views.

#### Pipeline status

The pipeline-side work is **complete for v1**. The genomic-regions viz session continues to iterate on top of this stable data layer. Natural follow-ups noted in the methods doc and earlier log entries:

1. Stratum lens picker (UI control to switch the PMI partner list across the 18 strata)
2. Module section in the card (read `region_modules` + `module_summary` per stratum)
3. Usage examples section (DuckDB query against `tokenized_corpus_chr16` for top files where this region is most distinctive)
4. Wire `n_files_total` color toggle to read from `region_cooccurrence_pmi`'s `corpus_baseline` marginals (deprecates orphan `region_stats.parquet`)
5. Out-of-corpus extension spike (client-side R2V via ONNX/transformers.js for arbitrary GWAS-hit lookup)

---

### 2026-04-28 — Tier 1 viz additions: focal lens + module color + senses bars

Three preliminary visualizations of the new pipeline elements, building on the proof-of-concept card:

#### (1) UMAP focal lens — partner overlays

When a region is picked, two new layers light up on the region UMAP:
- **Blue rings** = top-30 R2V kNN partners (synchronous lookup from `viz_chr16.knn_token_ids`)
- **Orange rings** = top-30 NPMI partners in `corpus_baseline` lens (async DuckDB query against `region_cooccurrence_pmi`)

Concentric blue+orange = robust grammar (both surfaces agree). Solo color = diagnostic (R2V-smoothed transitive vs raw corpus cooccurrence). The plan's "two parallel surfaces" methodological argument is now visible spatially, not just listed in the card.

Implementation: two new `vg.Selection.intersect({empty: true})` selections (`knnPartnersSel`, `npmiPartnersSel`); a bridge that updates them when `pickedTokenSel` changes (sync for kNN, await for NPMI); two new `vg.dot` layers in the regionUmap construction with distinct stroke colors. Layered before the picked-ring so the black ring stays on top.

#### (2) Module ID color toggle (10th option)

Region UMAP can now color by Leiden module ID from the `active_vs_repressive_pan_cell` stratum at γ=1.0. Module IDs cast to strings for categorical color via `tableau10` (cycles for 17 modules — acceptable visual).

Validated: Leiden modules ARE spatially coherent on the embedding even though the algorithm never saw UMAP coordinates. Different geographic regions of the UMAP correspond to different module IDs with reasonable boundaries. Tokens missing from the stratum (failed statistical floor) show as "—" gray.

The full corpus → NPMI graph → community detection chain → coherent UMAP geography is now one click away.

#### (3) Per-stratum activity bars (Section 4 of card — "Senses across contexts")

The previously-empty Section 4 now renders 18 horizontal bars, one per stratum, ordered L1 → L6. Each bar shows the focal region's marginal activity rate (`n_files_active / n_files_in_stratum`) as a percent, with bar color encoding the stratum level.

For the AG dELS default token: bars max out in `active_enhancers_pan_cell` and `ctcf_boundaries` at low single-digit percent (rare distal element); narrow strata like `erythroid_active` / `lymphoid_active` / `active_promoters_pan_cell` are grayed as "below floor" — biologically consistent.

Reads directly from `region_cooccurrence_pmi.parquet` per-stratum marginals. The card now answers "which biological contexts is this region active in?" at a glance.

#### Files touched

- `genomic-regions/src/index.md` — additions only; no breaking changes
  - Added `moduleIdByToken` materialization cell (DuckDB query for active_vs_repressive modules)
  - Augmented `classedChr16` with `module_id_avr` (categorical string)
  - Added `Module ID (active_vs_repressive)` to color radio + caption + colorChannel/Directives branches
  - Added `knnPartnersSel` + `npmiPartnersSel` Selections + bridge
  - Added two new `vg.dot` layers in `regionUmap` for partner overlays
  - Added `marginalsByStratum` query + `sensesHTML` block in `dictCard`
  - Updated UMAP caption to mention partner overlay legend
- `plans/2026-04-28-color-encoding-methods.md` — methods doc still accurate; module encoding documented
- `genomic-regions/scripts/check-module.mjs`, `check-focal-lens.mjs` — new test scripts

#### What's left from the menu

Tier 2: stratum lens picker, module catalogue page, kNN-vs-NPMI overlap visualization
Tier 3: featured-region tour, HBA1 static showcase, methods histograms
Tier 4: path explorer, stratum heatmap, out-of-corpus extension

---

### 2026-04-28 — Tier 2 viz additions: stratum lens picker + module catalogue

Two more preliminary visualizations of the new pipeline elements, building on the Tier 1 trio. Skipped the kNN-vs-NPMI overlap viz per user direction.

#### (4) Stratum lens picker

Single shared `currentLens` view (`Inputs.select` of all 18 strata) drives:
- The dictionary card's Section 5 (corpus-grammar partners) — async DuckDB query against `region_cooccurrence_pmi` filtered to the selected stratum
- The orange (NPMI) partner overlay on the region UMAP — partner-overlay bridge cell now depends on `currentLens` so it re-fires on lens change
- The module catalogue (Section 5) — filters module list to the selected lens

Default: `corpus_baseline` (broad, Goldilocks for most tokens).

The card's PMI section now also surfaces a **saturation warning** inline when the focal token's marginal in the selected lens ≥ 80% — text appears in red explaining that NPMI is near-flat across partners and the partner ordering is mostly arbitrary, with an explicit suggestion to switch to a broader lens. Direct in-UI honesty about a real methodological pitfall (the same one we hit during pipeline validation when active_promoters_pan_cell saturated for HBA1).

The "below floor" case also gets a graceful message when the focal token failed `min_files_active >= 5` in the selected lens — instead of an empty section, the user sees "this region is below the statistical floor in this lens; try a broader one."

Validated: switching from `corpus_baseline` (13 modules) to `active_vs_repressive_pan_cell` (17 modules) — both numbers match pipeline validation. The card's PMI list, the orange UMAP rings, and the module catalogue all re-render coherently.

#### (5) Module catalogue (Section 5)

A new section at the end of the page that lists all Leiden modules in the currently-selected lens at γ=1.0, sorted by size (largest first). Each row shows:
- Module ID
- Anchor token (region + cclass chip)
- Dominant class label
- Member count + class breakdown string

Color-coded left border by dominant class (yellow for dELS, orange for pELS, red for PLS, blue for CA-CTCF, pink for CA-H3K4me3). Click a row → `pickedTokenSel` updates to the module's anchor token → dictionary card re-renders for that anchor; smooth-scrolls back to Section 4.

Reads `module_summary.parquet`'s 599 rows once at load; filters in JS by current lens.

Validated: lens-switch produces correct module count per pipeline data (13 baseline → 17 active_vs_repressive → likely 7 in lineage broads → 266 in tf_bound). The catalogue is the natural complement to the per-region card — browse by module rather than by region.

#### Files touched

- `genomic-regions/src/index.md`
  - Added `currentLens` view (`Inputs.select` of 18 strata) at top of Section 4
  - Updated `dictCard` to use `currentLens` in cooc query and surface lens / saturation / below-floor states in copy
  - Updated partner-overlay bridge to depend on `currentLens` (orange UMAP ring follows the picker)
  - Added new `## 5. Browse the module catalogue` section: materialization cell + render cell with click-to-pick handlers
- `genomic-regions/scripts/check-catalogue.mjs`, `check-lens-switch.mjs`, `inspect-lens.mjs` — new test scripts

#### What's left from the menu

Tier 3: featured-region tour, HBA1 static showcase, methods histograms
Tier 4: path explorer, stratum heatmap, out-of-corpus extension
Skipped (per user): Tier 2 kNN-vs-NPMI overlap visualization

---

### 2026-04-29 — Tier 3 first item: HBA1 worked-example showcase

A pinned static showcase at the top of Section 4, before the lens picker and live card. Pulls the same data the live card uses but always renders for token 255786 (HBA1 promoter), with inline callouts explaining what to read into each section.

#### Design choices

- **"WORKED EXAMPLE" badge** in dark blue — visually distinguishes the showcase from the live card below.
- **Two-column layout** for compactness: left column = class soft profile + concept axes; right column = senses across contexts + top target evidence. Grammatical relations span both columns at the bottom.
- **Inline callouts** (blue-bordered, italic) follow each section, explaining what the values mean and what teaching point they make. Five callouts total: class, axes, senses, grammar, target.
- **AG badge highlighting** — partner chips for AG cluster members (chr16:218k-238k tokens) get a yellow "AG" badge appended, making cluster recovery visible at a glance.
- **Bottom summary line** — "What to read out of this:" — positions the showcase as a teaching artifact, then explicitly invites the user into the live card.

#### Content callouts (what each section is meant to teach)

1. **Class profile** — PLS at 84% confidence; the soft profile and the categorical SCREEN label agree. Mild pELS overlap is biologically reasonable given HBA1 sits adjacent to the α-globin enhancer cluster.
2. **Concept axes** — anchor = +0.67 (above typical PLS median +0.43; the embedding sees HBA1 as more promoter-like than average). K562 specificity = +0.07 (small but positive — real erythroid bias). Read together: "broadly accessible promoter, slightly more so in K562."
3. **Senses across contexts** — saturated at 100% in active_promoters_pan_cell (the case the live card's saturation warning catches); 96% in erythroid_active; 80% in lymphoid_active and hepatic_active. Below floor in polycomb_repressed_pan_cell. The transcriptional output is erythroid-specific; the chromatin footprint isn't.
4. **Grammatical relations** — kNN top-30 includes N AG cluster members; corpus_baseline NPMI top-30 includes M AG members; together they recover X of 13 AG siblings. The two views are complementary, not redundant: kNN finds smoothed-similarity partners, NPMI finds direct-cooccurrence partners. Their agreement on AG cluster siblings is the dictionary's grammar holding up to triangulation.
5. **Target evidence** — top "target" is AXIN1 (~57 datasets). Methodological caveat: HBA1 is a PLS, so target_gene reads as "TAD-mate genes" not "regulated targets." AXIN1, NPRL3, SNRNP25 etc. are TAD-mates; the 3D-chromatin data recovers the TAD without ever computing TADs explicitly.

#### Why this matters before the live card

The live card is interactive and renders different content depending on what the user clicks. Without an anchor example first, a user clicking around might not know:
- Which sections matter
- How to read the numbers
- What "good" looks like
- Where the methodological caveats live (PLS-vs-enhancer reading, saturation, below-floor)

The worked example gives readers a mental model before they start exploring. Once they've absorbed HBA1's entry, the live card becomes "I know what these sections are; what does *this other region* show?"

#### Files touched

- `genomic-regions/src/index.md` — new `hba1Showcase` cell at top of Section 4 (before lens picker); 200+ lines of HTML rendering logic. Same data sources as the live card; values stay accurate as parquets evolve.
- `genomic-regions/scripts/check-showcase.mjs` — screenshot test script.

#### What's left from the menu

Tier 3: featured-region tour, methods histograms
Tier 4: path explorer, stratum heatmap, out-of-corpus extension

---

### 2026-04-29 — Fudoki-inspired addition: module-as-sentence rendering

After a brainstorm on which Fudoki UI patterns would translate cleanly to genomics, the user picked the most-novel item from the menu and rejected the genome-track-style reading modes (genomic order isn't the right ordering for regulatory grammar — that's the project's own thesis).

The shipped feature: **clicking a row in the module catalogue (Section 5) now expands inline to show the module's members as a flowing stream of class-colored chips**, ordered by SCREEN class hierarchy (PLS → pELS → dELS → CA-H3K4me3 → CA-CTCF → unclassed) then by within-module eigenvector centrality. The Leiden community becomes a "regulatory sentence" — anchor + class-grouped co-members read out in grammatical order rather than genomic order.

#### Design choices

- **Anchor highlighted with ★** (highest within-module centrality). Visually distinct from other chips.
- **Class-colored left borders** on each chip — matches the existing chip vocabulary (partner chips in the card, partner overlays on the UMAP). Creates visual consistency: a pELS chip looks the same wherever it appears in the dictionary.
- **Each chip is clickable**. Click → updates `pickedTokenSel` to that member → live card re-renders with that member's entry → smooth-scrolls back to Section 4. Same handler infrastructure as the live UMAP click.
- **Toggle behavior**: clicking the same row collapses; clicking a different row collapses the previous one and expands the new one. Single-expansion-at-a-time.
- **Arrow indicator** ▸ rotates to ▾ when expanded.
- **Top-50 cap**: large modules (3,000+ members at γ=1.0) would be unreadable as a single stream. Cap at top-50 by centrality desc — the most "anchor-like" members. Footer text states "50 of N members shown."
- **Italic caveat below the chips**: "Order is grammatical (class hierarchy + centrality), not genomic" — keeps the methodological framing visible at the point where users might assume otherwise.

#### Why this is the Fudoki adaptation that earned its keep

The earlier brainstorm had options that translated Fudoki more literally — render BED files as colored token streams in genomic order, multi-document switcher, etc. The user rejected those because the project's whole argument is that genomic order is **not** grammatical order. A Fudoki-faithful BED-file reader would have been a tokenized track browser, undermining the dictionary's distributional-semantics thesis.

Module-as-sentence keeps the Fudoki-inspired UX (a clickable colored token stream rendered as a "sentence") while preserving the project's methodological commitment: the module is a corpus-derived community whose ordering is *grammatical* (class roles + centrality), not chromosomal. This is the dictionary's argument made directly visible — a "sentence" of regulatory regions whose structure was learned from cooccurrence, not position.

#### Files touched

- `genomic-regions/src/index.md` — replaced the catalogue row click handler with a richer one that toggles inline expansion. Added `renderModuleSentence` helper. ~100 lines net change.
- `genomic-regions/scripts/check-module-sentence.mjs` — new test script.

#### What's left from the menu

Tier 3: featured-region tour, methods histograms
Tier 4: path explorer (deferred, low value-to-effort), stratum heatmap, out-of-corpus extension
Skipped (per user): the rest of the Fudoki brainstorm (full BED reading mode, multi-document switcher, featured-interval reading panes, word-in-context view)
