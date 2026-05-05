---
date: 2026-04-26
status: in-progress
description: Migrate the working React/Vite "regulatory dictionary" viz from `mattress/` to a fresh Observable Framework repo, using Mosaic/vgplot for interactive brushing + cross-plot selection. Course deliverable target — submit a hosted static-site URL. Self-contained plan; designed to execute in a compacted/fresh session.
---

## Decisions locked in (2026-04-27 post-handoff data refresh)

8. **Bundle refresh from `mattress/public/data/dictionary/`** — 8 parquets, ~117 MB. Three new ones (`region_stats`, `region_cooccurrence`, `tokenized_corpus_chr16`) plus updated existing five. Featured files now 17 (was 21) — asymmetric grid, GM12878 missing H3K9me3.
9. **Step 5 = ego network** centered on a clicked/selected region, with a single context selector (`all / K562 / GM12878 / HepG2`). Small-multiples per-context view deferred as a polish step.
10. **Step 4 click-any-file** — promote from the 17-file dropdown to any of the 16,799 corpus files via lazy DuckDB-WASM query against `tokenized_corpus_chr16.parquet`. Accept the ~5s first-load cost; subsequent queries sub-second.
11. **Featured grid GM12878 gap** — derive rows from data (already does), drop the row entirely rather than render a placeholder. One-line note in copy.
12. **Step 3 kNN disclaimer** — one-sentence copy fix to distinguish `knn_*` (embedding similarity) from `region_cooccurrence` (corpus co-activation).

## Decisions locked in (2026-04-26 pre-execution Q&A)

1. **Repo state**: empty. Run the Framework init in step 1.
2. **Deploy target**: GitHub Pages. Site will be served at `https://sanghoonio.github.io/genomic-regions/`, so `observablehq.config.js` needs `root: "src"` (default) and the build needs to emit assets at the `/genomic-regions/` base path. Add `.github/workflows/deploy.yml` (Framework build → `actions/upload-pages-artifact` → `actions/deploy-pages`).
3. **Polish**: 1-week budget; polish regardless. Don't cut corners on narrative copy or styling.
4. **Click-to-select on file UMAP**: ship in v1 using Mosaic `vg.Selection` + a point-click handler on the file UMAP `vg.dot` mark. Selecting a file row drives the same `selectedFiles` selection that the brush populates; downstream plots filter or highlight off that single shared selection.
5. **`viz_chr16.parquet` size (29 MB)**: accept for v1. Revisit (column pruning or build-time data loader) once the full viz is finalized.
6. **vgplot vs plain Plot for Step 1/2**: prefer vgplot for cross-plot selection participation; fall back to plain `@observablehq/plot` if the universe-row + 21-file-row combined layout fights vgplot. Document the fallback if used.
7. **Scrollytelling**: dropped. Plain stacked markdown sections.

# Observable Framework migration — regulatory dictionary viz with Mosaic/vgplot

## Goal

Port the working 5-step narrative regulatory-genomics dictionary from `mattress/src/dictionary/*.tsx` to a clean Observable Framework project. Use Mosaic vgplot (`@uwdata/vgplot`) instead of plain `@observablehq/plot` so we get **interactive brushing and cross-plot selections for free**:

- Brush a region of the file UMAP → filters all the file rows in Steps 1/2
- Brush a region of the region UMAP → highlights tokens at all featured intervals
- Click a featured BED file → highlights its activated regions on the region UMAP
- Linked filtering by SCREEN class / assay across all plots

Deliverable: hosted static site (Vercel / GitHub Pages) URL, suitable for course submission.

## Repo

**Target repo**: https://github.com/sanghoonio/genomic-regions.git

Clone locally before executing:

```bash
cd /Users/sam/Documents/Work/
git clone https://github.com/sanghoonio/genomic-regions.git
cd genomic-regions
```

The fresh session should treat this clone as the working directory. Inspect what's already in the repo before initializing — if it's bare, run the Framework init (step 1 below). If something's already scaffolded, work with it rather than overwriting.

Why a separate repo (not part of `mattress/`):

- `mattress/` is a general genome-exploration tool with its own scope; the dictionary viz shouldn't entangle with that
- Framework projects have their own structure conventions (`src/`, `observablehq.config.js`, `npm run build`)
- Course submission benefits from a clean public repo + deployed-site URL pair

## Reference materials (paths the executing session will need)

### Working React port (port-from)

`/Users/sam/Documents/Work/ai-sandbox/workspaces/sam/bedbase/repos/mattress/src/dictionary/`

- `DictionaryApp.tsx` — main layout, interval picker, scrollytelling structure
- `data.ts` — DuckDB + Arrow conversion + parquet registration
- `colors.ts` — SCREEN class palette + assay palette
- `types.ts` — row-shape types (port loosely; Framework is JS not TS)
- `Step1Peaks.tsx` — peak rectangles per file at current interval
- `Step2Tokenize.tsx` — universe band + per-file token activations
- `Step3Embedding.tsx` — region UMAP scatter
- `Step4CrossLink.tsx` — file UMAP ↔ region UMAP highlight via dropdown
- `Step5Hypothesis.tsx` — placeholder

### Source data (parquets to copy in)

`/Users/sam/Documents/Work/spatial-region-features/genomic-dict/data/precomputed/`

| file | rows | size | what it has |
|---|---|---|---|
| `viz_chr16.parquet` | 35,934 | ~29 MB | per-region: token_id, region (`chr:start-end`), chrom, start, end, overlaps_screen, accession_hex, cclass, umap_x, umap_y, knn_token_ids (list), knn_distances (list), 18 bigwig mark/ATAC mean scalars |
| `viz_files.parquet` | 84,718 | ~3 MB | per-BED-file: id, name, description, umap_x, umap_y, assay, cell_line, cell_type, tissue, is_unlabeled |
| `featured_intervals.parquet` | 4 | <10 KB | curated chr16 narrative loci: interval_id, chrom, start, end, label, narrative_caption, n_universe_tokens, universe_token_ids (list) |
| `featured_files.parquet` | 21 | ~75 KB | per-featured-file: file_id, name, assay, cell_line, target, role (`featured`/`mystery`), n_chr16_active_tokens, chr16_active_token_ids (list) |
| `featured_tracks.parquet` | ~70 | <10 KB | per (file, interval): file_id, interval_id, n_active_tokens, active_token_ids (list), n_peaks, peak_starts (list), peak_ends (list) |

### Background plan

`/Users/sam/Documents/Work/spatial-region-features/plans/2026-04-22-regulatory-dictionary-viz.md` — full project context, especially the 2026-04-26 update describing the pretrained-model pivot and 5-step narrative.

## Tech stack

- **Observable Framework** (latest) — markdown-first, builds to a static site
- **Mosaic vgplot** (`@uwdata/vgplot`) — reactive coordinator + plot DSL with built-in brushing/selection
- **DuckDB-WASM** — included transitively via Mosaic
- **@observablehq/plot** — only if needed for marks vgplot doesn't expose cleanly (e.g., custom HTML for entry cards)

## Why Mosaic instead of plain Plot

Plain `Plot.plot()` is a static rendering primitive — to wire interactivity you hand-roll DOM event listeners and re-render on state changes (what we do in mattress with `useEffect` + `replaceChildren`). With Mosaic:

- `vg.Selection` is a shared filter state — multiple plots can subscribe
- `vg.intervalX(selection)` mark turns any plot into a brushable filter source
- Plots that read `filterBy: selection` automatically re-render when the brush changes
- DuckDB is the data layer — queries are reactive on selection

This is what makes the "brush the file UMAP, watch the region UMAP filter" interactions cheap to add.

## Section structure (port from mattress)

Mirror the 5-step narrative with the same intervals and labels.

### Section 0 — header + interval picker

- Title and intro paragraph (port from `DictionaryApp.tsx` header)
- `viewof currentIntervalId = Inputs.select(...)` — dropdown of 4 featured intervals
- Display interval label + chrom:start-end + narrative caption

### Section 1 — Step 1: peaks at the featured interval

- Per-file horizontal-bar plot, x = genomic coord clipped to interval, y = file label, color = assay
- Width-rendering rule: peaks ≥200 bp as `vg.rectX` (or `Plot.barX`), narrower omitted
- Same y-domain as Step 2 (all 21 featured files; rows align across plots)
- Tooltip shows file label + peak coords + width

### Section 2 — Step 2: universe + per-file token activation

- Single plot with universe row at the top of the y-domain (label `— UNIVERSE —`) and 21 file rows below
- Universe band: cCREs from `viz_chr16` filtered to current interval, colored by SCREEN class
- File rows: tokens activated by each file in the current interval (same SCREEN class colors)
- Boxes only (≥200 bp), `clip: true` so boundary-spanning items don't bleed into the y-axis margin
- Same x-axis as Step 1 → tokens align with peaks above

### Section 3 — Step 3: region UMAP (chr16, all 35,934 tokens)

- `vg.dot(viz_chr16, {x: "umap_x", y: "umap_y", fill: "cclass"})`
- Brushable: `vg.intervalXY(regionSelection)` — lets the user lasso a UMAP cluster
- Color by SCREEN class with the standard palette
- Outlined ring on tokens within the current interval (ties back to Step 2)
- Hover tooltip with region coord + class + token id

### Section 4 — Step 4: file UMAP ↔ region UMAP cross-link

- Side-by-side scatters
- File UMAP: `vg.dot(viz_files, {x, y, fill: "assay"})`, mystery files (is_unlabeled=true) as triangle outline
- Brushable on either side; the brush filters the *other* plot
- Selecting a single file (via dropdown for now; click-select can be added later) overlays that file's `chr16_active_token_ids` as larger/highlighted dots on the region UMAP

### Section 5 — Step 5: hypothesis generation

- Placeholder for now (see `plans/2026-04-26-region-interpretation-step5.md`)
- Two text-generation cells: kNN class aggregation for a selected region; file-UMAP-kNN aggregation for a selected mystery file

## Mosaic-specific interactivity to add (the why-Mosaic payoff)

These are the moves we *can't* get cleanly with plain Plot:

| interaction | implementation sketch |
|---|---|
| Brush region UMAP → highlight in the universe band | shared `vg.Selection` between Step 3 dot mark and Step 2 rect mark; Step 2 reads `filterBy: regionSel` to dim non-selected tokens |
| Brush file UMAP → restrict Step 1/2 file rows | shared `vg.Selection` filtering `featured_tracks` by `file_id ∈ selectedFiles` |
| Click a SCREEN class in the legend → keep only that class everywhere | Mosaic's legend is auto-interactive when bound to a selection |
| Cross-zoom: pan/zoom in Step 1 → zoom Step 2 | `domainX` linked via shared selection |

Reference: https://idl.uw.edu/mosaic/api/vgplot/inputs.html

## Migration mapping (mattress → Framework)

| in mattress (React/TS) | in Framework (Observable JS) |
|---|---|
| `useState`/`useMemo` | `mutable`/cell-level reactivity |
| `useEffect` + `useRef` + `replaceChildren` | cell returns DOM node directly; runtime handles re-render |
| `<Step1Peaks ...props />` | one cell per visual: `step1 = vg.plot(...)` |
| `useDictionaryData` (DuckDB load + Arrow conversion) | one cell per registered table: `regions = await db.query(...)` |
| `arrowToRows` helper (Arrow Vector → plain JS array, BigInt → Number) | **same function needed** — Mosaic's `coordinator.query` returns Arrow Tables; same Vector pitfall |
| `mutable currentIntervalId` (in DictionaryApp) | `viewof currentIntervalId = Inputs.select(...)` |
| Tailwind classes | `html` template literal with inline styles, or skip styling |

## Steps to execute (in order)

1. **Set up repo** — clone `https://github.com/sanghoonio/genomic-regions.git`. Check what's there:
   - If empty / no `package.json`: run `npx @observablehq/framework@latest create .` inside the clone, pick "Empty" template, accept defaults. Confirm `npm install && npm run dev` shows the default Framework landing page.
   - If something's already scaffolded: read it, work with the existing structure, don't overwrite without confirming.

2. **Add Mosaic** — `npm install @uwdata/vgplot @uwdata/mosaic-sql @uwdata/mosaic-core`. The `vg` namespace will be the main API.

3. **Drop in parquets** — copy from `spatial-region-features/genomic-dict/data/precomputed/` to `src/data/dictionary/` (a directory under Framework's `src/`). Use the `FileAttachment` API to load.

4. **Set up Mosaic coordinator + register tables** — top of `src/index.md`. Create vgplot's coordinator backed by DuckDB-WASM, then `CREATE TABLE` from each parquet. Provide the same `arrowToRows` helper for any cells that read raw rows.

5. **Port colors + types** — JS module `src/components/dictionary.js` (or inline) with `SCREEN_CLASS_COLORS`, `ASSAY_COLORS`, etc.

6. **Port section by section** — each section is a markdown block followed by a JS code fence. Reproduce Step 1 first (peaks) end-to-end, then Step 2, then Step 3, then Step 4, then Step 5.

7. **Add Mosaic interactivity** — once each step is rendering statically, weave in `vg.Selection` for cross-plot filtering. Test brushing.

8. **Narrative copy + polish** — port the markdown captions and any styling.

9. **Deploy** — `npm run build`, deploy the `dist/`:
   - GitHub Pages is most natural since the repo is on GitHub. Add a workflow at `.github/workflows/deploy.yml` that runs the Framework build on push to `main` and publishes `dist/` via `actions/upload-pages-artifact` + `actions/deploy-pages`.
   - Or Vercel: `vercel --prod` from the clone, point at `dist/` as the output directory.
   - Submit the deployed URL (`https://sanghoonio.github.io/genomic-regions/` or the Vercel URL) as the course deliverable.

## Caveats / gotchas (from porting experience)

- **Arrow Vector vs JS array.** DuckDB-WASM list columns return as Arrow `Vector` objects, not arrays. `Array.isArray(v)` is false; iterating with `for (let i; i < v.length; i++) v[i]` returns undefined. Convert with `Array.from(v)` or `v.toArray()`. Use the `arrowToRows` helper from `mattress/src/dictionary/data.ts:27`.

- **BigInt vs Number.** parquet `int64` becomes JS `BigInt`. Coerce to `Number` for arithmetic and Plot. The `arrowToRows` helper handles this.

- **Sub-pixel rectangles.** Genome-scale plots (chr16 windows of 20–500 kb) make small peaks (200–700 bp) sub-pixel-wide. Filter to peaks ≥200 bp; render only those. At 20 kb / 900 px = 22 bp/px, a 200 bp peak is ~9 px wide — visible.

- **Y-axis row alignment.** Compute the file-label list once at the top of the notebook (alphabetical sort of all featured + mystery file labels) and pass to all plots that share the y-axis. Don't let each plot derive its own labels — they'll diverge.

- **Clip overflow.** Peaks/tokens that span the interval boundary will visually leak into the y-axis margin without `clip: true` on the mark.

- **Tokenizer vocabulary mismatch.** The data already enforces this — `featured_files.parquet` and `featured_tracks.parquet` contain pretrained `databio/r2v-encode-hg38` token IDs that align with `viz_chr16.parquet`. Don't try to rebuild tokenization in the notebook; use the parquet IDs as-is.

## Featured intervals (locked in)

| id | window | label | narrative |
|---|---|---|---|
| `alpha_globin_genes` | chr16:218,000–238,000 (20 kb) | α-globin gene cluster | HBA1/HBA2/HBQ1, well-characterized active-chromatin domain. K562 active. |
| `ciita_promoters` | chr16:10,960,000–10,980,000 (20 kb) | CIITA alternative promoters | pI/pIII/pIV alternative promoters. GM12878 (B-cell) hits pIII. |
| `cdh1_promoter` | chr16:68,725,000–68,745,000 (20 kb) | CDH1 / E-cadherin promoter | TSS at 68.737M; HepG2 active. Methylation silences in tumors. |
| `fto_enhancer` | chr16:53,760,000–53,780,000 (20 kb) | FTO intron 1 enhancer | rs1421085 disrupts ARID5B → derepresses enhancer → activates IRX3/IRX5 ~1 Mb away |

## Color palettes (locked in)

```js
SCREEN_CLASS_COLORS = {
  PLS:          "#ff0000",
  pELS:         "#ffa700",
  dELS:         "#ffcd00",
  "CA-CTCF":    "#00b0f0",
  "CA-H3K4me3": "#ffaaaa",
  unclassed:    "#cccccc",
}
ASSAY_COLORS = {
  "ATAC-seq":         "#1f77b4",
  "DNase-seq":        "#2ca02c",
  "ChIP-seq":         "#9467bd",
  "TF ChIP-seq":      "#e377c2",
  "Histone ChIP-seq": "#ff7f0e",
}
```

## Out of scope for this plan

- **Stage 05** (sequence features: FIMO, ChromHMM, phastCons) — would enrich the entry card; not blocking the demo
- **Stage 06** (extrinsic annotations: nearest gene, GWAS, eQTL) — same, not blocking
- **Click-to-select on the file UMAP** — Mosaic's interval brush replaces this need; can implement single-file picking via dropdown or single-point click
- **Step 5 hypothesis generation** — keep the placeholder, defer; see `plans/2026-04-26-region-interpretation-step5.md`

## Estimated effort

3–5 hours of focused work in a fresh session, broken down:

- 30 min: repo init + dependencies + parquet copy
- 30 min: Mosaic coordinator + table registration + `arrowToRows`
- 60 min: Step 1 + Step 2 (the trickiest — y-axis alignment, single combined plot)
- 30 min: Step 3 region UMAP
- 60 min: Step 4 file UMAP + cross-link with Mosaic selections
- 30 min: Step 5 placeholder
- 30 min: narrative copy + intro/outro
- 30 min: deployment

## Status / next steps

**Draft.** No code written yet. Execute in a fresh Claude session by:

1. Reading this plan
2. Cloning the target repo: `git clone https://github.com/sanghoonio/genomic-regions.git`, inspect any existing scaffold before adding/overwriting
3. Reading the mattress source files listed above (port-from)
4. Reading any background context that's relevant (the 2026-04-22 plan, especially the 2026-04-26 update at the bottom)
5. Walking through the 9 steps in order

The mattress port is the source of truth for the visual decisions. The plan above is the destination spec. Quiz user on Mosaic-vs-plain-Plot tradeoffs and the data-loading approach before executing if there's any ambiguity.
