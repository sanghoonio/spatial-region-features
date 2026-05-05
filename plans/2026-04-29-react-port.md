---
date: 2026-04-29
status: draft
description: Fresh React (Vite + TS + Tailwind v4 + DaisyUI v5) build at the genomic-regions repo root, scoped to a curated set of featured intervals — not a port of the Observable Framework pages. Observable kept as a sibling backup.
---

# React build for genomic-regions viz

## Goal

Build a fresh React app at `/Users/sam/Documents/Work/genomic-regions/` that focuses the visualization on a curated set of featured intervals (same scope as Reference page Section 1 — α-globin, CIITA, CDH1, FTO, plus the candidate hubs from Canvas 1: 21.4 Mb, 88.2 Mb, 15.1 Mb). Don't port the Observable pages 1:1 — start from a clean slate that combines the most useful pieces (interval-scoped Section 1, region UMAP with dict card, chr16 distribution probe) into one focused demo.

Keep the existing Observable Framework version as a sibling directory `/Users/sam/Documents/Work/genomic-regions-observable/` (untouched copy).

Why fresh: assignment allows React; no value in maintaining bug-compatible parity with three Observable pages that overlap heavily; the scoped demo is more focused and easier to evaluate.

## Repository layout

**Before:**

```
genomic-regions/
├── src/                       # Observable Framework
│   ├── index.md
│   ├── canvas.md
│   ├── canvas2.md
│   ├── components/
│   └── data/dictionary/
├── scripts/
├── package.json
├── observablehq.config.ts
└── ...
```

**After:**

```
genomic-regions-observable/    # SIBLING — exact copy of current state, untouched
├── src/
├── package.json
└── ...

genomic-regions/                # React app at root, fresh
├── package.json               # Vite/React/TS/Tailwind/DaisyUI deps
├── vite.config.ts
├── tsconfig.json
├── index.html
├── src/
│   ├── main.tsx
│   ├── app.tsx
│   ├── app.css                # Tailwind + DaisyUI plugin imports
│   ├── lib/
│   ├── stores/
│   ├── hooks/
│   ├── components/
│   └── pages/
├── public/
│   └── data/dictionary/       # parquet files (moved from src/data/dictionary/)
└── scripts/                   # keep puppeteer scripts
```

## Tech stack

Locked to bedbase-ui's proven Mosaic+React stack (`/Users/sam/Documents/Work/ai-sandbox/workspaces/sam/bedbase/repos/bedbase-ui`), since they have all the integration patterns working in production.

- **Build**: Vite + TypeScript, `build.target: 'esnext'` (for top-level await needed by DuckDB-WASM init)
- **Vite plugins**: `@vitejs/plugin-react`, `@tailwindcss/vite`, `vite-plugin-wasm`, `vite-plugin-top-level-await`
- **UI**: **React 19** + Tailwind v4 + DaisyUI v5 + `@tailwindcss/typography` + `lucide-react`
- **Routing**: react-router-dom v7
- **State**: **React Context** (not Zustand — bedbase-ui's pattern works fine and avoids an extra dependency). One context per concern: `MosaicCoordinatorContext` (DuckDB-WASM), `SelectionContext` (picked token, brush, lens). `useState` for local UI.
- **UMAP rendering**: **`embedding-atlas/react`** — `EmbeddingViewMosaic` handles WebGL/WebGPU rendering, brush/click, lasso, and Mosaic Selection integration natively. Proven at >1M points; trivial for our 36K. Avoids hand-rolling ref-mounted vgplot for the UMAP itself.
- **Other plots** (Section 1 token bars, chr16 distribution strip): `vg.plot` + Observable Plot mounted via `useEffect` into ref'd divs (the standard bedbase-ui pattern for non-UMAP plots).
- **Mosaic**: `@uwdata/vgplot`, `@uwdata/mosaic-core`, `@uwdata/mosaic-sql` (helpers like `isIn` for predicate building)
- **Data**: DuckDB-WASM via `@uwdata/mosaic-core`'s `wasmConnector`

Styling discipline per `tailwind-ui-styling` skill: utility-first, DaisyUI semantic tokens, no raw color classes, no inline styles, no scoped CSS. Define a custom DaisyUI theme in `app.css` (bedbase-ui has a clean template).

## Scope of the demo

A single primary page with:

1. **Interval picker** (top): the 4 reference featured intervals + the 3 hub candidates (or a clean subset). DaisyUI `select` or button group.
2. **Section 1 — universe + per-file activations** for the picked interval. Tokens-only rendering (universe band + per-file token bars colored by SCREEN class).
3. **Region UMAP** (full chr16, faded background) with the interval's regions colored & clickable on top. Mosaic vgplot with the brush/click bridge from `canvas.md`.
4. **Dictionary card** (right of UMAP). Header + soft class profile + concept axis bars + kNN partners + NPMI partners (under selected lens).
5. **Stratum lens picker** (drives the card's NPMI partners + the orange ring overlay on the UMAP).
6. **chr16 partner-distribution strip** below the UMAP. 250 bins; bar height = universe density; bar color = kNN/NPMI partners per bin (toggle for source); vertical rule at picked token.

Not in scope (defer or drop):

- Module catalogue (Reference page Section 2 expansion table)
- File UMAP / intersection threshold (Reference page Section 3)
- HBA1 static showcase
- Brushable arbitrary-region UMAP (Canvas 2's no-anchor exploration mode)
- Continuous bigwig signal track (Reference page Section 1 alternative mode)

These can be added later if the scoped demo lands cleanly.

## Key technical decisions

### Reactivity model

Observable's reactive cells become explicit React state + memoization:

- Pure cells → `useMemo`
- Side-effect cells (event listeners on Mosaic selections) → `useEffect` hooks with proper cleanup
- `Mutable(value)` + Observable's reactive graph → Zustand store for cross-component shared state
- `view(Inputs.*)` shortcuts → React form components (DaisyUI `select`, `input range`, `radio` group)

### Mosaic integration in React (bedbase-ui pattern)

**For the UMAP**: `EmbeddingViewMosaic` from `embedding-atlas/react` accepts `coordinator`, `table`, `filter` (a `vg.Selection`), and color/size/category props directly. It handles render lifecycle, brush, click, and predicate updates internally. We don't manually mount/cleanup; React handles it as a normal component.

**For other plots** (Section 1 bars, chr16 distribution): mount via `useEffect` into ref'd divs:

```tsx
function ChrStrip() {
  const ref = useRef<HTMLDivElement>(null);
  const { coordinator } = useMosaicCoordinator();
  useEffect(() => {
    if (!ref.current) return;
    const plot = vg.plot(/* marks */);
    ref.current.appendChild(plot);
    return () => { ref.current?.removeChild(plot); };
  }, [coordinator, /* deps */]);
  return <div ref={ref} />;
}
```

Selections (`vg.Selection.intersect`, `vg.intervalXY`, `clausePoint`, `clausePoints`) are framework-agnostic — same code as Observable.

### React Strict Mode (bedbase-ui pattern)

Dev mode double-invokes effects. Mosaic plot setup creates SQL views and event listeners. Defensive patterns proven in bedbase-ui:

- **Coordinator singleton**: `coordinatorRef = useRef<vg.Coordinator | null>(null)` plus a lazy `getCoordinator()` that creates exactly one. Strict Mode's double-invoke of the provider's first useEffect doesn't create two coordinators.
- **Init guards**: `dataInitializedRef = useRef(false)` so re-running `initializeData()` is a no-op if already done.
- **Idempotent SQL**: `CREATE OR REPLACE TABLE/VIEW` everywhere; never `CREATE TABLE` without replacement.
- **Cleanup returns**: every `useEffect` that mounts DOM, adds listeners, or registers selection bridges must return a cleanup. The standard pattern from bedbase-ui's `embedding-plot.tsx` is the reference.

### Methodological cleanup folded in

While building fresh, apply the methodological corrections we identified:

- **Concept axes from BED, not bigwigs.** `embedding_features.py` (pipeline-side) gets rewritten to define case/control via peak counts under cell-line + mark file pools, not bigwig signal means. New `region_concept_axes.parquet` shipped.
- **Bigwigs only used for visualization.** No bigwig means in `viz_chr16` for embedding-side computations. Per-(cell, mark) signal columns can stay only if a continuous signal track is added later (currently out of scope).

## Migration phases

### Phase 0 — sibling backup (15 min)

```bash
cp -r /Users/sam/Documents/Work/genomic-regions/ /Users/sam/Documents/Work/genomic-regions-observable/
```

Verify the backup runs (`cd genomic-regions-observable && npm install && npm run dev`).

### Phase 1 — bootstrap (½ day)

In `genomic-regions/`:

1. Remove Observable-specific files: `observablehq.config.ts`, `src/*.md`, `src/components/`, `src/.observablehq/`, current root `package.json`, `package-lock.json`, `node_modules/`.
2. `npm create vite@latest . -- --template react-ts`
3. Install the bedbase-ui-aligned stack:

   ```
   npm install -D tailwindcss @tailwindcss/vite @tailwindcss/typography daisyui
   npm install -D vite-plugin-wasm vite-plugin-top-level-await
   npm install lucide-react react-router-dom
   npm install @uwdata/vgplot @uwdata/mosaic-core @uwdata/mosaic-sql @observablehq/plot
   npm install embedding-atlas
   ```

4. `vite.config.ts` mirrors bedbase-ui: `react()`, `tailwindcss()`, `wasm()`, `topLevelAwait()`, `build.target: 'esnext'`.
5. `app.css` per skill template + custom DaisyUI theme (steal bedbase-ui's structure, pick our own colors):

   ```css
   @import 'tailwindcss';
   @plugin "@tailwindcss/typography";
   @plugin "daisyui";
   @plugin "daisyui/theme" { name: "regions"; default: true; /* our palette */ }
   ```

6. Move parquets: `mv src/data/dictionary public/data/dictionary`.
7. App shell: DaisyUI nav bar + theme switcher + a "hello world" page that mounts DuckDB-WASM and runs a sanity SELECT.

### Phase 2 — Mosaic integration (½ day)

Crib bedbase-ui's `mosaic-coordinator-context.tsx` exactly — it's the reference implementation.

1. `MosaicCoordinatorContext` provider with the bedbase-ui pattern: `useRef` singleton, `dataInitializedRef` guard, lazy `getCoordinator()`, eager `useEffect(() => initializeData())` to register all parquet tables on mount.
2. `useMosaicCoordinator()` consumer hook.
3. Initial `initializeData()` runs `CREATE OR REPLACE TABLE` for each parquet (regions, files, intervals, etc.) — same shape as Observable's `tablesReady` cell.
4. Sanity test: a component that renders a one-line SELECT result confirms coordinator + tables are wired.

### Phase 3 — data layer (½ day)

1. `lib/duckdb.ts`: coordinator, table registry, parquet URL resolution.
2. `lib/data.ts`: parquet loading hooks (`useRegions`, `useFeaturedFiles`, `useIntervals`); lookup map builders.
3. `lib/colors.ts`: SCREEN class color domain/range, helper functions, ported from Observable's `dictionary.js`.
4. Stub data fetching against the static parquets in `public/data/dictionary/`.

### Phase 4 — interval-scoped Section 1 (1 day)

`<Section1Tokens>` component:

- Universe band (interval's universe regions, colored by SCREEN class)
- Per-file token bars (intersect each featured file's `chr16_active_token_ids` with the interval's universe set)
- DaisyUI styling: cards, semantic colors, tabular nums on count labels.

Validates the data layer + non-Mosaic Plot rendering pattern.

### Phase 5 — Region UMAP + dict card (1 day, less because embedding-atlas)

1. `<RegionUMAP>` using `EmbeddingViewMosaic` from `embedding-atlas/react`. Configure: `coordinator`, `table='regions_classed'`, color by `cclass`, custom highlight layer for the selected interval's regions, click handler updating the picked-token Selection. The brush/click/lasso interactions come for free — no manual ref-mount needed.
2. `<DictCard>` (plain React + DaisyUI cards) with sections: header, soft class profile bars, concept axis bars, kNN partners (chips), NPMI partners under current lens (chips).
3. `<LensPicker>` (DaisyUI select) driving the card's NPMI partners + UMAP overlay highlight.
4. Selection bridges (picked token, kNN partner overlay, NPMI partner overlay) wired through `SelectionContext`.

### Phase 6 — chr16 distribution strip (½ day)

`<ChrStrip>` component:

- 250 bins across chr16, bar height = universe density, bar color = kNN-or-NPMI partner count per bin (continuous YlOrRd scale).
- Source toggle (kNN vs NPMI).
- Vertical rule at picked token's chr16 midpoint.
- Tooltip on hover.

### Phase 7 — polish + parity check (½ day)

- DaisyUI theme switching (light/dark/auto).
- Keyboard navigation where reasonable (arrow keys for interval picker, esc to clear selection).
- Responsive layout: stacks on mobile (`md:` breakpoints).
- Side-by-side screenshot comparison vs Observable canvas page using puppeteer.
- Deployment story: Vercel / Netlify / GitHub Pages.

## Pipeline-side methodological cleanup

In parallel with the React build (or as a follow-up):

1. **Rewrite `embedding_features.py`** to compute concept axes from BED-derived peak counts (not bigwig signal means).
2. **Drop bigwig mean columns from `viz_chr16.parquet`** — saves ~16 MB.
3. **Re-run stages 08 (viz precompute), embedding_features** to regenerate `region_concept_axes.parquet`.

This yields a bigwig-free embedding-derivation pipeline; bigwigs only flow into Section 1 visualization (currently out of scope, but kept available).

## Risks & unknowns

1. **vgplot tooltip / hover lifecycle in React**. Mosaic's DOM event listeners need explicit cleanup. Validate early in phase 2.
2. **React Strict Mode double-rendering**. Idempotency for SQL exec + selection bridges. Plan: every `useEffect` with side effects gets a cleanup return.
3. **TypeScript types for Mosaic API**. Likely incomplete; expect some `as any` early.
4. **Parquet caching in Vite dev server**. May need cache headers tuned for parquet updates.
5. **Tailwind v4 + vgplot styling**. vgplot uses inline styles; should not conflict but verify in phase 2.
6. **DaisyUI v5 freshness**. Latest version has API changes vs v4; cross-check class names against current docs during phase 1.

## Effort estimate

**~4-5 working days** for the scoped demo (phases 0–7), revised down because:
- `embedding-atlas/react` removes most of phase 5's UMAP integration work.
- Cribbing bedbase-ui's `mosaic-coordinator-context.tsx` directly removes most of phase 2's design risk.

**+1-2 days** for the pipeline-side bigwig cleanup.
**+2-3 days** for Canvas 3 features (file UMAP brush + on-demand DuckDB NPMI), if pursued after the scoped demo lands.

## Reference implementations to crib from

- **`bedbase-ui/src/contexts/mosaic-coordinator-context.tsx`** — coordinator singleton, init guards, `CREATE OR REPLACE` patterns, Strict Mode safety.
- **`bedbase-ui/src/components/umap/embedding-plot.tsx`** — `EmbeddingViewMosaic` configuration, Selection composition, brush/click handling, custom point highlighting.
- **`bedbase-ui/vite.config.ts`** — Vite plugin set + `build.target` for top-level await.
- **`bedbase-ui/src/app.css`** — Tailwind v4 + DaisyUI v5 + custom theme template.

## Plan comprehension quiz topics (for next session)

When implementation starts, re-quiz on:

1. **Why are we scoping the demo to curated intervals instead of building the full chr16 brush flow first?** *Expected*: focused storytelling + simpler UX surface; the full-chr16 brush (canvas2's no-anchor mode) is more powerful but harder to demo cleanly. Add it after the scoped demo proves the React+Mosaic stack works.
2. **What's the riskiest part of the build?** *Expected*: Mosaic interactivity in React's lifecycle (Strict Mode + cleanup). Phase 2 validates this on a single component before scaling up.
3. **Why fold the bigwig-removal cleanup in alongside the React build instead of doing one then the other?** *Expected*: the React app needs to consume the new BED-derived concept axes anyway — sequencing them together avoids shipping a temporary bigwig-dependent React build that we'd then have to retrofit. Either order works; the integration matters.
