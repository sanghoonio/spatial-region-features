# Observable notebook scaffold — Regulatory Dictionary (5-step narrative)

This scaffold organizes the notebook around the 5-step narrative arc: BED files →
tokenization → embedding → file/region cross-link → interpretation.

**Before pasting cells**: create a notebook on observablehq.com and upload these
parquet files (paperclip icon → upload):

- `viz_chr16.parquet` — region dictionary (35,934 chr16 tokens; UMAP, kNN, marks)
- `viz_files.parquet` — experiment dictionary (84,718 files; file UMAP)
- `featured_intervals.parquet` — 4 narrative intervals with universe slices
- `featured_files.parquet` — ~21 curated files with chr16 token activations
- `featured_tracks.parquet` — per-(file, interval) peak coords and active tokens

Sections below are pasted in order. Markdown cells go in Markdown cells; JS in regular cells.

---

## Section 1 — Setup

### Cell 1 — title (Markdown)

```md
# A Dictionary of Regulatory Genomics

Most of the human genome is non-coding. Most genetic variants linked to disease
sit in non-coding DNA. To interpret what those variants mean, we need a way to
read what regulatory machinery is at any genomic location.

This notebook walks through one approach: treat regulatory elements as **words**
in a vocabulary, treat experimental measurements (BED files) as **documents**,
and use the same kind of distributional learning that works for natural
language to build a *dictionary of regulatory genomics*. We focus on chromosome
16 of the human genome, with four narrative-anchoring loci.

The five steps:

1. **Here are some BED files. What's in them?**
2. **We tokenize them against a fixed universe.**
3. **The universe organizes itself by co-occurrence — embedding space.**
4. **Each experiment uses a slice of the dictionary.**
5. **Hypothesizing about what we don't know.**
```

### Cell 2 — DuckDB handle

```javascript
db = DuckDBClient.of({
  regions: FileAttachment("viz_chr16.parquet"),
  files: FileAttachment("viz_files.parquet"),
  featured_intervals: FileAttachment("featured_intervals.parquet"),
  featured_files: FileAttachment("featured_files.parquet"),
  featured_tracks: FileAttachment("featured_tracks.parquet"),
})
```

### Cell 3 — region table

```javascript
regions = await db.query(`
  SELECT
    token_id, region, chrom, "start", "end",
    overlaps_screen, accession_hex, cclass,
    umap_x, umap_y,
    knn_token_ids, knn_distances,
    "K562__ATAC__mean" AS k562_atac_mean,
    "GM12878__ATAC__mean" AS gm12878_atac_mean,
    "HepG2__ATAC__mean" AS hepg2_atac_mean,
    "K562__H3K4me1__mean" AS k562_h3k4me1,
    "K562__H3K4me3__mean" AS k562_h3k4me3,
    "K562__H3K27ac__mean" AS k562_h3k27ac,
    "K562__H3K27me3__mean" AS k562_h3k27me3,
    "K562__H3K9me3__mean" AS k562_h3k9me3,
    "GM12878__H3K4me1__mean" AS gm_h3k4me1,
    "GM12878__H3K4me3__mean" AS gm_h3k4me3,
    "GM12878__H3K27ac__mean" AS gm_h3k27ac,
    "GM12878__H3K27me3__mean" AS gm_h3k27me3,
    "GM12878__H3K9me3__mean" AS gm_h3k9me3,
    "HepG2__H3K4me1__mean" AS hepg2_h3k4me1,
    "HepG2__H3K4me3__mean" AS hepg2_h3k4me3,
    "HepG2__H3K27ac__mean" AS hepg2_h3k27ac,
    "HepG2__H3K27me3__mean" AS hepg2_h3k27me3,
    "HepG2__H3K9me3__mean" AS hepg2_h3k9me3
  FROM regions
`)
```

### Cell 4 — files, intervals, featured tables

```javascript
files = await db.query(`SELECT * FROM files`)
```

```javascript
featuredIntervals = await db.query(`SELECT * FROM featured_intervals`)
```

```javascript
featuredFiles = await db.query(`SELECT * FROM featured_files`)
```

```javascript
featuredTracks = await db.query(`SELECT * FROM featured_tracks`)
```

### Cell 5 — color palette + helpers

```javascript
classColors = ({
  "PLS":         "#ff0000",
  "pELS":        "#ffa700",
  "dELS":        "#ffcd00",
  "CA-CTCF":     "#00b0f0",
  "CA-H3K4me3":  "#ffaaaa",
  null:          "#cccccc",
})
```

```javascript
assayColors = ({
  "ATAC-seq":          "#1f77b4",
  "DNase-seq":         "#2ca02c",
  "ChIP-seq":          "#9467bd",
  "TF ChIP-seq":       "#e377c2",
  "Histone ChIP-seq":  "#ff7f0e",
})
```

```javascript
// Token id → row lookup, used by step 2 and step 4 to resolve coords.
tokenLookup = new Map(regions.map(r => [r.token_id, r]))
```

### Cell 6 — shared mutable state

```javascript
mutable currentIntervalId = "fto_irx3"
```

```javascript
mutable highlightedFileId = null
```

```javascript
mutable selectedToken = null
```

### Cell 7 — interval picker (viewof)

```javascript
viewof currentIntervalIdInput = Inputs.select(
  featuredIntervals.map(i => i.interval_id),
  {
    label: "Featured interval",
    value: currentIntervalId,
    format: id => featuredIntervals.find(i => i.interval_id === id).label,
  }
)
```

```javascript
{
  // sync the input back to the mutable
  mutable currentIntervalId = currentIntervalIdInput;
}
```

```javascript
currentInterval = featuredIntervals.find(i => i.interval_id === currentIntervalId)
```

---

## Section 2 — Step 1: What's in a BED file?

### Cell 8 — narrative (Markdown)

```md
## Step 1 — Reading the raw experimental data

A BED file is the standard format for genomic experimental results: a list of
intervals where some kind of activity was detected. ATAC-seq files contain
intervals where chromatin is accessible. ChIP-seq files contain intervals where
a specific protein binds. DNase-seq files contain intervals of nuclease
hypersensitivity.

Below is the **${currentInterval.label}** locus on chromosome 16. The narrative
hook for this region: *${currentInterval.narrative_caption}*.

Each row is one ENCODE experiment. Rectangles mark places that experiment
called peaks. Notice how different experiments emphasize different regions —
they're measuring different molecular events.
```

### Cell 9 — interval coordinate display (Markdown helper)

```md
**Interval**: ${currentInterval.chrom}:${currentInterval.start.toLocaleString()}–${currentInterval.end.toLocaleString()} (${((currentInterval.end - currentInterval.start) / 1000).toFixed(1)} kb)
```

### Cell 10 — peak tracks (Step 1 visual)

```javascript
peakTracksData = {
  // Flatten featured_tracks into one row per peak
  const rows = []
  for (const t of featuredTracks) {
    if (t.interval_id !== currentIntervalId) continue
    const file = featuredFiles.find(f => f.file_id === t.file_id)
    if (!file) continue
    for (let i = 0; i < t.peak_starts.length; i++) {
      rows.push({
        file_id: t.file_id,
        file_label: `${file.cell_line || "—"} · ${file.target || file.assay}`,
        role: file.role,
        assay: file.assay,
        peak_start: t.peak_starts[i],
        peak_end: t.peak_ends[i],
      })
    }
  }
  return rows
}
```

```javascript
peakTracksPlot = Plot.plot({
  width: 900,
  height: Math.max(300, 22 * new Set(peakTracksData.map(d => d.file_label)).size),
  marginLeft: 200,
  x: {
    domain: [currentInterval.start, currentInterval.end],
    label: currentInterval.chrom,
    grid: true,
    tickFormat: d => (d / 1e6).toFixed(2) + "M",
  },
  y: { label: "experiment", domain: [...new Set(peakTracksData.map(d => d.file_label))] },
  color: {
    domain: Object.keys(assayColors),
    range: Object.values(assayColors),
    legend: true,
  },
  marks: [
    Plot.rect(peakTracksData, {
      x1: "peak_start",
      x2: "peak_end",
      y: "file_label",
      fill: "assay",
      stroke: d => d.role === "mystery" ? "black" : null,
      strokeWidth: 1,
      title: d => `${d.file_label}\npeak: ${d.peak_start}-${d.peak_end}`,
    }),
  ],
})
```

**TODO** for the user:
- Decide whether to show all 21 files or filter to a more digestible 5–10 per
  interval.
- Add tissue grouping if files-per-tissue is uneven.
- Style the mystery files (currently outlined in black) however you want them
  to stand out.

---

## Section 3 — Step 2: Tokenizing against a universe

### Cell 11 — narrative (Markdown)

```md
## Step 2 — Projecting onto a fixed vocabulary

The peaks above are continuous coordinates — they don't share boundaries
between experiments, and the same regulatory region might have slightly
different peak boundaries in different files. To build a *dictionary*, we need
a fixed vocabulary.

We project each file onto a **universe** of regulatory regions — the SCREEN
candidate cis-regulatory elements (cCREs), which integrates evidence across
many ENCODE experiments. Each file becomes a **binary presence vector**: for
each universe region, did this experiment overlap it?

Now look at the same locus, but with the experimental peaks (faded) and the
universe regions overlaid. Each row shows where a given file *activates* the
universe. The same data — but now in a shared vocabulary.
```

### Cell 12 — universe + binary tokenization plot (Step 2 visual)

```javascript
binaryTokensData = {
  // For each (file, interval) record, expand active tokens into per-token rows
  const rows = []
  for (const t of featuredTracks) {
    if (t.interval_id !== currentIntervalId) continue
    const file = featuredFiles.find(f => f.file_id === t.file_id)
    if (!file) continue
    for (const tid of t.active_token_ids) {
      const region = tokenLookup.get(tid)
      if (!region) continue
      rows.push({
        file_id: t.file_id,
        file_label: `${file.cell_line || "—"} · ${file.target || file.assay}`,
        token_id: tid,
        cclass: region.cclass,
        token_mid: (region.start + region.end) / 2,
      })
    }
  }
  return rows
}
```

```javascript
universeRowData = {
  // Universe tokens within the current interval, for the top "universe" track
  const interval = currentInterval
  return regions
    .filter(r =>
      r.chrom === interval.chrom &&
      r.end > interval.start &&
      r.start < interval.end
    )
    .map(r => ({
      token_id: r.token_id,
      cclass: r.cclass,
      mid: (r.start + r.end) / 2,
    }))
}
```

```javascript
binaryTokensPlot = {
  // Two stacked panels: universe band on top, file rows below.
  const fileLabels = [...new Set(binaryTokensData.map(d => d.file_label))]
  return html`
    <div>
      ${Plot.plot({
        width: 900,
        height: 60,
        marginLeft: 200,
        x: {
          domain: [currentInterval.start, currentInterval.end],
          axis: null,
        },
        y: { domain: ["universe"], label: null },
        color: {
          domain: Object.keys(classColors).filter(k => k !== "null"),
          range: Object.keys(classColors).filter(k => k !== "null").map(k => classColors[k]),
          legend: true,
        },
        marks: [
          Plot.tickX(universeRowData, {
            x: "mid",
            y: () => "universe",
            stroke: "cclass",
            strokeWidth: 1.5,
            title: d => `cCRE class: ${d.cclass || "(none)"}`,
          }),
        ],
      })}
      ${Plot.plot({
        width: 900,
        height: Math.max(220, 22 * fileLabels.length),
        marginLeft: 200,
        x: {
          domain: [currentInterval.start, currentInterval.end],
          label: currentInterval.chrom,
          grid: true,
          tickFormat: d => (d / 1e6).toFixed(2) + "M",
        },
        y: { label: "experiment", domain: fileLabels },
        color: { type: "categorical", legend: false },
        marks: [
          Plot.dot(binaryTokensData, {
            x: "token_mid",
            y: "file_label",
            fill: "cclass",
            r: 2.2,
            opacity: 0.85,
            title: d => `token ${d.token_id} (${d.cclass || "no class"})`,
          }),
        ],
      })}
    </div>
  `
}
```

**TODO**:
- Currently shows two stacked plots (universe band + file rows). User can rework
  into a single plot if cleaner.
- Consider showing the original peaks faded as a backdrop, with tokens overlaid
  as dots — the visual transition from continuous to discrete is the pedagogical
  point.
- Color tokens by SCREEN class (current default), or alternatively by file's
  assay if step 2 is about "the experiment's activation."

---

## Section 4 — Step 3: From local to global — embedding space

### Cell 13 — narrative (Markdown)

```md
## Step 3 — A learned grammar of the genome

So far we've looked at one locus. But the universe has 35,934 regions on
chromosome 16 alone — and across the whole genome, ~1 million. With 84,000+
BED files in the corpus, every region has a co-occurrence pattern: which
*other* regions does it tend to be active alongside?

A natural-language-style model (Region2Vec, derived from word2vec) takes that
co-occurrence and builds a 100-dimensional embedding for each region. Regions
that show up in similar files together end up close in embedding space.

Project to 2D, color by SCREEN regulatory class, and you can see the embedding
*has learned the dictionary*: promoters cluster together, enhancers cluster
together, CTCF-bound insulators cluster together. The grammar is emergent.
```

### Cell 14 — region UMAP (Step 3 visual)

```javascript
regionsData = regions.filter(d => d.cclass !== null)  // hide unclassed by default
```

```javascript
regionsPlot = Plot.plot({
  width: 900,
  height: 600,
  color: {
    domain: Object.keys(classColors).filter(k => k !== "null"),
    range: Object.keys(classColors).filter(k => k !== "null").map(k => classColors[k]),
    legend: true,
  },
  marks: [
    Plot.dot(regionsData, {
      x: "umap_x", y: "umap_y",
      fill: "cclass",
      r: 1.2,
      opacity: 0.6,
      title: d => `${d.region}\n${d.cclass}`,
    }),
    // Highlight the tokens within the current interval, to tie back to step 2
    Plot.dot(regionsData.filter(d => {
      return d.chrom === currentInterval.chrom &&
             d.end > currentInterval.start &&
             d.start < currentInterval.end
    }), {
      x: "umap_x", y: "umap_y",
      stroke: "black",
      r: 2,
      strokeWidth: 1,
    }),
  ],
})
```

**TODO**:
- The current interval's tokens are outlined to tie steps 2 and 3 visually. User
  can change this to a different visual encoding.
- Could add `Inputs.toggle` to switch between class-coloring vs. tissue-activity
  coloring (which would re-emphasize a different facet of the embedding).
- Consider hover-and-pin interaction so users can keep tooltip context.

---

## Section 5 — Step 4: Each experiment is a partial vocabulary

### Cell 15 — narrative (Markdown)

```md
## Step 4 — Cross-linking experiments and regions

Every BED file activates a *subset* of the universe. So an experiment is a
"sentence" that uses some words from the dictionary and not others.

The two scatters below show this. On the left: the **experiment dictionary** —
84,000 BED files, positioned by their file-level R2V embedding (their overall
"vocabulary fingerprint"). On the right: the **region dictionary** from Step 3.

Click a file in the experiment dictionary, and the regions it activates light
up in the region dictionary. You can see that experiments aren't random
samples — they cluster, and the regions they activate cluster too.
```

### Cell 16 — file UMAP (left side of step 4)

```javascript
viewof filterAssay = Inputs.select(
  ["(all)", "DNase-seq", "ATAC-seq", "ChIP-seq", "TF ChIP-seq", "Histone ChIP-seq"],
  { label: "Filter by assay", value: "(all)" }
)
```

```javascript
filesPlot = {
  const data = files.filter(d => filterAssay === "(all)" || d.assay === filterAssay)
  const plot = Plot.plot({
    width: 600,
    height: 500,
    color: {
      domain: Object.keys(assayColors),
      range: Object.values(assayColors),
      legend: true,
    },
    marks: [
      Plot.dot(data.filter(d => !d.is_unlabeled), {
        x: "umap_x", y: "umap_y",
        fill: "assay",
        r: 1.5,
        opacity: 0.5,
        title: d => `${d.assay}\n${d.cell_line}\n${d.name}`,
      }),
      Plot.dot(data.filter(d => d.is_unlabeled), {
        x: "umap_x", y: "umap_y",
        stroke: "black",
        fill: "white",
        r: 5,
        symbol: "triangle",
        title: d => `[mystery] ${d.name}`,
      }),
      // Highlight currently selected (featured) file
      Plot.dot(data.filter(d => d.id === highlightedFileId), {
        x: "umap_x", y: "umap_y",
        stroke: "magenta",
        strokeWidth: 2,
        r: 8,
        fill: "none",
      }),
    ],
  })
  return plot
}
```

### Cell 17 — featured-file picker (Step 4 click stand-in)

```javascript
// The "click on file UMAP" interaction is hard in Plot. As a stand-in, expose
// a dropdown of the featured + mystery files. Wiring real clicks is a TODO.
viewof pickedFileId = Inputs.select(
  [null, ...featuredFiles.map(f => f.file_id)],
  {
    label: "Highlight a featured experiment",
    format: id => {
      if (id === null) return "(none)"
      const f = featuredFiles.find(x => x.file_id === id)
      return `${f.role === "mystery" ? "[mystery] " : ""}${f.cell_line || "—"} · ${f.target || f.assay}`
    },
  }
)
```

```javascript
{ mutable highlightedFileId = pickedFileId; }
```

### Cell 18 — region UMAP with file-highlight overlay (right side of step 4)

```javascript
regionsPlotWithHighlight = {
  const baseData = regions.filter(d => d.cclass !== null)
  const featuredFile = featuredFiles.find(f => f.file_id === highlightedFileId)
  const activeSet = new Set(featuredFile ? featuredFile.chr16_active_token_ids : [])
  const inactive = baseData.filter(d => !activeSet.has(d.token_id))
  const active = baseData.filter(d => activeSet.has(d.token_id))

  return Plot.plot({
    width: 600,
    height: 500,
    color: {
      domain: Object.keys(classColors).filter(k => k !== "null"),
      range: Object.keys(classColors).filter(k => k !== "null").map(k => classColors[k]),
      legend: true,
    },
    marks: [
      Plot.dot(inactive, {
        x: "umap_x", y: "umap_y",
        fill: "cclass",
        r: 1,
        opacity: highlightedFileId ? 0.1 : 0.6,
      }),
      Plot.dot(active, {
        x: "umap_x", y: "umap_y",
        fill: "cclass",
        stroke: "black",
        strokeWidth: 0.4,
        r: 2.5,
        opacity: 0.95,
      }),
    ],
    title: featuredFile
      ? `${featuredFile.role === "mystery" ? "[mystery] " : ""}${featuredFile.cell_line || "—"} · ${featuredFile.target || featuredFile.assay} — ${activeSet.size} active tokens`
      : "Pick a file to highlight its activated regions",
  })
}
```

**TODO**:
- Real click-to-select on `filesPlot` (Plot's `Plot.pointer` is hover-only;
  click needs custom DOM event handlers — there are recipes online).
- Side-by-side layout for cells 16 and 18.
- Consider a "diff against another experiment" mode: pick two files, show
  the symmetric difference of activated tokens.

---

## Section 6 — Step 5: Hypothesizing about the unknown

### Cell 19 — narrative (Markdown)

```md
## Step 5 — Interpretation beyond what's labeled

The dictionary's class labels (PLS, pELS, dELS, ...) come from SCREEN's
integrative ENCODE analysis. They cover a lot but they're coarse — and they
don't tell us, for any given region, *what biological context this region
operates in*.

Two open-ended questions the dictionary helps us ask:

1. **For an unlabeled region**: what are its embedding-space neighbors? What
   do *they* look like? Can we build a hypothesis from neighborhood structure?
2. **For a mystery experiment**: which experiments live near it in
   experiment-embedding space? What kind of measurement is it likely to be?

This is the dictionary as a *generative tool* — not just lookup, but
interpolation between known things to reason about unknown ones.
```

### Cell 20 — region hypothesis generator (placeholder)

```javascript
hypothesisForToken = {
  if (!selectedToken) {
    return html`<em>Click or pick a region to hypothesize about its function.</em>`
  }

  const ids = selectedToken.knn_token_ids || []
  const neighbors = ids.map(tid => tokenLookup.get(tid)).filter(Boolean)

  const classCounts = neighbors.reduce((acc, n) => {
    const k = n.cclass ?? "(none)"
    acc[k] = (acc[k] ?? 0) + 1
    return acc
  }, {})

  const dominant = Object.entries(classCounts).sort((a, b) => b[1] - a[1])[0]

  return html`
    <div style="padding:1em; border:1px solid #ddd; border-radius:6px">
      <h4>Hypothesis for ${selectedToken.region}</h4>
      <p>Of the 30 nearest neighbors in embedding space:</p>
      <ul>${Object.entries(classCounts).map(([k, n]) =>
        html`<li><strong>${k}</strong>: ${n}</li>`
      )}</ul>
      <p><em>Best-guess class: <strong>${dominant?.[0]}</strong>
      (${Math.round(100 * dominant[1] / neighbors.length)}% of neighbors).</em></p>
      <!-- TODO: enrich with mark-profile averaging, dominant tissue, motif overlap -->
    </div>
  `
}
```

### Cell 21 — file hypothesis generator (placeholder)

```javascript
hypothesisForMystery = {
  if (!highlightedFileId) {
    return html`<em>Highlight a mystery file to hypothesize about it.</em>`
  }
  // TODO: compute file-UMAP kNN among labeled files,
  // aggregate dominant assay / cell_line / tissue,
  // surface as a templated paragraph.
  return html`
    <div style="padding:1em; border:1px dashed #999">
      <em>Hypothesis-generation cell — TODO. See
      <code>plans/2026-04-26-region-interpretation-step5.md</code>.</em>
    </div>
  `
}
```

---

## Suggested order of attack (notebook author POV)

1. Paste cells 1–7. Confirm interval picker shows 4 entries.
2. Paste step 1 (cells 8–10). The peak tracks should render for each interval as you swap.
3. Paste step 2 (cells 11–12). Universe band + binary token rows. This is the most visually-dense step; budget time to refine.
4. Paste step 3 (cells 13–14). Region UMAP with current-interval tokens outlined.
5. Paste step 4 (cells 15–18). File picker drives the highlight on the region UMAP.
6. Paste step 5 (cells 19–21) as placeholders for now.
7. Polish: shared layout, consistent typography, scroll structure, intro/outro narrative.

## Things I left as deliberate placeholders for you

- The actual click-to-select mechanic on the file UMAP. Plot's recipe is to
  attach a DOM event listener to the returned SVG; needs hand-rolled hit testing.
- Coordinated zoom / pan between scatters.
- Style of the "mystery" files vs labeled files.
- Whether to show step 1's peak tracks and step 2's binary tracks side-by-side
  vs. stacked vs. as a toggle.
- Featured-entry deep dive at FTO/IRX3 — would replace the generic interval
  picker for that one with a curated narration.
- Visual transitions between steps (Observable supports CSS, scrollytelling
  via observablehq/scrolly).
