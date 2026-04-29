# genomic-dict

Interactive "regulatory genomics dictionary" — course data-viz project.

See [`../plans/2026-04-22-regulatory-dictionary-viz.md`](../plans/2026-04-22-regulatory-dictionary-viz.md) for the plan.

## Layout

```
genomic-dict/
├── data/           # gitignored — corpus BEDs, universe files, annotations, precomputed viz artifacts
├── pipeline/       # numbered Python stages (00_curate_corpus.py, 01_prepare_universe.py, ...)
├── models/         # trained R2V checkpoints
└── notebook/       # Observable source or link
```

Dependencies are managed at the workspace root (`../pyproject.toml`).
