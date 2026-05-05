# scratch

Two sandboxes for getting hands-on with the file types and operations that the
planned experiments depend on. Each subdir is self-contained.

- **[chr22_demo/](chr22_demo/)** — the original walkthrough. K562 ATAC /
  H3K27ac / H3K4me1 bigwigs sliced to chr22, plus the cCRE BED for the same
  chromosome, and example scripts for `pyBigWig`, `deeptools`, and `gtars`.
  Useful for learning bigwig-plus-BED mechanics at a portable scale.

- **[ccre_examples/](ccre_examples/)** — curated playground for the planned
  pilots. Selects 30 cCREs (10 PLS / 10 pELS / 10 dELS) stratified by
  tissue-specificity and CRISPRi-validation status, then pulls bigwigs and
  peak BEDs for six modalities across K562, GM12878, HepG2. Each cCRE gets a
  track-stack plot so you can eyeball continuous vs. binary signal inside a
  single interval before scaling up.

Both share the workspace root's `uv` environment (`~/Documents/Work/spatial-region-features/pyproject.toml`).
