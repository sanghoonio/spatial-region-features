# AI Usability Guidelines: genomic-dict pipeline

**Date:** 2026-04-22
**Project type:** Multi-stage data pipeline (numbered Python scripts under `genomic-dict/pipeline/`)
**Stack:** Python 3.13, uv-managed env, invoked via `python genomic-dict/pipeline/NN_name.py ...` or `uv run`
**Consumers:** AI agent (Claude) orchestrating stages via `Bash` + human (Sam) reviewing outputs between stages

Not a published library, not a packaged CLI — but the scripts have CLI-like properties and must be AI-orchestratable. Library principles mostly don't apply; CLI principles largely do.

## User-stated requirements (load-bearing)

1. **Numbered scripts** — each pipeline stage is `NN_<verb>.py` where `NN` is a two-digit order.
2. **Every script emits an AI-ingestible JSON summary** — machine-readable record of what happened, written next to the stage's outputs.
3. **Central config file for params** — one source of truth for pipeline parameters; scripts read it, optionally overridden by CLI flags.

These three are load-bearing because they're what lets the AI orchestrate the pipeline without re-reading the code each time.

---

## Applicable principles (concrete for this project)

### P1 — Structured output (JSON summary per stage)

**Why it matters:** between stages, Claude must read what happened and decide what to do next. A printed log is unparseable; a JSON summary is a contract.

**Concrete requirements:**
- Every stage writes `<output_dir>/<NN>_<name>.summary.json` alongside its primary outputs.
- Summary schema (minimum):
  ```json
  {
    "stage": "01_curate_corpus",
    "status": "success",
    "timestamp": "2026-04-22T14:30:00Z",
    "duration_seconds": 123.4,
    "inputs": [{"path": "...", "sha256": "...", "size_bytes": 12345}],
    "outputs": [{"path": "...", "sha256": "...", "size_bytes": 12345, "record_count": 20000}],
    "config_used": { /* resolved config section for this stage */ },
    "config_file": "genomic-dict/config.toml",
    "config_hash": "...",
    "warnings": [],
    "metrics": { /* stage-specific numerical results */ },
    "versions": {"python": "3.13.1", "genomic_dict": "0.1.0", "geniml": "0.8.4"}
  }
  ```
- Primary data outputs go to files (parquet / BED / joblib), not stdout.
- Scripts may also print a short human-readable summary to stderr (not stdout), gated by `--quiet`.

**Example:**
```python
summary = {
    "stage": "01_curate_corpus",
    "status": "success",
    ...,
    "metrics": {"input_files": 50384, "output_files": 20000, "filter_rejection_rate": 0.603},
}
Path(out_dir / "01_curate_corpus.summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
```

### P2 — Config interchangeability (central TOML, CLI overrides)

**Why it matters:** one source of truth prevents drift between what Claude reports and what scripts actually ran. CLI flag overrides let Claude run single-stage reruns without editing the config file.

**Concrete requirements:**
- Single central config at `genomic-dict/config.toml`.
- Structure: `[paths]` for shared paths, `[stage_NN_name]` section per stage with that stage's params.
- Every script loads `config.toml`, reads its own `[stage_NN_name]` section, then applies CLI flag overrides.
- Key names match between config and CLI: config `min_region_count` ↔ CLI `--min-region-count`.
- `--dump-config` flag on every script prints the full resolved config (post-override) as JSON to stdout.
- **Config file is the methods section.** Check it in to git; treat edits as scientific decisions.

**Example config:**
```toml
[paths]
data_dir = "genomic-dict/data"
models_dir = "genomic-dict/models"

[stage_01_curate_corpus]
min_tissue_label = true
min_assay_label = true
target_file_count = 20000
seed = 42

[stage_05_train_r2v]
vocab_source = "screen_ccre"
embedding_dim = 100
window_size = 5
epochs = 50
seed = 42
```

### P3 — Composable steps (numbered stages, runnable independently)

**Why it matters:** this is already your stated requirement — numbered scripts. The principle sharpens it: each stage must be **independently re-runnable** given only upstream outputs, not only after a full pipeline run from stage 00.

**Concrete requirements:**
- Stage `NN` reads its inputs from filesystem paths declared in config (not from previous stage's in-memory state).
- Stage `NN` writes its outputs to filesystem paths declared in config.
- Running `python NN_foo.py` in isolation works, assuming upstream output files exist.
- No "run all" monolithic wrapper for Shape A — Claude invokes stages individually and inspects summaries between them. A thin `run_pipeline.sh` (or Snakefile, later) is fine once the contract stabilizes, but not required and definitely not the primary interface.

### P4 — Determinism by default

**Why it matters:** the R2V training and UMAP precompute are stochastic. Without fixed seeds, rerunning a stage gives different outputs, which defeats stage-level caching and makes Panel 3's diagnostic numbers irreproducible.

**Concrete requirements:**
- Every stochastic stage reads `seed` from its config section (default: `42`).
- Sort outputs before writing (e.g., BED files sorted by coord; token lists sorted by token id).
- No timestamp-suffixed filenames. Timestamps go in the summary JSON's `timestamp` field, not in paths.
- UMAP: fix `random_state`. R2V: fix seed in `Region2VecExModel` train call.

### P5 — Machine-readable errors + exit codes

**Why it matters:** when a stage fails, Claude has to decide whether to retry, change a param, or escalate. An exit code + structured error JSON says this clearly.

**Concrete requirements:**
- Exit codes:
  - `0` — success
  - `1` — bad input (upstream file missing, malformed, empty)
  - `2` — bad config / params (validation failure)
  - `3` — runtime failure (network, disk, unexpected exception)
- On failure, write a partial summary with `"status": "error"` and an `"error": {...}` block:
  ```json
  {
    "status": "error",
    "error": {
      "class": "bad_input",
      "message": "Corpus file 01_curate_corpus/corpus.parquet not found",
      "suggestion": "Run stage 01 first, or set [paths].corpus_file to an existing file.",
      "retryable": false
    }
  }
  ```
- Exceptions propagate with tracebacks to stderr; structured JSON goes to stdout or the summary file.

### P7 — Self-documenting parameters

**Why it matters:** Claude should be able to invoke a stage cold without reading its source. `--help` is the script's contract.

**Concrete requirements:**
- `argparse` or `click` with explicit `type=`, `default=`, and helpful `help=` on every flag.
- `--help` lists every CLI flag with its type, default, and one-line purpose.
- `--dump-config` prints the resolved config as JSON (see P2).
- Script module docstring at the top explains: what this stage does, what it reads, what it writes, typical runtime.

### P9 — Provenance in every summary

**Why it matters:** when Claude is reading back summaries to decide next moves, provenance tells Claude whether a summary is from *this* config or a stale earlier run.

**Concrete requirements:**
- Every summary includes: `timestamp`, `config_file`, `config_hash`, `versions` of key deps, `git_commit` of the workspace (if clean).
- Input and output file entries include `sha256` so Claude can detect if an upstream file changed between its summary and this stage's input.

---

## Principles that do NOT apply (and why)

- **P6 (no interactive prompts)** — scripts are non-interactive by construction; no affordance for prompts exists.
- **P8 (pipe-friendly I/O)** — this is a multi-stage batch pipeline, not a filter-chain tool. Stages read and write files via config paths; stdin/stdout streaming is not the interaction model.
- **Library principles (S1–S10)** — this project is not a library; there is no public API. Skip.

**Hygiene items worth stealing from library principles anyway:**
- **Type hints (S1)** — annotate public functions in `pipeline/_common.py` helpers.
- **Actionable error messages (S4)** — the P5 error JSON's `"suggestion"` field is the concrete form of this.

---

## Shared infrastructure (implied by the guidelines)

The guidelines imply a small shared helper module. I'll propose it as a design call for you to approve before writing:

**`genomic-dict/pipeline/_common.py`** — shared utilities:
- `load_config(path: str, stage: str) -> dict` — load TOML, return `[paths]` merged with `[stage_<NAME>]`
- `parse_cli_overrides(parser, config_section: dict)` — CLI flags override config keys
- `write_summary(out_dir: Path, stage: str, status: str, ...) -> None` — standard summary JSON writer with auto-provenance
- `sha256_file(path: Path) -> str`
- `versions_dict() -> dict` — collect key dep versions for provenance

Every numbered script imports from `_common` so the contract is enforced by shared code rather than by convention.

---

## Pre-flight checklist (apply before running any pipeline stage the first time)

- [ ] Script reads config from `genomic-dict/config.toml` and respects `--dump-config`
- [ ] CLI flags follow kebab-case and match snake_case config keys
- [ ] `--help` lists every flag with type/default/help
- [ ] Script writes `<NN>_<name>.summary.json` next to its outputs
- [ ] Summary includes inputs (with sha256), outputs (with sha256), config_used, timestamp, versions
- [ ] Stochastic operations use a config-supplied `seed`
- [ ] Outputs are sorted deterministically
- [ ] On failure, exit code is 1/2/3 (not `sys.exit(1)` for everything) and a partial summary with `"status": "error"` is written
- [ ] No output filename depends on a timestamp or random value

---

## Resolved decisions (as of 2026-04-22)

- **Config format: YAML** (`genomic-dict/config.yaml`) — consistency with `aging/` and `copd/` analyses.
- **No CLI framework.** Scripts take exactly one arg, `--config <path>`, via stdlib `argparse`. All other parameters flow from the YAML. If you want to try different params, edit the YAML and rerun.
- **Shared `_common.py` at `pipeline/scripts/_common.py`** provides `stage_start(...)` (loads config, resolves stage section, creates results dir, announces to stderr) and `write_summary(...)` (emits JSON with inputs/outputs + sha256, config_used, versions, git_commit, duration).
- **Cluster target:** Rivanna (`sp5fd@login.hpc.virginia.edu:/home/sp5fd/spatial-region-features`), allocation `shefflab`. SLURM wrappers in `pipeline/slurm/` only for cluster-bound stages (currently planned: stage 05 R2V training).
- **Directory layout** (mirrors `aging/`):
  ```
  genomic-dict/
  ├── config.yaml                  # methods section; git-tracked
  ├── pipeline/
  │   ├── scripts/                 # NN_*.py + _common.py
  │   └── slurm/                   # sbatch wrappers for cluster stages
  ├── data/                        # gitignored — corpus, universe, annotations, hf cache
  ├── results/                     # per-stage NN_*/ with summary.json + small artifacts (git-tracked)
  ├── models/                      # trained checkpoints (gitignored except .gitkeep)
  └── logs/                        # SLURM logs (gitignored)
  ```
- **Sync/pull scripts** at workspace root: `sync.sh` (push code to Rivanna; excludes `.venv`, `data/`, model binaries) and `pull-results.sh` (pull `*.json`/`*.md`/`*.html` from `results/**` + `models/` + SLURM logs back).
