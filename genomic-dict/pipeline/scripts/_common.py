"""Shared utilities for genomic-dict pipeline stages.

Every pipeline stage script should follow this pattern:

    from _common import stage_start, write_summary, file_record, PROJECT_ROOT

    def main() -> None:
        ctx = stage_start("00_inspect_metadata", __doc__)
        try:
            # ... do work ...
            write_summary(ctx, "success", outputs=[...], metrics={...})
        except Exception as e:
            write_summary(ctx, "error", error={
                "class": type(e).__name__,
                "message": str(e),
                "retryable": True,
            })
            raise

    if __name__ == "__main__":
        main()

This module encapsulates the contract: argparse --config, YAML config load,
resolved stage section, results dir creation, summary JSON emission with
provenance (inputs/outputs sha256, config hash, versions, git commit).
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

# Directory constants ---------------------------------------------------------
# _common.py lives at genomic-dict/pipeline/scripts/_common.py.
PROJECT_ROOT = Path(__file__).resolve().parents[2]      # genomic-dict/
REPO_ROOT = Path(__file__).resolve().parents[3]         # spatial-region-features/
DEFAULT_CONFIG = PROJECT_ROOT / "config.yaml"

# Packages whose versions we record for provenance in every summary.
TRACKED_PACKAGES: tuple[str, ...] = (
    "bbconf", "bedboss", "geniml", "huggingface_hub", "polars",
    "pybedtools", "pymemesuite", "pybigwig", "gtars", "umap-learn",
    "pyyaml", "pandas", "numpy", "scipy", "scikit-learn",
)


# File / config helpers -------------------------------------------------------

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_record(path: Path, record_count: int | None = None) -> dict[str, Any]:
    """Compact dict describing a file, for inputs/outputs in a summary."""
    p = Path(path)
    try:
        rel = p.resolve().relative_to(REPO_ROOT)
        path_str = str(rel)
    except ValueError:
        path_str = str(p)
    exists = p.exists()
    rec: dict[str, Any] = {
        "path": path_str,
        "exists": exists,
        "size_bytes": p.stat().st_size if exists else None,
        "sha256": sha256_file(p) if exists else None,
    }
    if record_count is not None:
        rec["record_count"] = record_count
    return rec


def load_config(config_path: Path) -> dict[str, Any]:
    with open(config_path) as f:
        return yaml.safe_load(f) or {}


def resolve_stage_cfg(cfg: dict[str, Any], stage_name: str) -> dict[str, Any]:
    """Return the stage's section merged with top-level paths/tissues/etc for convenience."""
    stages = cfg.get("stages", {}) or {}
    stage = stages.get(stage_name, {}) or {}
    top_level_shared = {
        k: v for k, v in cfg.items()
        if k not in ("stages", "project", "slurm")
    }
    return {**top_level_shared, **stage}


def config_hash(cfg: dict[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(cfg, sort_keys=True, default=str).encode()
    ).hexdigest()[:16]


# Provenance helpers ----------------------------------------------------------

def versions() -> dict[str, str]:
    out: dict[str, str] = {"python": sys.version.split()[0]}
    for pkg in TRACKED_PACKAGES:
        try:
            out[pkg] = importlib.metadata.version(pkg)
        except importlib.metadata.PackageNotFoundError:
            out[pkg] = "not-installed"
    return out


def git_commit() -> str | None:
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True, cwd=REPO_ROOT,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "diff", "--quiet"],
            cwd=REPO_ROOT,
        ).returncode != 0
        return f"{sha}-dirty" if dirty else sha
    except Exception:
        return None


# Stage context + lifecycle ---------------------------------------------------

@dataclass
class StageContext:
    name: str
    config_path: Path
    cfg: dict[str, Any]              # full config
    stage_cfg: dict[str, Any]        # resolved stage section + shared top-level keys
    results_dir: Path                # genomic-dict/results/<name>/
    start_time: float

    def path(self, key: str) -> Path:
        """Project-absolute path from a config 'paths' entry."""
        paths = self.stage_cfg.get("paths", {}) or {}
        rel = paths.get(key)
        if rel is None:
            raise KeyError(f"paths.{key} not in config {self.config_path}")
        return PROJECT_ROOT / rel


def stage_start(stage_name: str, docstring: str | None = None) -> StageContext:
    """Parse --config, load YAML, resolve stage section, create results dir, announce to stderr."""
    parser = argparse.ArgumentParser(
        description=(docstring or stage_name).strip().splitlines()[0] if docstring else stage_name,
    )
    parser.add_argument(
        "--config", type=Path, default=DEFAULT_CONFIG,
        help=f"Path to genomic-dict config.yaml (default: {DEFAULT_CONFIG})",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    stage_cfg = resolve_stage_cfg(cfg, stage_name)

    results_dir_rel = (cfg.get("paths", {}) or {}).get("results_dir", "results")
    results_dir = PROJECT_ROOT / results_dir_rel / stage_name
    results_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== {stage_name} ===", file=sys.stderr)
    if docstring:
        print(docstring.strip(), file=sys.stderr)
        print("", file=sys.stderr)
    print(f"config:  {args.config}", file=sys.stderr)
    print(f"results: {results_dir}", file=sys.stderr)
    print("", file=sys.stderr)

    return StageContext(
        name=stage_name,
        config_path=args.config,
        cfg=cfg,
        stage_cfg=stage_cfg,
        results_dir=results_dir,
        start_time=time.time(),
    )


def write_summary(
    ctx: StageContext,
    status: str,
    *,
    inputs: list[dict[str, Any]] | None = None,
    outputs: list[dict[str, Any]] | None = None,
    metrics: dict[str, Any] | None = None,
    warnings: list[str] | None = None,
    error: dict[str, Any] | None = None,
) -> Path:
    """Write results/<stage>/summary.json with full provenance."""
    summary: dict[str, Any] = {
        "stage": ctx.name,
        "status": status,
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "duration_seconds": round(time.time() - ctx.start_time, 2),
        "inputs": inputs or [],
        "outputs": outputs or [],
        "metrics": metrics or {},
        "warnings": warnings or [],
        "error": error,
        "config_file": str(ctx.config_path),
        "config_hash": config_hash(ctx.cfg),
        "config_used": ctx.stage_cfg,
        "versions": versions(),
        "git_commit": git_commit(),
    }
    path = ctx.results_dir / "summary.json"
    path.write_text(json.dumps(summary, indent=2, sort_keys=True, default=str))
    print(f"wrote {path.relative_to(REPO_ROOT)}", file=sys.stderr)
    return path
