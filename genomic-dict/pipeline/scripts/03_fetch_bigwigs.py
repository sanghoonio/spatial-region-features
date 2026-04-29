"""Stage 03 — fetch ENCODE bigwigs for 3 tissues × (5 histone marks + ATAC).

For each (biosample, assay, target) triple in stages.03_fetch_bigwigs.targets:
  1. Query the ENCODE search API for matching released GRCh38 experiments.
  2. Pick the top experiment; fetch its full file list.
  3. Filter to file_format=bigWig, output_type="fold change over control",
     assembly=GRCh38, status=released. Pick the most recent date_created.
  4. Download to data/annotations/bigwigs/<biosample>__<target-or-ATAC>.bigWig.
  5. Record ENCFF accession, URL, sha256, size in the output manifest.

The script is resumable: previously-downloaded bigwigs are kept; previously-
resolved accessions are reused from data/annotations/bigwig_manifest.yaml.

Reads:  (nothing upstream — queries ENCODE portal)
Writes:
  data/annotations/bigwigs/*.bigWig              — downloaded bigwigs (~1–3 GB each)
  data/annotations/bigwig_manifest.yaml          — resolved ENCFF accessions + sha256
  results/03_fetch_bigwigs/summary.json          — AI-ingestible summary
"""
from __future__ import annotations

import hashlib
import sys
import time
from pathlib import Path
from typing import Any

import requests
import yaml

sys.path.insert(0, str(Path(__file__).parent))
from _common import (  # noqa: E402
    PROJECT_ROOT, REPO_ROOT, file_record, stage_start, write_summary,
)


def target_key(biosample: str, assay: str, target: str | None) -> str:
    """Stable filesystem-friendly key for a (biosample, assay, target) triple."""
    tag = target if target else "ATAC"
    return f"{biosample}__{tag}"


def search_experiment(
    portal: str, biosample: str, assay: str, target: str | None, assembly: str,
) -> dict[str, Any] | None:
    """Query ENCODE search API; return the top matching experiment JSON or None."""
    params = {
        "type": "Experiment",
        "assay_title": assay,
        "biosample_ontology.term_name": biosample,
        "status": "released",
        "assembly": assembly,
        "format": "json",
        "limit": "5",
    }
    if target:
        params["target.label"] = target
    r = requests.get(
        f"{portal}/search/",
        params=params,
        headers={"Accept": "application/json"},
        timeout=30,
    )
    r.raise_for_status()
    graph = r.json().get("@graph", [])
    return graph[0] if graph else None


def pick_bigwig(
    portal: str, experiment_accession: str, preferred_output_type: str, assembly: str,
) -> dict[str, Any] | None:
    """Return the chosen bigwig file JSON for an experiment, or None."""
    r = requests.get(
        f"{portal}/experiments/{experiment_accession}/",
        headers={"Accept": "application/json"},
        timeout=30,
    )
    r.raise_for_status()
    files = r.json().get("files", [])
    matches = [
        f for f in files
        if f.get("file_format") == "bigWig"
        and f.get("output_type") == preferred_output_type
        and f.get("assembly") == assembly
        and f.get("status") == "released"
    ]
    if not matches:
        return None
    matches.sort(key=lambda f: f.get("date_created", ""), reverse=True)
    return matches[0]


def sha256_stream(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download_streaming(url: str, dest: Path) -> tuple[str, int]:
    """Stream-download url to dest; return (sha256, size). Idempotent.

    Deletes any stale .part file on entry. Connect timeout 10s; per-chunk
    read timeout 300s (needed for ~1 GB downloads with occasional TCP stalls).
    """
    if dest.exists():
        return sha256_stream(dest), dest.stat().st_size
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()
    h = hashlib.sha256()
    size = 0
    with requests.get(url, stream=True, timeout=(10, 300)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as out:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if not chunk:
                    continue
                out.write(chunk)
                h.update(chunk)
                size += len(chunk)
    tmp.rename(dest)
    return h.hexdigest(), size


def load_existing_manifest(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text()) or {}
    return data.get("entries", {}) or {}


def save_manifest(path: Path, entries: dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    doc = {"entries": entries}
    path.write_text(yaml.safe_dump(doc, sort_keys=True, width=120))


def main() -> None:
    ctx = stage_start("03_fetch_bigwigs", __doc__)
    sc = ctx.stage_cfg

    portal = sc["encode_portal"].rstrip("/")
    assembly = sc["assembly"]
    preferred = sc["preferred_output_type"]
    targets = sc["targets"]

    annotations_dir = PROJECT_ROOT / sc["paths"]["annotations_dir"]
    bigwig_dir = annotations_dir / "bigwigs"
    bigwig_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = annotations_dir / "bigwig_manifest.yaml"

    existing = load_existing_manifest(manifest_path)
    entries: dict[str, dict] = dict(existing)

    outputs: list[dict] = []
    warnings: list[str] = []

    try:
        # Phase 1 — resolve all target accessions from ENCODE (fast, no download).
        # Persist after each resolution so an interrupted run can resume.
        print("Phase 1: resolving accessions", file=sys.stderr)
        for t in targets:
            key = target_key(t["biosample"], t["assay"], t.get("target"))
            prev = entries.get(key)
            if prev and prev.get("file_accession"):
                continue  # already resolved

            print(f"  {key}: resolving...", file=sys.stderr)
            try:
                exp = search_experiment(portal, t["biosample"], t["assay"], t.get("target"), assembly)
            except Exception as e:
                warnings.append(f"{key}: search failed — {e}")
                continue
            if not exp:
                warnings.append(f"{key}: no experiment found")
                continue
            exp_acc = exp["accession"]

            try:
                bw = pick_bigwig(portal, exp_acc, preferred, assembly)
            except Exception as e:
                warnings.append(f"{key}: file listing failed for {exp_acc} — {e}")
                continue
            if not bw:
                warnings.append(f"{key}: no matching bigwig in {exp_acc}")
                continue

            file_acc = bw["accession"]
            href = bw.get("href", f"/files/{file_acc}/@@download/{file_acc}.bigWig")
            url = f"{portal}{href}"
            local_path = bigwig_dir / f"{key}.bigWig"

            entries[key] = {
                "biosample": t["biosample"],
                "assay": t["assay"],
                "target": t.get("target"),
                "experiment_accession": exp_acc,
                "file_accession": file_acc,
                "assembly": assembly,
                "output_type": preferred,
                "url": url,
                "local_path": str(local_path.relative_to(REPO_ROOT)),
                "expected_size_bytes": int(bw.get("file_size") or 0),
                "date_created": bw.get("date_created"),
                # sha256 and size_bytes filled in during download phase.
            }
            save_manifest(manifest_path, entries)
            time.sleep(0.25)  # be polite to the portal

        # Phase 2 — download each resolved bigwig. Skip if already on disk.
        print("\nPhase 2: downloading bigwigs", file=sys.stderr)
        for t in targets:
            key = target_key(t["biosample"], t["assay"], t.get("target"))
            entry = entries.get(key)
            if not entry or not entry.get("file_accession"):
                continue  # couldn't resolve in phase 1

            local_path = bigwig_dir / f"{key}.bigWig"
            expected_size = entry.get("expected_size_bytes") or 0

            if local_path.exists() and (expected_size == 0 or local_path.stat().st_size == expected_size):
                # Already downloaded (possibly from a prior run).
                if not entry.get("sha256"):
                    entry["sha256"] = sha256_stream(local_path)
                    entry["size_bytes"] = local_path.stat().st_size
                    save_manifest(manifest_path, entries)
                print(f"  {key}: cached {entry['file_accession']}", file=sys.stderr)
                outputs.append(file_record(local_path))
                continue

            print(
                f"  {key}: downloading {entry['file_accession']} ({expected_size/1e6:.0f} MB)...",
                file=sys.stderr,
            )
            try:
                sha, size = download_streaming(entry["url"], local_path)
            except Exception as e:
                warnings.append(f"{key}: download failed — {e}")
                continue

            entry["sha256"] = sha
            entry["size_bytes"] = int(size)
            save_manifest(manifest_path, entries)
            outputs.append(file_record(local_path))

        outputs.append(file_record(manifest_path))

        n_requested = len(targets)
        n_resolved = sum(
            1 for t in targets
            if target_key(t["biosample"], t["assay"], t.get("target")) in entries
            and entries[target_key(t["biosample"], t["assay"], t.get("target"))].get("file_accession")
        )
        total_bytes = sum(int(v.get("size_bytes") or 0) for v in entries.values())

        metrics = {
            "n_targets_requested": n_requested,
            "n_resolved": n_resolved,
            "n_missing": n_requested - n_resolved,
            "total_bytes_downloaded": total_bytes,
            "total_gb_downloaded": round(total_bytes / 1e9, 2),
            "per_target": {
                target_key(t["biosample"], t["assay"], t.get("target")): {
                    "biosample": t["biosample"],
                    "assay": t["assay"],
                    "target": t.get("target"),
                    "file_accession": (entries.get(target_key(t["biosample"], t["assay"], t.get("target"))) or {}).get("file_accession"),
                    "size_bytes": (entries.get(target_key(t["biosample"], t["assay"], t.get("target"))) or {}).get("size_bytes"),
                }
                for t in targets
            },
        }

        status = "success" if n_resolved == n_requested and not warnings else "partial"
        write_summary(
            ctx, status,
            outputs=outputs,
            metrics=metrics,
            warnings=warnings,
        )
        print(
            f"\n  resolved {n_resolved}/{n_requested} bigwigs; "
            f"{total_bytes/1e9:.1f} GB on disk",
            file=sys.stderr,
        )
        if warnings:
            print(f"  {len(warnings)} warnings — see summary.json", file=sys.stderr)

    except Exception as e:
        # Preserve whatever we've resolved so far.
        save_manifest(manifest_path, entries)
        write_summary(
            ctx, "error",
            outputs=outputs,
            warnings=warnings,
            error={
                "class": type(e).__name__,
                "message": str(e),
                "suggestion": (
                    "Check network to encodeproject.org. "
                    "Script is idempotent — rerun to resume from cached state."
                ),
                "retryable": True,
            },
        )
        raise


if __name__ == "__main__":
    main()
