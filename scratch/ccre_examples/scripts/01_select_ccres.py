"""
Select 30 cCREs for the playground.

Stratification:
  PLS:  5 ubiquitous (active in K562 + GM12878 + HepG2) + 5 K562-specific
  pELS: 5 ubiquitous + 5 K562-specific
  dELS: 5 Fulco-2019-positive + 5 Fulco-2019-negative  (falls back to
        ubiquitous + K562-specific if data/fulco_2019_crispri.tsv is absent)

Inputs (auto-fetched into data/_cache/ on first run):
  - ENCODE hg38 V3 cCRE registry (ENCFF420VPZ.bed.gz)
  - K562 / GM12878 / HepG2 ATAC-seq fold-change-over-control bigwigs
    (ENCSR IDs at top of file; fold-change-over-control accession is resolved
    at runtime via the ENCODE REST API)

User-provided (optional):
  - data/fulco_2019_crispri.tsv  (Fulco et al. 2019 Supplementary Table 6a,
    columns: chrom, start, end, significant; see ../README.md)

Outputs:
  - data/ccres.bed             30 intervals, BED6 (name = cCRE accession)
  - data/ccres_metadata.csv    full metadata per selected cCRE
"""

from __future__ import annotations

import gzip
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig


# --- paths ---

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"
CACHE = DATA / "_cache"
FULCO_TSV = DATA / "fulco_2019_crispri.tsv"


# --- sources ---

REGISTRY_ACCESSION = "ENCFF420VPZ"  # hg38 V3 biosample-agnostic cCRE registry

# Canonical ENCODE ATAC-seq experiments per tissue. Verify on the ENCODE portal
# if a download errors; edit these IDs if you want a different experiment.
ATAC_EXPERIMENTS: dict[str, str] = {
    "K562":    "ENCSR483RKN",
    "GM12878": "ENCSR095QNB",
    "HepG2":   "ENCSR042AWH",
}

TISSUES = list(ATAC_EXPERIMENTS.keys())

# Engreitz lab ENCODE-E2G benchmark: Fulco 2019 CRISPRi results pre-lifted to
# GRCh38 and re-harmonized with other CRISPRi screens. We filter to the Fulco
# 2019 subset and use the `Regulated` column (significant + negative effect on
# target gene = canonical enhancer hit).
ENGREITZ_BENCHMARK_URL = (
    "https://raw.githubusercontent.com/EngreitzLab/CRISPR_comparison/main/"
    "resources/crispr_data/EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz"
)


# --- knobs ---

SEED = 42
N_PER_STRATUM = 5
ACTIVITY_PERCENTILE = 80.0  # cCRE is "active in tissue T" iff mean ATAC FC in
                            # that cCRE's interval is above the 80th percentile
                            # across all cCREs, computed per tissue.


# --- fetching ---

def curl(url: str, dest: Path) -> Path:
    """Idempotent curl download."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[curl] {url}")
    subprocess.run(["curl", "-fsSL", "-o", str(dest), url], check=True)
    size_mb = dest.stat().st_size / 1e6
    print(f"[done] {dest.name} ({size_mb:.1f} MB)")
    return dest


def api_get(url: str) -> dict:
    """GET a JSON document from the ENCODE portal."""
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def fetch_ccre_registry() -> pd.DataFrame:
    """Download the hg38 biosample-agnostic cCRE registry and return as a DataFrame."""
    url = f"https://www.encodeproject.org/files/{REGISTRY_ACCESSION}/@@download/{REGISTRY_ACCESSION}.bed.gz"
    dest = CACHE / f"{REGISTRY_ACCESSION}.bed.gz"
    curl(url, dest)

    # Registry BED columns (10):
    # chrom, start, end, ccre_accession, score, strand,
    # thickStart, thickEnd, itemRgb, ccre_class
    with gzip.open(dest, "rt") as f:
        df = pd.read_csv(
            f, sep="\t", header=None,
            usecols=[0, 1, 2, 3, 9],
            names=["chrom", "start", "end", "ccre_accession", "ccre_class"],
            dtype={"chrom": str, "ccre_accession": str, "ccre_class": str},
        )
    print(f"[read] {len(df):,} cCREs from {dest.name}")
    print("       class counts:")
    print(df["ccre_class"].value_counts().to_string())
    return df


def resolve_atac_bigwig(encsr_id: str) -> tuple[str, str]:
    """
    Look up the fold-change-over-control bigwig for an ATAC-seq experiment.
    Falls back to signal-p-value if no fold-change bigwig is released.

    Returns (download_url, output_type_used).
    """
    url = f"https://www.encodeproject.org/experiments/{encsr_id}/?format=json"
    print(f"[api ] {url}")
    exp = api_get(url)
    files = exp.get("files", [])

    def pick(output_type: str):
        cands = [
            f for f in files
            if f.get("output_type") == output_type
            and f.get("file_format") == "bigWig"
            and f.get("status") == "released"
            and f.get("assembly") == "GRCh38"
        ]
        if not cands:
            return None
        # prefer the ENCODE-marked default if present
        defaults = [f for f in cands if f.get("preferred_default")]
        return (defaults or cands)[0]

    chosen = pick("fold change over control") or pick("signal p-value")
    if chosen is None:
        raise RuntimeError(
            f"{encsr_id}: no released fold-change or signal-p-value bigwig on GRCh38"
        )
    acc = chosen["accession"]
    ot = chosen["output_type"]
    dl = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bigWig"
    print(f"[api ] {encsr_id}: {acc}  ({ot})")
    return dl, ot


def fetch_atac_bigwigs() -> dict[str, Path]:
    """Download the fold-change-over-control bigwig for each tissue's ATAC experiment."""
    paths: dict[str, Path] = {}
    for tissue, encsr in ATAC_EXPERIMENTS.items():
        url, _ = resolve_atac_bigwig(encsr)
        acc = url.rsplit("/", 1)[-1].removesuffix(".bigWig")
        dest = CACHE / f"{tissue}_ATAC_{acc}.bigWig"
        curl(url, dest)
        paths[tissue] = dest
    return paths


# --- scoring ---

def score_ccres_per_tissue(
    ccres: pd.DataFrame, bw_paths: dict[str, Path]
) -> pd.DataFrame:
    """Compute mean ATAC fold-change signal per cCRE per tissue. O(n_cCREs × n_tissues)."""
    cache = CACHE / "ccre_atac_scores.csv.gz"
    if cache.exists():
        print(f"[score] loading cached scores from {cache.name} "
              f"(delete to force re-score)")
        return pd.read_csv(cache)

    scores = ccres[["ccre_accession"]].copy()
    for tissue, path in bw_paths.items():
        print(f"[score] {tissue}: {path.name}")
        bw = pyBigWig.open(str(path))
        chroms_in_bw = set(bw.chroms().keys())
        n = len(ccres)
        vals = np.zeros(n, dtype=np.float32)
        missing_chroms: set[str] = set()

        chroms = ccres["chrom"].to_numpy()
        starts = ccres["start"].to_numpy()
        ends = ccres["end"].to_numpy()

        for i in range(n):
            c = chroms[i]
            if c not in chroms_in_bw:
                missing_chroms.add(c)
                continue
            try:
                m = bw.stats(c, int(starts[i]), int(ends[i]), type="mean")
                v = m[0] if m else None
                vals[i] = float(v) if v is not None else 0.0
            except RuntimeError:
                vals[i] = 0.0
            if i % 200_000 == 0 and i > 0:
                print(f"        {i:,}/{n:,} scored...")
        bw.close()
        if missing_chroms:
            print(f"        skipped chroms not in bigwig: {sorted(missing_chroms)[:5]}...")
        scores[f"atac_mean_{tissue}"] = vals

    scores.to_csv(cache, index=False, compression="gzip")
    print(f"[score] cached to {cache.name}")
    return scores


def classify_activity(scores: pd.DataFrame) -> pd.DataFrame:
    """Add per-tissue active/inactive booleans + two derived strata columns."""
    for tissue in TISSUES:
        col = f"atac_mean_{tissue}"
        thr = float(np.percentile(scores[col], ACTIVITY_PERCENTILE))
        scores[f"active_{tissue}"] = scores[col] > thr
        n_active = int(scores[f"active_{tissue}"].sum())
        print(f"[thr  ] {tissue}: threshold={thr:.3f}  n_active={n_active:,}")
    active_cols = [f"active_{t}" for t in TISSUES]
    scores["active_all3"] = scores[active_cols].all(axis=1)
    scores["active_K562_only"] = (
        scores["active_K562"]
        & ~scores["active_GM12878"]
        & ~scores["active_HepG2"]
    )
    print(f"[strata] ubiquitous (all 3)  n={int(scores['active_all3'].sum()):,}")
    print(f"[strata] K562-specific       n={int(scores['active_K562_only'].sum()):,}")
    return scores


# --- Fulco 2019 overlap ---

def _fetch_fulco_from_engreitz() -> None:
    """Download Engreitz's GRCh38-lifted benchmark and filter to Fulco 2019 rows."""
    src = CACHE / "EPCrisprBenchmark_training_K562_GRCh38.tsv.gz"
    curl(ENGREITZ_BENCHMARK_URL, src)
    bench = pd.read_csv(src, sep="\t", compression="gzip")
    is_fulco = bench["Reference"].str.contains("Fulco et al., 2019", na=False, regex=False)
    fulco = bench.loc[is_fulco].copy()
    out = pd.DataFrame({
        "chrom": fulco["chrom"],
        "start": fulco["chromStart"].astype(int),
        "end": fulco["chromEnd"].astype(int),
        "gene": fulco["measuredGeneSymbol"],
        "significant": fulco["Regulated"].astype(bool),
    })
    FULCO_TSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(FULCO_TSV, sep="\t", index=False)
    n_pos = int(out["significant"].sum())
    print(f"[write] {FULCO_TSV.name}  ({len(out):,} pairs, {n_pos:,} Regulated)")


def load_fulco() -> pd.DataFrame | None:
    if not FULCO_TSV.exists():
        print(f"[fetch] {FULCO_TSV.name} not local — pulling from Engreitz benchmark")
        try:
            _fetch_fulco_from_engreitz()
        except Exception as e:
            print(f"[warn ] Fulco auto-fetch failed: {e}")
            print("        Skipping validation stratification.")
            return None

    df = pd.read_csv(FULCO_TSV, sep="\t")
    required = {"chrom", "start", "end", "significant"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(
            f"{FULCO_TSV.name} missing columns: {sorted(missing)}. "
            f"Have: {sorted(df.columns)}."
        )
    df["significant"] = df["significant"].astype(bool)
    print(f"[read ] {len(df):,} Fulco enhancer-gene pairs "
          f"({int(df['significant'].sum()):,} Regulated)")
    return df[["chrom", "start", "end", "significant"]]


def tag_ccres_by_fulco(
    ccres: pd.DataFrame, fulco: pd.DataFrame
) -> pd.DataFrame:
    """For each cCRE, does it overlap any Fulco test, and any significant one?"""
    import pyranges as pr

    ccre_pr = pr.PyRanges(
        ccres[["chrom", "start", "end", "ccre_accession"]].rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
    )
    fulco_pr = pr.PyRanges(
        fulco.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
    )
    joined = ccre_pr.join(fulco_pr, how=None).df  # inner join — rows per overlap
    if len(joined) == 0:
        print("[fulco] no cCREs overlap any Fulco test "
              "(coordinate-assembly mismatch? Fulco is hg19 by default)")
        out = ccres[["ccre_accession"]].copy()
        out["fulco_tested"] = False
        out["fulco_positive"] = False
        return out

    agg = (
        joined.groupby("ccre_accession")
        .agg(fulco_tested=("significant", "size"),
             fulco_positive=("significant", "any"))
        .reset_index()
    )
    agg["fulco_tested"] = agg["fulco_tested"] > 0

    out = ccres[["ccre_accession"]].merge(agg, on="ccre_accession", how="left")
    out["fulco_tested"] = out["fulco_tested"].fillna(False).astype(bool)
    out["fulco_positive"] = out["fulco_positive"].fillna(False).astype(bool)

    print(f"[fulco] {int(out['fulco_tested'].sum()):,} cCREs overlap a Fulco test; "
          f"{int(out['fulco_positive'].sum()):,} overlap a significant test")
    return out


# --- selection ---

def _sample(df: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    if len(df) == 0:
        return df
    take = min(n, len(df))
    if take < n:
        print(f"[warn ] requested {n}, only {take} available — taking what we have")
    return df.sample(n=take, random_state=seed).reset_index(drop=True)


def select_ccres(
    registry: pd.DataFrame,
    scores: pd.DataFrame,
    fulco_tags: pd.DataFrame | None,
) -> pd.DataFrame:
    df = registry.merge(scores, on="ccre_accession", how="inner")
    if fulco_tags is not None:
        df = df.merge(fulco_tags, on="ccre_accession", how="left")
        df["fulco_tested"] = df["fulco_tested"].fillna(False).astype(bool)
        df["fulco_positive"] = df["fulco_positive"].fillna(False).astype(bool)
    else:
        df["fulco_tested"] = False
        df["fulco_positive"] = False

    picks: list[pd.DataFrame] = []

    for klass in ("PLS", "pELS"):
        sub = df[df["ccre_class"] == klass]
        ubiq = _sample(sub[sub["active_all3"]], N_PER_STRATUM, SEED)
        ubiq = ubiq.assign(stratum="ubiquitous", validation_status="untested")
        k562 = _sample(sub[sub["active_K562_only"]], N_PER_STRATUM, SEED + 1)
        k562 = k562.assign(stratum="K562_specific", validation_status="untested")
        picks += [ubiq, k562]

    dels = df[df["ccre_class"] == "dELS"]
    if fulco_tags is not None and int(dels["fulco_tested"].sum()) > 0:
        pos = _sample(dels[dels["fulco_positive"]], N_PER_STRATUM, SEED)
        pos = pos.assign(stratum="K562_tested", validation_status="Fulco_positive")
        neg_pool = dels[dels["fulco_tested"] & ~dels["fulco_positive"]]
        neg = _sample(neg_pool, N_PER_STRATUM, SEED + 1)
        neg = neg.assign(stratum="K562_tested", validation_status="Fulco_negative")
        picks += [pos, neg]
    else:
        print("[info ] no Fulco overlap available for dELS — falling back to "
              "tissue-specificity stratification for dELS too.")
        ubiq = _sample(dels[dels["active_all3"]], N_PER_STRATUM, SEED)
        ubiq = ubiq.assign(stratum="ubiquitous", validation_status="untested")
        k562 = _sample(dels[dels["active_K562_only"]], N_PER_STRATUM, SEED + 1)
        k562 = k562.assign(stratum="K562_specific", validation_status="untested")
        picks += [ubiq, k562]

    return pd.concat(picks, ignore_index=True)


# --- output ---

def write_outputs(selected: pd.DataFrame) -> None:
    DATA.mkdir(parents=True, exist_ok=True)

    bed = selected[["chrom", "start", "end", "ccre_accession"]].copy()
    bed["score"] = 0
    bed["strand"] = "."
    bed_path = DATA / "ccres.bed"
    bed.to_csv(bed_path, sep="\t", header=False, index=False)
    print(f"[write] {bed_path}  ({len(bed)} rows)")

    meta_cols = [
        "ccre_accession", "chrom", "start", "end",
        "ccre_class", "stratum", "validation_status",
        "atac_mean_K562", "atac_mean_GM12878", "atac_mean_HepG2",
        "active_K562", "active_GM12878", "active_HepG2",
        "active_all3", "active_K562_only",
        "fulco_tested", "fulco_positive",
    ]
    meta_path = DATA / "ccres_metadata.csv"
    selected[meta_cols].to_csv(meta_path, index=False)
    print(f"[write] {meta_path}  ({len(selected)} rows)")


# --- main ---

def main() -> int:
    print(f"=== Fetching cCRE registry ===")
    registry = fetch_ccre_registry()

    print(f"\n=== Fetching ATAC bigwigs ({', '.join(TISSUES)}) ===")
    bw_paths = fetch_atac_bigwigs()

    print(f"\n=== Scoring cCREs ===")
    scores = score_ccres_per_tissue(registry, bw_paths)

    print(f"\n=== Classifying activity (top {100 - ACTIVITY_PERCENTILE:.0f}% per tissue) ===")
    scores = classify_activity(scores)

    print(f"\n=== Loading Fulco 2019 CRISPRi table (optional) ===")
    fulco = load_fulco()
    fulco_tags = tag_ccres_by_fulco(registry, fulco) if fulco is not None else None

    print(f"\n=== Selecting cCREs ===")
    selected = select_ccres(registry, scores, fulco_tags)

    print(f"\n=== Selection summary ===")
    summary = (
        selected.groupby(["ccre_class", "stratum", "validation_status"])
        .size()
        .reset_index(name="n")
    )
    print(summary.to_string(index=False))
    print(f"\ntotal: {len(selected)} cCREs")

    write_outputs(selected)
    return 0


if __name__ == "__main__":
    sys.exit(main())
