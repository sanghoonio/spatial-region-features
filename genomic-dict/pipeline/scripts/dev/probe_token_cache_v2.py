"""V2 probe — find the universe_id under tokens.zarr and verify load_bed_tokens works.

Run on Rivanna with the lab env sourced:
    source /project/shefflab/rivanna_config/env.sh
    uv run python genomic-dict/pipeline/scripts/probe_token_cache_v2.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

HF_MODEL = "databio/r2v-encode-hg38"
SAMPLE_BED_IDS = [
    "00694b547e9941b2e24e4fe2d6de240d",
    "5cecd5b72c6837a1b3e64eb2d9bdc4a4",
]


def banner(t: str) -> None:
    print()
    print("=" * 70)
    print(f"  {t}")
    print("=" * 70)


def main() -> int:
    if not os.environ.get("BBCLIENT_CACHE"):
        print("ERROR: BBCLIENT_CACHE not set.", file=sys.stderr)
        return 1

    cache_root = Path(os.environ["BBCLIENT_CACHE"])
    zarr_path = cache_root / "tokens.zarr"

    banner("1. Inspect tokens.zarr keys (these are universe_ids)")
    import zarr

    z = zarr.open(str(zarr_path), mode="r")
    print(f"   zarr root: {zarr_path}")
    print(f"   type: {type(z).__name__}")

    # Top-level keys = candidate universe IDs
    universes = list(z.group_keys()) if hasattr(z, "group_keys") else []
    arrays = list(z.array_keys()) if hasattr(z, "array_keys") else []
    print(f"   {len(universes)} group_keys, {len(arrays)} array_keys")
    for k in universes[:30]:
        sub = z[k]
        n_kids = (
            len(list(sub.array_keys())) + len(list(sub.group_keys()))
            if hasattr(sub, "array_keys")
            else 0
        )
        print(f"     group:  {k!r:50s}  ~{n_kids} children")
    for k in arrays[:30]:
        a = z[k]
        print(f"     array:  {k!r:50s}  shape={a.shape} dtype={a.dtype}")

    if not universes and not arrays:
        print("   (zarr root is empty?)")
        return 0

    banner("2. For each universe_id, count cached BED files + sample-check")
    from geniml.bbclient import BBClient

    bbc = BBClient()

    for uni in universes[:10]:  # limit if many
        sub = z[uni]
        # Children are bed_ids
        bed_keys = list(sub.array_keys()) + list(sub.group_keys())
        print(f"\n   universe: {uni!r}")
        print(f"     cached bed count: {len(bed_keys):,}")
        if bed_keys[:3]:
            print(f"     sample bed keys: {bed_keys[:3]}")

        # Try load_bed_tokens with this universe + each sample
        for bed_id in SAMPLE_BED_IDS:
            try:
                arr = bbc.load_bed_tokens(bed_id, uni)
                print(f"     ✓ load_bed_tokens('{bed_id[:12]}...', {uni!r}) "
                      f"→ shape={arr.shape} dtype={arr.dtype}")
                # Show first few values
                vals = list(arr[: min(8, arr.shape[0])])
                print(f"       first values: {vals}")
            except Exception as e:
                print(f"     ✗ load_bed_tokens('{bed_id[:12]}...', {uni!r}) "
                      f"→ {type(e).__name__}: {e}")

    banner("3. Cross-check: are our 84k corpus IDs in any cached universe?")
    # Read manifest, take 50 random IDs, check coverage under each candidate universe
    import polars as pl

    manifest_path = (
        Path(__file__).resolve().parents[2]
        / "data"
        / "corpus"
        / "manifest.parquet"
    )
    if not manifest_path.exists():
        print(f"   (manifest not found at {manifest_path}; skipping)")
        return 0

    sample_ids = (
        pl.read_parquet(manifest_path)
        .select("id")
        .sample(n=50, seed=42)
        ["id"]
        .to_list()
    )

    for uni in universes[:10]:
        sub = z[uni]
        cached_set = set(sub.array_keys()) | set(sub.group_keys())
        present = sum(1 for x in sample_ids if x in cached_set)
        print(f"   universe {uni!r}: {present}/50 sample manifest IDs cached")

    return 0


if __name__ == "__main__":
    sys.exit(main())
