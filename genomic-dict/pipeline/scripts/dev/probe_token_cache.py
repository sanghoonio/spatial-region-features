"""Probe whether the pretrained R2V universe's tokenizations are already
cached in bbcache. If yes, stage 11 can read from cache (minutes instead of hours).

Run on Rivanna with the lab env sourced:
    source /project/shefflab/rivanna_config/env.sh
    uv run python genomic-dict/pipeline/scripts/probe_token_cache.py

What it does:
  1. Loads BBClient and the pretrained Region2VecExModel.
  2. Inspects BBClient + its underlying caches for any "tokens" surface.
  3. Tries common token-lookup API names against a sample BED ID from our
     stage-01 manifest, with the pretrained universe's identifier.
  4. Reports findings: API name (if any), sample token count, cache size.
"""
from __future__ import annotations

import os
import sys
import inspect
from pathlib import Path

# We use the pretrained model name as configured.
HF_MODEL = "databio/r2v-encode-hg38"

# Sample file IDs to test against — pulled from the stage 01 manifest
# (these are real corpus IDs; any of them should be in bbcache).
SAMPLE_BED_IDS = [
    "00694b547e9941b2e24e4fe2d6de240d",
    "5cecd5b72c6837a1b3e64eb2d9bdc4a4",
    "126d5e7e0651e36f2a4c5f3c8f5c1e1d",
]


def banner(title: str) -> None:
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def main() -> int:
    if not os.environ.get("BBCLIENT_CACHE"):
        print("ERROR: BBCLIENT_CACHE not set. Source the lab env first.", file=sys.stderr)
        return 1

    banner("1. Loading BBClient + pretrained model")
    from geniml.bbclient import BBClient
    from geniml.region2vec.main import Region2VecExModel

    bbc = BBClient()
    print(f"   bbcache root: {os.environ['BBCLIENT_CACHE']}")
    print(f"   loading {HF_MODEL}...")
    model = Region2VecExModel(model_path=HF_MODEL)
    print(f"   loaded; tokenizer type = {type(model.tokenizer).__name__}")

    # --- 2. What attributes / methods does BBClient expose around tokens?
    banner("2. BBClient surface (token-related members)")
    for name in sorted(dir(bbc)):
        if name.startswith("_"):
            continue
        if "tok" in name.lower() or "cache" in name.lower():
            attr = getattr(bbc, name)
            kind = "method" if callable(attr) else "attr"
            try:
                sig = str(inspect.signature(attr)) if callable(attr) else ""
            except (ValueError, TypeError):
                sig = ""
            print(f"   {name:30s} ({kind}) {sig}")

    # --- 3. Tokenizer / universe identifier
    banner("3. Probing pretrained-tokenizer / universe identifier")
    tok = model.tokenizer
    candidates = ["universe_id", "universe", "name", "config", "vocab_size"]
    for c in candidates:
        if hasattr(tok, c):
            try:
                val = getattr(tok, c)
                if callable(val):
                    val = val()
                preview = str(val)[:200]
                print(f"   tokenizer.{c} = {preview}")
            except Exception as e:
                print(f"   tokenizer.{c} -> error: {e}")

    # --- 4. Inspect the bedfile cache subfolder for any "tokens"-like dirs
    banner("4. bbcache directory layout")
    cache_root = Path(os.environ["BBCLIENT_CACHE"])
    for entry in sorted(cache_root.iterdir())[:30]:
        is_dir = entry.is_dir()
        size = ""
        if not is_dir:
            try:
                size = f"  {entry.stat().st_size:>12,} bytes"
            except OSError:
                pass
        print(f"   {entry.name:40s}{'  (dir)' if is_dir else ''}{size}")

    # --- 5. Try plausible token-lookup API names against a sample
    banner("5. Trying token-lookup APIs against sample bed IDs")
    universe_attr_for_lookup = None
    for c in ("universe_id", "universe", "name"):
        if hasattr(tok, c):
            v = getattr(tok, c)
            if not callable(v) and v:
                universe_attr_for_lookup = v
                print(f"   universe value used for lookup: {v!r}")
                break
    if universe_attr_for_lookup is None:
        print("   no obvious universe identifier on tokenizer; trying with HF_MODEL string")
        universe_attr_for_lookup = HF_MODEL

    candidate_methods = [
        "seek_tokens",
        "get_tokens",
        "load_tokens",
        "fetch_tokens",
        "tokens",
    ]
    for method_name in candidate_methods:
        if not hasattr(bbc, method_name):
            continue
        method = getattr(bbc, method_name)
        if not callable(method):
            continue
        print(f"\n   trying bbc.{method_name}(universe, bed_id) ...")
        for bed_id in SAMPLE_BED_IDS:
            try:
                result = method(universe_attr_for_lookup, bed_id)
                length = len(result) if hasattr(result, "__len__") else "?"
                preview = str(result)[:80] if length else "(empty)"
                print(f"     {bed_id[:12]}...  → length={length}  {preview}")
            except Exception as e:
                print(f"     {bed_id[:12]}...  → ERROR ({type(e).__name__}): {e}")

    # --- 6. Last resort: time how long it takes to tokenize one BED from scratch.
    banner("6. From-scratch tokenization timing (control)")
    import time

    sample_id = SAMPLE_BED_IDS[0]
    try:
        t0 = time.time()
        rs = bbc.load_bed(sample_id)
        t1 = time.time()
        toks = list(model.tokenizer.encode(model.tokenizer.tokenize(rs)))
        t2 = time.time()
        print(f"   load_bed({sample_id[:12]}...): {t1 - t0:.2f}s")
        print(f"   tokenize+encode: {t2 - t1:.2f}s")
        print(f"   total tokens emitted: {len(toks):,}")
        print(f"   per-file budget (84k files × {t2 - t0:.2f}s): "
              f"{(t2 - t0) * 84_698 / 60:.0f} min")
    except Exception as e:
        print(f"   from-scratch test failed: {e}")

    banner("Recap")
    print(
        "   If any candidate method in section 5 returned a token list,\n"
        "   stage 11 can read from cache. If everything errored, the lab\n"
        "   pipeline doesn't expose a public token-read API and we should\n"
        "   inspect the cache directory by hand (section 4) to look for\n"
        "   cached parquet/json files keyed by universe + bed id.\n"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
