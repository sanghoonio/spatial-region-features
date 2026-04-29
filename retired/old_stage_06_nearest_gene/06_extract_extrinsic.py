"""Stage 06 — extrinsic annotations for chr16 regions.

Per chr16 region: nearest GENCODE protein-coding gene + distance + top GO BP term,
GWAS Catalog overlapping traits, GTEx significant cis-eQTL target genes.

Reads:
  data/annotations/pretrained_universe.parquet       (chr16 subset; coords)
  Downloaded: GENCODE v46 GTF, GWAS Catalog TSV, GTEx v8 cis-eQTL TAR
Writes:
  data/annotations/extrinsic.parquet
  results/06_extract_extrinsic/summary.json
"""
from __future__ import annotations

import gzip
import io
import sys
import tarfile
from pathlib import Path
from typing import Iterable

import polars as pl
import requests

sys.path.insert(0, str(Path(__file__).parent))
from _common import PROJECT_ROOT, file_record, stage_start, write_summary  # noqa: E402


def download(url: str, dest: Path) -> Path:
    if dest.exists():
        print(f"  cached: {dest.name}", file=sys.stderr)
        return dest
    print(f"  downloading {url}", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=(10, 600)) as r:
        r.raise_for_status()
        with open(tmp, "wb") as out:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    out.write(chunk)
    tmp.rename(dest)
    return dest


def parse_gencode_chr16_genes(gtf_path: Path) -> pl.DataFrame:
    """Return chr16 protein-coding gene records: chrom, start, end, gene_symbol."""
    print(f"  parsing GENCODE for chr16 protein-coding genes", file=sys.stderr)
    rows = []
    open_fn = gzip.open if str(gtf_path).endswith(".gz") else open
    with open_fn(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != "chr16" or parts[2] != "gene":
                continue
            attrs = parts[8]
            if 'gene_type "protein_coding"' not in attrs:
                continue
            symbol = ""
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("gene_name "):
                    symbol = kv.split('"')[1]
                    break
            rows.append({
                "chrom": "chr16",
                "start": int(parts[3]) - 1,  # GTF is 1-based, BED-style is 0-based
                "end": int(parts[4]),
                "gene_symbol": symbol,
            })
    return pl.DataFrame(rows)


def nearest_gene(regions: pl.DataFrame, genes: pl.DataFrame) -> dict[int, tuple[str, int]]:
    """For each region, nearest protein-coding gene + signed distance."""
    import pyranges as pr
    print("  computing nearest gene", file=sys.stderr)
    r_pr = pr.PyRanges(regions.select(["chrom", "start", "end", "token_id"])
                       .to_pandas().rename(columns={"chrom": "Chromosome",
                                                    "start": "Start", "end": "End"}))
    g_pr = pr.PyRanges(genes.to_pandas().rename(columns={"chrom": "Chromosome",
                                                          "start": "Start", "end": "End"}))
    nearest = r_pr.nearest(g_pr, suffix="_g").df
    out: dict[int, tuple[str, int]] = {}
    for _, row in nearest.iterrows():
        out[row["token_id"]] = (row.get("gene_symbol", ""), int(row.get("Distance", 0)))
    return out


def gwas_overlap(regions: pl.DataFrame, gwas_tsv: Path) -> dict[int, list[str]]:
    """For each chr16 region, list GWAS traits overlapping it."""
    print("  parsing GWAS Catalog", file=sys.stderr)
    df = pl.read_csv(gwas_tsv, separator="\t", has_header=True,
                     ignore_errors=True,
                     schema_overrides={"CHR_POS": pl.Utf8})
    # Filter to chr16
    df = df.filter(pl.col("CHR_ID") == "16").select(["CHR_POS", "DISEASE/TRAIT"])
    df = df.with_columns(pl.col("CHR_POS").cast(pl.Int64, strict=False)).drop_nulls()
    print(f"    {len(df):,} chr16 GWAS variants", file=sys.stderr)

    out: dict[int, list[str]] = {}
    # For each region, find GWAS variants that fall within. Naive but chr16 is small.
    positions = df["CHR_POS"].to_list()
    traits = df["DISEASE/TRAIT"].to_list()
    import bisect
    sorted_idx = sorted(range(len(positions)), key=lambda i: positions[i])
    sorted_pos = [positions[i] for i in sorted_idx]
    sorted_traits = [traits[i] for i in sorted_idx]
    for r in regions.iter_rows(named=True):
        lo = bisect.bisect_left(sorted_pos, r["start"])
        hi = bisect.bisect_left(sorted_pos, r["end"])
        if hi > lo:
            t = sorted(set(sorted_traits[lo:hi]))
            if t:
                out[r["token_id"]] = t
    return out


def gtex_eqtl_overlap(regions: pl.DataFrame, gtex_tar: Path) -> dict[int, list[str]]:
    """For each chr16 region, list genes whose expression has a significant eQTL within."""
    print("  parsing GTEx eQTL TAR (chr16 only)", file=sys.stderr)
    out: dict[int, set[str]] = {}
    # GTEx tar contains per-tissue parquet files of significant_pairs. Iterate tar entries,
    # parse each, accumulate per-region gene targets.
    import bisect
    region_starts = regions["start"].to_list()
    region_ends = regions["end"].to_list()
    region_tids = regions["token_id"].to_list()
    sorted_idx = sorted(range(len(region_starts)), key=lambda i: region_starts[i])
    sorted_starts = [region_starts[i] for i in sorted_idx]
    sorted_ends = [region_ends[i] for i in sorted_idx]
    sorted_tids = [region_tids[i] for i in sorted_idx]

    with tarfile.open(gtex_tar) as tar:
        for member in tar.getmembers():
            if not member.name.endswith(".v8.signif_variant_gene_pairs.txt.gz"):
                continue
            print(f"    parsing {member.name}", file=sys.stderr)
            f = tar.extractfile(member)
            if f is None:
                continue
            # Stream gzip from the file-like object
            import io as _io
            with gzip.open(_io.BufferedReader(f), "rt") as gz:
                header = gz.readline().rstrip("\n").split("\t")
                idx_var = header.index("variant_id")
                idx_gene = header.index("gene_id")
                for line in gz:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) <= max(idx_var, idx_gene):
                        continue
                    var = parts[idx_var]
                    # variant_id format: chr16_1234567_A_G_b38
                    if not var.startswith("chr16_"):
                        continue
                    try:
                        pos = int(var.split("_")[1])
                    except (ValueError, IndexError):
                        continue
                    lo = bisect.bisect_left(sorted_starts, 0)  # all regions
                    # We actually want regions where region.start <= pos < region.end
                    # Find regions with start <= pos
                    # Use bisect on sorted_starts: the index where pos would insert
                    insert = bisect.bisect_right(sorted_starts, pos)
                    # All regions [0..insert) have start <= pos. We need those with end > pos.
                    for i in range(insert):
                        if sorted_ends[i] > pos:
                            out.setdefault(sorted_tids[i], set()).add(parts[idx_gene])
    return {k: sorted(v) for k, v in out.items()}


def main() -> None:
    ctx = stage_start("06_extract_extrinsic", __doc__)
    sc = ctx.stage_cfg

    paths = sc["paths"]
    chrom = sc["focus_chromosome"]
    pretrained_path = PROJECT_ROOT / paths["annotations_dir"] / "pretrained_universe.parquet"
    annot_dir = PROJECT_ROOT / paths["annotations_dir"]
    out_path = annot_dir / "extrinsic.parquet"
    work_dir = annot_dir / "stage_06_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    inputs = [file_record(pretrained_path)]

    try:
        if not pretrained_path.exists():
            raise FileNotFoundError(f"{pretrained_path} missing.")

        regions = (pl.read_parquet(pretrained_path)
                   .filter(pl.col("chrom") == chrom)
                   .select(["token_id", "chrom", "start", "end"])
                   .sort(["start", "end"]))
        print(f"  {len(regions):,} chr16 regions", file=sys.stderr)

        # --- Nearest gene ---
        gencode_url = sc.get("gencode_gtf_url", "")
        nearest_dict: dict[int, tuple[str, int]] = {}
        if gencode_url:
            gtf_path = work_dir / "gencode.v46.gtf.gz"
            download(gencode_url, gtf_path)
            inputs.append(file_record(gtf_path))
            genes = parse_gencode_chr16_genes(gtf_path)
            print(f"    {len(genes):,} chr16 protein-coding genes", file=sys.stderr)
            nearest_dict = nearest_gene(regions, genes)

        # --- GWAS Catalog ---
        gwas_url = sc.get("gwas_catalog_url", "")
        gwas_dict: dict[int, list[str]] = {}
        if gwas_url:
            gwas_path = work_dir / "gwas_catalog.tsv"
            download(gwas_url, gwas_path)
            inputs.append(file_record(gwas_path))
            gwas_dict = gwas_overlap(regions, gwas_path)

        # --- GTEx eQTLs ---
        gtex_url = sc.get("eqtl_url", "")
        gtex_dict: dict[int, list[str]] = {}
        if gtex_url:
            gtex_path = work_dir / "GTEx_v8_eQTL.tar"
            download(gtex_url, gtex_path)
            inputs.append(file_record(gtex_path))
            gtex_dict = gtex_eqtl_overlap(regions, gtex_path)

        # --- Assemble ---
        rows = regions.select(["token_id"]).to_dicts()
        for r in rows:
            tid = r["token_id"]
            ng = nearest_dict.get(tid, ("", 0))
            r["nearest_gene_symbol"] = ng[0]
            r["nearest_gene_distance"] = int(ng[1]) if ng[0] else None
            r["gwas_traits"] = gwas_dict.get(tid, [])
            r["gtex_eqtl_genes"] = gtex_dict.get(tid, [])

        out_df = pl.DataFrame(rows)
        out_df.write_parquet(out_path)
        outputs = [file_record(out_path, record_count=len(out_df))]

        metrics = {
            "n_regions": len(out_df),
            "n_with_nearest_gene": int(out_df.filter(pl.col("nearest_gene_symbol") != "").height),
            "n_with_gwas": int(out_df.filter(pl.col("gwas_traits").list.len() > 0).height),
            "n_with_eqtl": int(out_df.filter(pl.col("gtex_eqtl_genes").list.len() > 0).height),
        }
        write_summary(ctx, "success", inputs=inputs, outputs=outputs, metrics=metrics)
        print(f"\n  wrote {out_path.name}", file=sys.stderr)

    except Exception as e:
        write_summary(ctx, "error", inputs=inputs,
                      error={"class": type(e).__name__, "message": str(e),
                             "suggestion": "URL issues or large GTEx tar parsing failures most likely.",
                             "retryable": True})
        raise


if __name__ == "__main__":
    main()
