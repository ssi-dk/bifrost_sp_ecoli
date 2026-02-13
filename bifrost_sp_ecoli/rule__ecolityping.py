#!/usr/bin/env python3
import argparse
import os
import sys
from io import StringIO
from typing import Dict, List, Tuple, Optional
import pandas as pd
import yaml

CORE_PREFIXES = ["stx", "wzx", "wzy", "wzt", "wzm", "flic", "fli", "fl", "eae", "ehxa"]

FINAL_DETAIL_SPECS: List[Tuple[str, str]] = [
    ("toxin_details", "toxin"),
    ("O_type_details", "Otype"),
    ("H_type_details", "Htype"),
    ("type_details", "adhesin_virulence"),
    ("other_details", "other"),
]

# ------------------------- Parsing & Threshold resolution ------------------------- #

def parse_gene_from_template(template: str) -> Tuple[str, str]:
    NA = "-"
    parts = template.split("__") if "__" in template else template.split("_")
    gene = parts[1] if len(parts) >= 2 else parts[0]
    allele = parts[2] if len(parts) >= 3 else NA
    return gene, allele


def resolve_threshold_key_for_gene(gene: str, thresholds: Dict[str, List[float]]) -> str:
    gene_str = str(gene)
    gene_l = gene_str.lower()

    if gene_str in thresholds:
        return gene_str

    for k in thresholds:
        if k.lower() == gene_l:
            return k

    candidates = [
        k for k in thresholds
        if k.lower() != "other" and gene_l.startswith(k.lower())
    ]
    if candidates:
        return max(candidates, key=lambda k: (len(k), k.lower()))

    for k in thresholds:
        if k.lower() == "other":
            return k

    raise ValueError(f"No threshold for gene '{gene_str}' and no 'other' fallback in config.")

# ------------------------- Step 1: Read YAML gene config ------------------------- #

def read_geneconfig(config_path: str, organism_key: str) -> Dict[str, List[float]]:
    if yaml is None:
        raise RuntimeError("Missing dependency 'pyyaml'. Install with: pip install pyyaml")

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, "r") as fh:
        cfg = yaml.safe_load(fh)

    if not isinstance(cfg, dict):
        raise ValueError("Config YAML must be a mapping at the top level (dict).")

    if organism_key not in cfg:
        raise ValueError(
            f"Organism key '{organism_key}' not found in config. Available: {list(cfg.keys())}"
        )

    section = cfg[organism_key]
    if not isinstance(section, dict):
        raise ValueError(f"Config section for '{organism_key}' must be a mapping of gene -> [cov,id,depth].")

    thresholds: Dict[str, List[float]] = {}
    for gene, vals in section.items():
        if not isinstance(gene, str):
            raise ValueError(f"Gene name must be a string. Got: {gene!r}")
        if not (isinstance(vals, list) and len(vals) == 3):
            raise ValueError(f"Threshold for gene '{gene}' must be a 3-item list: [cov,id,depth]. Got: {vals!r}")
        thresholds[gene] = [float(vals[0]), float(vals[1]), float(vals[2])]

    return thresholds

# ------------------------- Step 2: Read + Extract KMA .res ------------------------- #

def process_kma_res(res_path: str) -> pd.DataFrame:
    if not os.path.exists(res_path):
        raise FileNotFoundError(f"KMA .res file not found: {res_path}")

    with open(res_path, "r") as fh:
        lines = fh.readlines()

    header_idx: Optional[int] = None
    for i, line in enumerate(lines):
        if line.strip().startswith("#Template"):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("Could not find header line starting with '#Template' in the .res file.")

    content = "".join(lines[header_idx:])
    df = pd.read_csv(StringIO(content), sep=r"\s+", engine="python")

    required_cols = {"#Template", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in .res: {sorted(missing)}. Found: {df.columns.tolist()}")

    parsed = df["#Template"].apply(parse_gene_from_template)
    df["gene"] = parsed.apply(lambda x: x[0])
    df["allele"] = parsed.apply(lambda x: x[1])

    return df[["gene", "allele", "Template_Coverage", "Query_Identity", "Depth"]].copy()


def filter_kma_res(df: pd.DataFrame, thresholds: Dict[str, List[float]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"filter_kma_res: input df missing columns: {sorted(missing)}")

    work = df.copy()

    matched_keys: List[str] = [
        resolve_threshold_key_for_gene(str(g), thresholds)
        for g in work["gene"].tolist()
    ]
    work["matched_key"] = matched_keys

    th_df = pd.DataFrame.from_dict(
        thresholds, orient="index", columns=["th_cov", "th_id", "th_depth"]
    )
    work = work.join(th_df, on="matched_key")

    cov = work["Template_Coverage"].astype(float)
    ide = work["Query_Identity"].astype(float)
    dep = work["Depth"].astype(float)

    pass_mask = (cov >= work["th_cov"]) & (ide >= work["th_id"]) & (dep >= work["th_depth"])
    work["status"] = ["PASS" if x else "FAIL" for x in pass_mask.tolist()]

    out_cols = [
        "gene", "allele", "matched_key",
        "Template_Coverage", "Query_Identity", "Depth",
        "th_cov", "th_id", "th_depth", "status"
    ]
    pass_df = work.loc[pass_mask, out_cols].copy()
    fail_df = work.loc[~pass_mask, out_cols].copy()

    return pass_df, fail_df

# ------------------------- Helpers for details ------------------------- #

def collect_alleles(df: pd.DataFrame, prefix: str) -> List[str]:
    genes_lower = df["gene"].astype(str).str.lower()
    alleles = df["allele"].astype(str)

    mask = genes_lower.str.startswith(prefix.lower())
    out: List[str] = []
    for a in alleles[mask].tolist():
        if a and a != "-":
            out.append(a)
    return out


def build_details_from_df(df: pd.DataFrame, prefixes: List[str], include_gene: bool = False) -> str:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"build_details_from_df: df missing columns: {sorted(missing)}")

    if df.empty:
        return "-"

    genes = df["gene"].astype(str)
    genes_lower = genes.str.lower()

    mask = pd.Series([False] * len(df), index=df.index)
    for p in prefixes:
        mask = mask | genes_lower.str.startswith(p.lower())

    if not mask.any():
        return "-"

    subset = df.loc[mask, ["gene", "allele", "Template_Coverage", "Query_Identity", "Depth"]]
    parts: List[str] = []
    for row in subset.itertuples(index=False):
        gene = str(row.gene)
        allele = str(row.allele)
        cov = float(row.Template_Coverage)
        qid = float(row.Query_Identity)
        dep = float(row.Depth)

        if include_gene:
            parts.append(f"{gene}_{allele}_{cov:.2f}_{qid:.2f}_{dep:.2f}")
        else:
            parts.append(f"{allele}_{cov:.2f}_{qid:.2f}_{dep:.2f}")

    return "-" if len(parts) == 0 else ";".join(parts)

# ------------------------- Typing functions ------------------------- #

def determine_stx_subtype(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)

    if missing_pass:
        raise ValueError(f"determine_stx_subtype: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_stx_subtype: fail_df missing columns: {sorted(missing_fail)}")

    stx1_alleles = collect_alleles(pass_df, "stx1")
    stx2_alleles = collect_alleles(pass_df, "stx2")

    stx1_col = "-" if len(stx1_alleles) == 0 else ";".join(stx1_alleles)
    stx2_col = "-" if len(stx2_alleles) == 0 else ";".join(stx2_alleles)

    stx1_unique = sorted(set(stx1_alleles))
    stx2_unique = sorted(set(stx2_alleles))

    toxin_parts: List[str] = []
    if len(stx1_unique) == 1:
        toxin_parts.append(stx1_unique[0])
    if len(stx2_unique) == 1:
        toxin_parts.append(stx2_unique[0])

    toxin_call = "-" if len(toxin_parts) == 0 else ";".join(toxin_parts)

    toxin_pass_details = build_details_from_df(pass_df, ["stx"], include_gene=False)
    toxin_fail_details = build_details_from_df(fail_df, ["stx"], include_gene=False)

    toxin_details = f"pass:{toxin_pass_details}|fail:{toxin_fail_details}"

    stx_out = pd.DataFrame([{
        "sample_id": sample_id,
        "stx1": stx1_col,
        "stx2": stx2_col,
        "Toxin": toxin_call,
        "toxin_pass_details": toxin_pass_details,
        "toxin_fail_details": toxin_fail_details,
        "toxin_details": toxin_details,
    }])

    return stx_out


def determine_O_type(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_O_type: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_O_type: fail_df missing columns: {sorted(missing_fail)}")

    wzx_alleles = collect_alleles(pass_df, "wzx")
    wzy_alleles = collect_alleles(pass_df, "wzy")
    wzt_alleles = collect_alleles(pass_df, "wzt")
    wzm_alleles = collect_alleles(pass_df, "wzm")

    wzx_col = "-" if len(wzx_alleles) == 0 else ";".join(wzx_alleles)
    wzy_col = "-" if len(wzy_alleles) == 0 else ";".join(wzy_alleles)
    wzt_col = "-" if len(wzt_alleles) == 0 else ";".join(wzt_alleles)
    wzm_col = "-" if len(wzm_alleles) == 0 else ";".join(wzm_alleles)

    wzx_unique = sorted(set(wzx_alleles))
    wzy_unique = sorted(set(wzy_alleles))
    wzt_unique = sorted(set(wzt_alleles))
    wzm_unique = sorted(set(wzm_alleles))

    O_candidates: List[str] = []

    if len(wzx_unique) == 1 and len(wzy_unique) == 1 and wzx_unique[0] == wzy_unique[0]:
        O_candidates.append(wzx_unique[0])

    if len(wzt_unique) == 1 and len(wzm_unique) == 1 and wzt_unique[0] == wzm_unique[0]:
        O_candidates.append(wzt_unique[0])

    O_type = "-"
    if O_candidates and len(set(O_candidates)) == 1:
        O_type = O_candidates[0]

    prefixes = ["wzx", "wzy", "wzt", "wzm"]
    O_type_pass_details = build_details_from_df(pass_df, prefixes, include_gene=True)
    O_type_fail_details = build_details_from_df(fail_df, prefixes, include_gene=True)
    O_type_details = f"pass:{O_type_pass_details}|fail:{O_type_fail_details}"

    out = pd.DataFrame([{
        "sample_id": sample_id,
        "wzx": wzx_col,
        "wzy": wzy_col,
        "wzt": wzt_col,
        "wzm": wzm_col,
        "O_type": O_type,
        "O_type_pass_details": O_type_pass_details,
        "O_type_fail_details": O_type_fail_details,
        "O_type_details": O_type_details,
    }])

    return out


def determine_H_type(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_H_type: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_H_type: fail_df missing columns: {sorted(missing_fail)}")

    fliC_alleles = collect_alleles(pass_df, "fliC")

    if len(fliC_alleles) > 0:
        use_prefixes = ["fliC"]
        fl_alleles = fliC_alleles
    else:
        use_prefixes = ["fl"]
        fl_alleles = collect_alleles(pass_df, "fl")

    flagellin = "-" if len(fl_alleles) == 0 else ";".join(fl_alleles)

    unique = sorted(set(fl_alleles))
    if len(unique) == 1:
        H_type = unique[0]
    else:
        H_type = "-"

    H_type_pass_details = build_details_from_df(pass_df, use_prefixes, include_gene=True)
    H_type_fail_details = build_details_from_df(fail_df, use_prefixes, include_gene=True)
    H_type_details = f"pass:{H_type_pass_details}|fail:{H_type_fail_details}"

    out = pd.DataFrame([{
        "sample_id": sample_id,
        "flagellin": flagellin,
        "H_type": H_type,
        "H_type_pass_details": H_type_pass_details,
        "H_type_fail_details": H_type_fail_details,
        "H_type_details": H_type_details,
    }])

    return out


def determine_adhesin_virulence(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_adhesin_virulence: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_adhesin_virulence: fail_df missing columns: {sorted(missing_fail)}")

    eae_positive = collect_alleles(pass_df, "eae")
    ehxA_positive = collect_alleles(pass_df, "ehxA")

    Adhesin = "positive" if eae_positive else "-"
    Virulence = "positive" if ehxA_positive else "-"

    prefixes = ["eae", "ehxA"]
    type_pass_details = build_details_from_df(pass_df, prefixes, include_gene=True)
    type_fail_details = build_details_from_df(fail_df, prefixes, include_gene=True)
    type_details = f"pass:{type_pass_details}|fail:{type_fail_details}"

    out = pd.DataFrame([{
        "sample_id": sample_id,
        "Adhesin": Adhesin,
        "Virulence": Virulence,
        "type_pass_details": type_pass_details,
        "type_fail_details": type_fail_details,
        "type_details": type_details,
    }])

    return out

# ------------------------- Final merge ------------------------- #

def is_missing_detail_value(v: object) -> bool:
    if v is None:
        return True
    if isinstance(v, float) and pd.isna(v):
        return True
    s = str(v).strip()
    return s == "" or s == "-"


def build_verbose_from_detail_columns(row: pd.Series) -> str:
    parts: List[str] = []
    for col, label in FINAL_DETAIL_SPECS:
        if col in row and not is_missing_detail_value(row[col]):
            parts.append(f"{label}:{row[col]}")
    return "-" if len(parts) == 0 else "||".join(parts)


def build_details_excluding_prefixes(df: pd.DataFrame, excluded_prefixes: List[str]) -> str:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"build_details_excluding_prefixes: df missing columns: {sorted(missing)}")

    if df.empty:
        return "-"

    genes_lower = df["gene"].astype(str).str.lower()

    mask_excl = pd.Series([False] * len(df), index=df.index)
    for p in excluded_prefixes:
        mask_excl = mask_excl | genes_lower.str.startswith(p.lower())

    subset = df.loc[~mask_excl, ["gene", "allele", "Template_Coverage", "Query_Identity", "Depth"]]
    if subset.empty:
        return "-"

    parts: List[str] = []
    for row in subset.itertuples(index=False):
        gene = str(row.gene)
        allele = str(row.allele)
        cov = float(row.Template_Coverage)
        qid = float(row.Query_Identity)
        dep = float(row.Depth)
        parts.append(f"{gene}_{allele}_{cov:.2f}_{qid:.2f}_{dep:.2f}")

    return "-" if len(parts) == 0 else ";".join(parts)


def determine_other(pass_df: pd.DataFrame, fail_df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    required = {"gene", "allele", "Template_Coverage", "Query_Identity", "Depth"}
    missing_pass = required - set(pass_df.columns)
    missing_fail = required - set(fail_df.columns)
    if missing_pass:
        raise ValueError(f"determine_other: pass_df missing columns: {sorted(missing_pass)}")
    if missing_fail:
        raise ValueError(f"determine_other: fail_df missing columns: {sorted(missing_fail)}")

    genes_lower = pass_df["gene"].astype(str).str.lower()
    mask_excl = pd.Series([False] * len(pass_df), index=pass_df.index)
    for p in CORE_PREFIXES:
        mask_excl = mask_excl | genes_lower.str.startswith(p)

    other_pass_genes = sorted(set(pass_df.loc[~mask_excl, "gene"].astype(str).tolist()))
    Other = "-" if len(other_pass_genes) == 0 else ";".join(other_pass_genes)

    other_pass_details = build_details_excluding_prefixes(pass_df, CORE_PREFIXES)
    other_fail_details = build_details_excluding_prefixes(fail_df, CORE_PREFIXES)
    other_details = f"pass:{other_pass_details}|fail:{other_fail_details}"

    return pd.DataFrame([{
        "sample_id": sample_id,
        "Other": Other,
        "other_pass_details": other_pass_details,
        "other_fail_details": other_fail_details,
        "other_details": other_details,
    }])


def merge_final_outputs(
    stx_out: pd.DataFrame,
    o_out: pd.DataFrame,
    h_out: pd.DataFrame,
    av_out: pd.DataFrame,
    other_out: pd.DataFrame,
) -> pd.DataFrame:
    final_df = stx_out.merge(o_out, on="sample_id", how="outer")
    final_df = final_df.merge(h_out, on="sample_id", how="outer")
    final_df = final_df.merge(av_out, on="sample_id", how="outer")
    final_df = final_df.merge(other_out, on="sample_id", how="outer")

    final_df["toxin_details_combined"] = final_df.apply(build_verbose_from_detail_columns, axis=1)

    o = final_df["O_type"].apply(lambda v: "-" if is_missing_detail_value(v) else str(v).strip())
    h = final_df["H_type"].apply(lambda v: "-" if is_missing_detail_value(v) else str(v).strip())
    final_df["sero_serotype_finder"] = o + ";" + h

    drop_cols: List[str] = []
    for c in final_df.columns:
        if c.endswith("_pass_details") or c.endswith("_fail_details"):
            drop_cols.append(c)

    for col, _label in FINAL_DETAIL_SPECS:
        if col in final_df.columns:
            drop_cols.append(col)

    for c in ["O_type", "H_type", "other", "Other", "verbose"]:
        if c in final_df.columns:
            drop_cols.append(c)

    drop_cols = sorted(set(drop_cols))
    if drop_cols:
        final_df = final_df.drop(columns=drop_cols)

    final_df = final_df.rename(columns={"toxin_details_combined": "toxin_details"})

    keep = ["sample_id", "Toxin", "sero_serotype_finder", "Adhesin", "Virulence", "toxin_details"]
    for c in keep:
        if c not in final_df.columns:
            final_df[c] = "-"

    final_df = final_df[keep].copy()
    final_df = final_df.fillna("-")

    return final_df

# ------------------------- Organism normalization ------------------------- #

def normalize_organism_key(user_organism: str) -> str:
    org = user_organism.strip()

    ecoli_aliases = {"Escherichia", "Escherichia coli", "E. coli", "E.coli"}
    shigella_aliases = {
        "Shigella",
        "Shigella sonnei", "S. sonnei", "S.sonnei",
        "Shigella flexneri", "S. flexneri", "S.flexneri",
        "Shigella boydii", "S. boydii", "S.boydii",
        "Shigella dysenteriae", "S. dysenteriae", "S.dysenteriae",
    }

    if org in ecoli_aliases:
        return "Escherichia"
    if org in shigella_aliases:
        return "Shigella"

    raise ValueError(f"Unrecognized organism '{user_organism}'.")

# ------------------------- CLI ------------------------- #

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Filter KMA .res for E. coli/Shigella typing and produce PASS/FAIL + summary TSVs. '--output' is a prefix."
    )
    p.add_argument("--KMA_res", required=True, help="Input KMA .res file")
    p.add_argument("--config", required=True, help="YAML config with per-organism gene thresholds")
    p.add_argument("--organism", required=True, help="Organism string (e.g., 'Escherichia coli', 'Shigella sonnei')")
    p.add_argument("--sample_id", required=True, help="Sample identifier")
    p.add_argument("--store_tmp", action="store_true", help="Write intermediate *_stx.tsv, *_O.tsv, *_H.tsv, *_adhesin_virulence.tsv, *_other.tsv")
    p.add_argument("--verbose", action="store_true", help="Write final merged *_final.tsv")
    p.add_argument("--output", required=True, help="Output prefix (directory/prefix; we append _pass.tsv, _fail.tsv, etc.)")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    organism_key = normalize_organism_key(args.organism)
    thresholds = read_geneconfig(args.config, organism_key)

    extracted_df = process_kma_res(args.KMA_res)
    pass_df, fail_df = filter_kma_res(extracted_df, thresholds)

    prefix = args.output
    out_dir = os.path.dirname(prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    pass_path = f"{prefix}_pass.tsv"
    fail_path = f"{prefix}_fail.tsv"

    pass_df.to_csv(pass_path, sep="\t", index=False)
    fail_df.to_csv(fail_path, sep="\t", index=False)

    print(f"Organism input: {args.organism!r} -> config key: {organism_key!r}", file=sys.stderr)
    print(f"Loaded {len(thresholds)} thresholds from {args.config}", file=sys.stderr)
    print(f"Read {len(extracted_df)} hits from {args.KMA_res}", file=sys.stderr)
    print(f"PASS: {len(pass_df)} -> {pass_path}", file=sys.stderr)
    print(f"FAIL: {len(fail_df)} -> {fail_path}", file=sys.stderr)

    stx_out = determine_stx_subtype(pass_df, fail_df, args.sample_id)
    o_out = determine_O_type(pass_df, fail_df, args.sample_id)
    h_out = determine_H_type(pass_df, fail_df, args.sample_id)
    av_out = determine_adhesin_virulence(pass_df, fail_df, args.sample_id)
    other_out = determine_other(pass_df, fail_df, args.sample_id)

    if args.store_tmp:
        stx_path = f"{prefix}_stx.tsv"
        stx_out.to_csv(stx_path, sep="\t", index=False)
        print(f"STX: wrote {stx_path}", file=sys.stderr)

        o_path = f"{prefix}_O.tsv"
        o_out.to_csv(o_path, sep="\t", index=False)
        print(f"O: wrote {o_path}", file=sys.stderr)

        h_path = f"{prefix}_H.tsv"
        h_out.to_csv(h_path, sep="\t", index=False)
        print(f"H: wrote {h_path}", file=sys.stderr)

        av_path = f"{prefix}_adhesin_virulence.tsv"
        av_out.to_csv(av_path, sep="\t", index=False)
        print(f"Adhesin/Virulence: wrote {av_path}", file=sys.stderr)

        other_path = f"{prefix}_other.tsv"
        other_out.to_csv(other_path, sep="\t", index=False)
        print(f"Other: wrote {other_path}", file=sys.stderr)
    else:
        print("store_tmp not set: skipping *_stx.tsv, *_O.tsv, *_H.tsv, *_adhesin_virulence.tsv, *_other.tsv", file=sys.stderr)

    if args.verbose:
        final_df = merge_final_outputs(stx_out, o_out, h_out, av_out, other_out)
        final_path = f"{prefix}_final.tsv"
        final_df.to_csv(final_path, sep="\t", index=False)
        print(f"FINAL: wrote {final_path}", file=sys.stderr)
    else:
        print("verbose not set: skipping *_final.tsv", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

