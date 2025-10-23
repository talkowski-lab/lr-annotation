#!/usr/bin/env python3

import argparse
import gzip
from collections import defaultdict
from typing import Dict, Tuple, List

import pysam
import pandas as pd


def is_af_field(key: str) -> bool:
    lower_key = key.lower()
    return lower_key == "af" or lower_key.startswith("af_") or lower_key.endswith("_af")


def normalize_af_value(value):
    if isinstance(value, tuple) or isinstance(value, list):
        return value[0] if value else None
    return value


def normalize_af_field(field: str) -> frozenset:
    parts = field.lower().replace("af_", "").replace("_af", "").split("_")
    normalized_parts = set()
    for part in parts:
        if part == "male":
            normalized_parts.add("xy")
        elif part == "female":
            normalized_parts.add("xx")
        else:
            normalized_parts.add(part)
    return frozenset(normalized_parts)


def parse_truth_info_string(info_str: str) -> Dict[str, str]:
    items = info_str.strip().split(";") if info_str else []
    kv = {}
    for it in items:
        if "=" in it:
            k, v = it.split("=", 1)
            kv[k] = v
        else:
            kv[it] = ""
    return kv


def load_truth_info(tsv_gz_paths):
    truth = {}
    for p in tsv_gz_paths:
        if p is None:
            continue
        with gzip.open(p, "rt") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) != 2:
                    continue
                vid, info = parts
                truth[vid] = parse_truth_info_string(info)
    return truth


def parse_vep_header_line(header_line: str) -> Tuple[str, List[str]]:
    line = header_line.strip()
    lower = line.lower()
    if "id=csq" in lower:
        id_key = "CSQ"
    elif "id=vep" in lower:
        id_key = "VEP"
    else:
        raise ValueError("VEP/CSQ ID not found in header line")

    fmt_part = line.split("Format:")[-1].strip().strip('">')
    fmt_fields = [f.strip().lower() for f in fmt_part.split("|")]
    return id_key, fmt_fields


def get_eval_vep_header_from_vcf(vcf_path: str) -> Tuple[str, List[str]]:
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf.header.records:
            if rec.key == "INFO" and rec.get("ID"):
                vep_id = rec.get("ID")
                if vep_id and vep_id.lower() in ["vep", "csq"]:
                    desc = rec.get("Description", "")
                    fmt = desc.split("Format:")[-1].strip().strip('"').lower()
                    return vep_id, [f.strip() for f in fmt.split("|")]
    raise ValueError("Could not find eval VEP/CSQ header in VCF")


def get_vep_annotations(
    info_dict, vep_key: str, indices: Dict[int, str]
) -> Dict[str, str]:
    vep_string = ""
    if vep_key in info_dict:
        v = info_dict.get(vep_key, "")
        if isinstance(v, (tuple, list)):
            vep_string = v[0] if v else ""
        else:
            vep_string = v
    fields = str(vep_string).lower().split("|") if vep_string else []
    ann = {}
    for idx, cat in indices.items():
        ann[cat] = fields[idx] if idx < len(fields) and fields[idx] else "N/A"
    return ann


def write_summary_stats(final_vcf_path: str, contig: str, output_path: str):
    total_variants = 0
    match_counts = defaultdict(int)
    with pysam.VariantFile(final_vcf_path) as vcf_in:
        for record in vcf_in:
            total_variants += 1
            if "gnomAD_V4_match" in record.info:
                match_type = record.info["gnomAD_V4_match"]
                match_counts[match_type] += 1
            else:
                match_counts["UNMATCHED"] += 1
    total_matched = total_variants - match_counts["UNMATCHED"]
    summary_stats = {
        "contig": contig,
        "total_variants": total_variants,
        "total_matched": total_matched,
        "total_unmatched": match_counts["UNMATCHED"],
        "percent_matched": (
            (total_matched / total_variants * 100) if total_variants > 0 else 0.0
        ),
        "percent_unmatched": (
            (match_counts["UNMATCHED"] / total_variants * 100)
            if total_variants > 0
            else 0.0
        ),
    }
    for match_type, count in match_counts.items():
        if match_type != "UNMATCHED":
            summary_stats[f"{str(match_type).lower()}_count"] = count
            summary_stats[f"{str(match_type).lower()}_percent"] = (
                (count / total_variants * 100) if total_variants > 0 else 0.0
            )
    with open(output_path, "w") as f:
        f.write("\t".join(summary_stats.keys()) + "\n")
        f.write("\t".join(str(v) for v in summary_stats.values()) + "\n")


def write_summary_table(
    final_vcf_path: str,
    truth_variants: Dict[str, Dict[str, str]],
    vep_keys: Tuple[str, str, Dict[int, str], Dict[int, str]],
    output_path: str,
):
    all_rows = []
    eval_vep_key, truth_vep_key, eval_indices, truth_indices = vep_keys

    with pysam.VariantFile(final_vcf_path) as vcf_in:
        for record in vcf_in:
            row = {
                "eval_variant_id": record.id,
                "match_status": False,
                "truth_variant_id": ".",
            }
            if "gnomAD_V4_match_ID" in record.info:
                match_id = record.info["gnomAD_V4_match_ID"]
                if match_id in truth_variants:
                    row["match_status"] = True
                    row["truth_variant_id"] = match_id

                    truth_info = truth_variants[match_id]
                    eval_af_pairs = {
                        normalize_af_field(k): normalize_af_value(v)
                        for k, v in record.info.items()
                        if is_af_field(k)
                    }
                    truth_af_pairs_raw = {
                        k: v for k, v in truth_info.items() if is_af_field(k)
                    }
                    truth_af_pairs = {}
                    for k, v in truth_af_pairs_raw.items():
                        try:
                            if isinstance(v, str) and "," in v:
                                v = v.split(",")[0]
                            truth_af_pairs[normalize_af_field(k)] = float(v)
                        except Exception:
                            continue
                    common_af_keys = set(eval_af_pairs.keys()) & set(
                        truth_af_pairs.keys()
                    )
                    for af_key_set in common_af_keys:
                        af_key_str = "_".join(sorted(list(af_key_set)))
                        row[f"{af_key_str}_eval"] = eval_af_pairs[af_key_set]
                        row[f"{af_key_str}_truth"] = truth_af_pairs[af_key_set]

                    eval_annos = get_vep_annotations(
                        record.info, eval_vep_key, eval_indices
                    )
                    truth_annos = get_vep_annotations(
                        truth_info, truth_vep_key, truth_indices
                    )
                    for category in set(eval_indices.values()) | set(
                        truth_indices.values()
                    ):
                        row[f"{category}_eval"] = eval_annos.get(category, "N/A")
                        row[f"{category}_truth"] = truth_annos.get(category, "N/A")

            all_rows.append(row)

    df = pd.DataFrame(all_rows)
    cols_to_drop = []
    prefixes = set()
    for col in df.columns:
        if col.endswith("_eval") or col.endswith("_truth"):
            prefixes.add(col.rsplit("_", 1)[0])
    for prefix in prefixes:
        ec = f"{prefix}_eval"
        tc = f"{prefix}_truth"
        if ec in df.columns and tc in df.columns:
            eval_empty = (df[ec].isna() | (df[ec] == "N/A")).all()
            truth_empty = (df[tc].isna() | (df[tc] == "N/A")).all()
            if eval_empty and truth_empty:
                cols_to_drop.extend([ec, tc])
    df.drop(columns=cols_to_drop, inplace=True, errors="ignore")
    df.to_csv(output_path, sep="\t", index=False, na_rep=".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--final_vcf", required=True)
    ap.add_argument("--truth_tsv_snv", required=True)
    ap.add_argument("--truth_tsv_sv", required=True)
    ap.add_argument("--truth_vep_header", required=True)
    args = ap.parse_args()

    truth_info = load_truth_info([args.truth_tsv_snv, args.truth_tsv_sv])

    with open(args.truth_vep_header, "r") as f:
        truth_header_line = f.readline()
    truth_vep_key, truth_fields = parse_vep_header_line(truth_header_line)
    eval_vep_key, eval_fields = get_eval_vep_header_from_vcf(args.final_vcf)

    common_categories = set(eval_fields) & set(truth_fields)
    eval_indices = {
        i: cat for i, cat in enumerate(eval_fields) if cat in common_categories
    }
    truth_indices = {
        i: cat for i, cat in enumerate(truth_fields) if cat in common_categories
    }

    summary_table_path = f"{args.prefix}.benchmark_summary.tsv"
    write_summary_table(
        args.final_vcf,
        truth_info,
        (eval_vep_key, truth_vep_key, eval_indices, truth_indices),
        summary_table_path,
    )

    summary_stats_path = f"{args.prefix}.summary_stats.tsv"
    write_summary_stats(args.final_vcf, args.contig, summary_stats_path)


if __name__ == "__main__":
    main()
