#!/usr/bin/env python3

import argparse
import gzip
from typing import Dict, Tuple, List, Set

import pysam
import pandas as pd
from collections import defaultdict, Counter


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
                    fmt_fields = [f.strip() for f in fmt.split("|")]
                    return vep_id, fmt_fields
    raise ValueError("Could not find eval VEP/CSQ header in VCF")


def get_case_insensitive(info_dict: Dict[str, object], key: str):
    target = key.lower()
    for k, v in info_dict.items():
        try:
            if str(k).lower() == target:
                return v
        except Exception:
            continue
    return None


def extract_vep_annotations(
    info_dict: Dict[str, object], vep_key: str, indices: Dict[int, str]
) -> Dict[str, str]:
    vep_string = ""
    v = get_case_insensitive(info_dict, vep_key)
    if v is not None:
        if isinstance(v, tuple) or isinstance(v, list):
            vep_string = v[0] if v else ""
        else:
            vep_string = v
    fields = str(vep_string).lower().split("|") if vep_string else []
    annos = {}
    for idx, cat in indices.items():
        if idx < len(fields) and fields[idx]:
            annos[cat] = fields[idx]
        else:
            annos[cat] = "N/A"
    return annos


def load_truth_info(tsv_gz_paths: List[str]) -> Dict[str, Dict[str, str]]:
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


def load_shard_matches(path: str) -> Tuple[Dict[str, str], Set[str]]:
    mapping = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            eval_id, truth_id = parts[0], parts[1]
            mapping[eval_id] = truth_id
    return mapping, set(mapping.keys())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--final_vcf", required=True)
    ap.add_argument("--matched_shard_tsv", required=True)
    ap.add_argument("--truth_tsv_snv", required=True)
    ap.add_argument("--truth_tsv_sv", required=True)
    ap.add_argument("--truth_vep_header", required=True)
    ap.add_argument("--shard_label", required=True)
    ap.add_argument(
        "--skip_vep_categories", help="Comma-separated list of VEP categories to skip"
    )
    args = ap.parse_args()

    af_out_path = (
        f"{args.prefix}.{args.contig}.shard_{args.shard_label}.af_pairs.tsv.gz"
    )
    vep_out_path = (
        f"{args.prefix}.{args.contig}.shard_{args.shard_label}.vep_pairs.tsv.gz"
    )

    id_map, shard_eval_ids = load_shard_matches(args.matched_shard_tsv)
    truth_info = load_truth_info([args.truth_tsv_snv, args.truth_tsv_sv])

    with open(args.truth_vep_header, "r") as f:
        truth_header_line = f.readline()
    truth_vep_key, truth_vep_fields = parse_vep_header_line(truth_header_line)

    eval_vep_key, eval_vep_fields = get_eval_vep_header_from_vcf(args.final_vcf)

    # Parse skip_vep_categories
    skip_categories = set()
    if args.skip_vep_categories:
        skip_categories = {
            cat.strip().lower() for cat in args.skip_vep_categories.split(",")
        }
        print(f"Skipping VEP categories: {', '.join(sorted(skip_categories))}")

    # Filter out skipped categories
    common_categories = set(eval_vep_fields) & set(truth_vep_fields)
    common_categories = {cat for cat in common_categories if cat not in skip_categories}

    eval_indices = {
        i: cat for i, cat in enumerate(eval_vep_fields) if cat in common_categories
    }
    truth_indices = {
        i: cat for i, cat in enumerate(truth_vep_fields) if cat in common_categories
    }

    af_rows = []
    vep_counts: Dict[str, Counter] = defaultdict(Counter)

    with pysam.VariantFile(args.final_vcf) as vcf_in:
        for rec in vcf_in:
            if rec.id not in shard_eval_ids:
                continue
            if (
                "gnomAD_V4_match" not in rec.info
                or "gnomAD_V4_match_ID" not in rec.info
            ):
                continue
            truth_id = id_map.get(rec.id)
            if not truth_id:
                continue
            truth_rec_info = truth_info.get(truth_id)
            if not truth_rec_info:
                continue

            eval_af_pairs = {
                normalize_af_field(k): normalize_af_value(v)
                for k, v in rec.info.items()
                if is_af_field(k)
            }
            truth_af_pairs_raw = {
                k: v for k, v in truth_rec_info.items() if is_af_field(k)
            }
            truth_af_pairs = {}
            for k, v in truth_af_pairs_raw.items():
                if "," in v:
                    v = v.split(",")[0]
                truth_af_pairs[normalize_af_field(k)] = float(v)

            for key_set, eval_val in eval_af_pairs.items():
                if key_set in truth_af_pairs and eval_val is not None:
                    e = float(eval_val)
                    t = float(truth_af_pairs[key_set])
                    if e > 0 and t > 0:
                        af_rows.append(
                            {
                                "af_key": "_".join(sorted(list(key_set))),
                                "eval_af": e,
                                "truth_af": t,
                            }
                        )

            eval_ann = extract_vep_annotations(rec.info, eval_vep_key, eval_indices)
            truth_ann = extract_vep_annotations(
                truth_rec_info, truth_vep_key, truth_indices
            )
            for cat in common_categories:
                eval_val = eval_ann.get(cat, "N/A")
                truth_val = truth_ann.get(cat, "N/A")
                # if eval_val == 'N/A' and truth_val == 'N/A':
                #     continue
                vep_counts[cat][(eval_val, truth_val)] += 1

    if af_rows:
        df_af = pd.DataFrame(af_rows)
        with gzip.open(af_out_path, "wt") as f:
            df_af.to_csv(f, sep="\t", index=False)
    else:
        with gzip.open(af_out_path, "wt") as f:
            f.write("af_key\teval_af\ttruth_af\n")

    rows = []
    for cat, ctr in vep_counts.items():
        for (e, t), c in ctr.items():
            rows.append({"category": cat, "eval": e, "truth": t, "count": c})
    df_vep = pd.DataFrame(rows, columns=["category", "eval", "truth", "count"])
    with gzip.open(vep_out_path, "wt") as f:
        df_vep.to_csv(f, sep="\t", index=False)


if __name__ == "__main__":
    main()
