import argparse
import os
import tarfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import re
import subprocess
from pysam import VariantFile
import tempfile
import numpy as np

# --- Helper Functions ---
def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def get_vep_annotations(info_dict, vep_key, vep_format):
    vep_string = info_dict.get(vep_key, '')
    if not vep_string: return []
    annotations = []
    vep_fields = vep_format.split('|')
    for annotation_blob in vep_string.split(','):
        values = annotation_blob.split('|')
        annotations.append(dict(zip(vep_fields, values)))
    return annotations

def normalize_af_field(field):
    parts = field.lower().split('_')
    normalized_parts = set()
    for part in parts:
        if part == 'male': normalized_parts.add('xy')
        elif part == 'female': normalized_parts.add('xx')
        elif part != 'af': normalized_parts.add(part)
    return frozenset(normalized_parts)

def is_af_field(key):
    lower_key = key.lower()
    return lower_key.startswith('af_') or lower_key.endswith('_af')

def run_truvari(vcf_eval, vcf_truth, pct_seq, output_dir, ref_fasta):
    cmd = ["truvari", "bench", "-b", vcf_truth, "-c", vcf_eval, "-f", ref_fasta, "-o", output_dir, "--pctseq", str(pct_seq), "--passonly"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return os.path.join(output_dir, "tp-comp.vcf.gz")

def parse_closest_bed(bed_path):
    """Parse the output of the R script to find matches."""
    matches = {}
    with open(bed_path) as f:
        for line in f:
            if 'query_svid' in line: continue # Skip header
            fields = line.strip().split('\t')
            if len(fields) > 1:
                query_id = fields[0]
                truth_id = fields[1]
                # The R script output indicates no match with "NA"
                if truth_id != "NA":
                    matches[query_id] = ("BEDTOOLS_CLOSEST", truth_id)
    return matches

# --- Main Annotation and Benchmarking Logic ---
def benchmark_annotations(vcf_eval_path, vcf_truth_path, sv_truth_path, closest_bed_path, contig, prefix, create_benchmarks, ref_fasta_path):
    
    # --- Phase 1: Exact Matching ---
    print("Starting Phase 1: Exact Matching")
    truth_variants = {}
    with VariantFile(vcf_truth_path) as vcf_in:
        for rec in vcf_in:
            for alt in rec.alts:
                truth_variants[(rec.chrom, rec.pos, rec.ref, alt)] = rec

    matches = {}
    unmatched_eval_records = []
    eval_records_map = {}
    with VariantFile(vcf_eval_path) as vcf_in:
        for rec in vcf_in:
            eval_records_map[rec.id] = rec
            is_matched = False
            for alt in rec.alts:
                key = (rec.chrom, rec.pos, rec.ref, alt)
                if key in truth_variants:
                    matches[rec.id] = ("EXACT", truth_variants[key])
                    is_matched = True
                    del truth_variants[key]
                    break
            if not is_matched:
                unmatched_eval_records.append(rec)
    print(f"Found {len(matches)} exact matches.")

    # --- Phase 2: Truvari Matching ---
    print("Starting Phase 2: Truvari Matching")
    sv_unmatched_records = [rec for rec in unmatched_eval_records if 'SVLEN' in rec.info and abs(rec.info['SVLEN'][0]) >= 10]
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".vcf") as f:
        with VariantFile(vcf_eval_path) as vcf_in:
            with VariantFile(f.name, 'w', header=vcf_in.header) as vcf_out:
                for rec in sv_unmatched_records: vcf_out.write(rec)
        sv_unmatched_vcf = f.name
    
    non_snv_truth_records = {rec.id: rec for rec in truth_variants.values() if rec.info.get('SVTYPE') != 'SNV'}
    
    if sv_unmatched_records and non_snv_truth_records:
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".vcf") as f:
            with VariantFile(vcf_truth_path) as vcf_in:
                with VariantFile(f.name, 'w', header=vcf_in.header) as vcf_out:
                    for rec in non_snv_truth_records.values(): vcf_out.write(rec)
            remaining_truth_vcf = f.name
            
        for pct_seq in [0.9, 0.7, 0.5]:
            print(f"Running Truvari with pctseq={pct_seq}")
            truvari_dir = f"truvari_{pct_seq}"
            tp_comp_vcf = run_truvari(sv_unmatched_vcf, remaining_truth_vcf, pct_seq, truvari_dir, ref_fasta_path)
            
            newly_matched_eval_ids = set()
            if os.path.exists(tp_comp_vcf):
                with VariantFile(tp_comp_vcf) as vcf_in:
                    for rec in vcf_in:
                        base_id = rec.info['MatchId'].split('.')[1]
                        matches[rec.id] = (f"TRUVARI_{pct_seq}", non_snv_truth_records[base_id])
                        newly_matched_eval_ids.add(rec.id)
                        if base_id in non_snv_truth_records:
                            del non_snv_truth_records[base_id]

            sv_unmatched_records = [rec for rec in sv_unmatched_records if rec.id not in newly_matched_eval_ids]
            if not sv_unmatched_records or not non_snv_truth_records: break
            
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".vcf") as f:
                with VariantFile(vcf_truth_path) as vcf_in:
                    with VariantFile(f.name, 'w', header=vcf_in.header) as vcf_out:
                        for rec in non_snv_truth_records.values(): vcf_out.write(rec)
                remaining_truth_vcf = f.name

    # --- Phase 3: Bedtools Closest ---
    print("Starting Phase 3: Bedtools Closest")
    bedtools_matches = parse_closest_bed(closest_bed_path)
    sv_truth_map = {rec.id: rec for rec in VariantFile(sv_truth_path)}
    for rec in sv_unmatched_records:
        if rec.id in bedtools_matches:
            truth_id = bedtools_matches[rec.id][1]
            if truth_id in sv_truth_map:
                matches[rec.id] = ("BEDTOOLS_CLOSEST", sv_truth_map[truth_id])
    print(f"Found {len(bedtools_matches)} potential Bedtools Closest matches.")

    # --- Final Annotation Step ---
    print("Starting final annotation.")
    annotated_vcf_path = f"{prefix}.annotated.vcf"
    with VariantFile(vcf_eval_path) as vcf_in:
        header = vcf_in.header
        header.info.add('gnomAD_V4_match', '.', 'String', 'Matching status in gnomAD v4')
        with VariantFile(annotated_vcf_path, 'w', header=header) as vcf_out:
            for rec in vcf_in:
                if rec.id in matches:
                    rec.info['gnomAD_V4_match'] = matches[rec.id][0]
                vcf_out.write(rec)
    subprocess.run(f"bgzip -f {annotated_vcf_path} && tabix -p vcf {annotated_vcf_path}.gz", shell=True, check=True)

    # --- Optional Benchmarking ---
    if create_benchmarks and matches:
        # ... (Full benchmarking logic with AF plots and VEP heatmaps will be implemented here)
        print("Benchmarking logic not fully implemented in this step.")

def main():
    parser = argparse.ArgumentParser(description="Benchmark and annotate SNV/SV annotations.")
    # ... (arguments as before)
    # This is a stub for the final script structure
    pass

if __name__ == "__main__":
    main() 