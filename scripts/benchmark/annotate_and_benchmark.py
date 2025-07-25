#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import os
import tarfile
import subprocess
import numpy as np
import bgzip
from collections import defaultdict

def parse_info(info_str):
    """Parse VCF INFO string into a dictionary."""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def get_vep_format(vcf_path, vep_key_pattern="VEP"):
    """Extract VEP format from VCF header."""
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf.header.records:
            if record.key == 'INFO' and record.get('ID') and record.get('ID').upper() in ['VEP', 'CSQ']:
                desc = record.get('Description')
                if 'Format: ' in desc:
                    return desc.split('Format: ')[-1].strip('"')
    return None

def get_vep_annotations(info_dict, vep_key, vep_format):
    """Extract VEP annotations from INFO dictionary."""
    if not vep_format or vep_key not in info_dict:
        return []
    
    vep_string = info_dict.get(vep_key, '')
    annotations = []
    vep_fields = vep_format.split('|')
    for annotation_blob in vep_string.split(','):
        values = annotation_blob.split('|')
        annotations.append(dict(zip(vep_fields, values)))
    return annotations

def normalize_af_field(field):
    """Normalize AF field to a common format for comparison."""
    parts = field.lower().replace('af_', '').replace('_af', '').split('_')
    normalized_parts = set()
    for part in parts:
        if part == 'male':
            normalized_parts.add('xy')
        elif part == 'female':
            normalized_parts.add('xx')
        else:
            normalized_parts.add(part)
    return frozenset(normalized_parts)

def is_af_field(key):
    """Check if a key is an allele frequency field."""
    lower_key = key.lower()
    return lower_key.startswith('af_') or lower_key.endswith('_af')

def plot_af_correlation(df, pop, output_dir):
    """Generate and save a scatter plot for allele frequencies."""
    plt.figure(figsize=(10, 10))
    sns.scatterplot(data=df, x='eval_af', y='truth_af')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Eval VCF AF (log scale)')
    plt.ylabel('Truth VCF AF (log scale)')
    plt.title(f'Allele Frequency Correlation for {pop}')
    plt.grid(True, which="both", ls="--")
    
    # Calculate R^2
    r, _ = pearsonr(df['eval_af'], df['truth_af'])
    r_squared = r**2
    plt.text(0.05, 0.95, f'$R^2 = {r_squared:.4f}$', transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top')
             
    plot_path = os.path.join(output_dir, f"{pop}.png")
    plt.savefig(plot_path)
    plt.close()

def plot_vep_heatmap(df, column, output_dir):
    """Generate and save a heatmap for VEP annotation concordance."""
    concordance_matrix = pd.crosstab(df['eval_vep'], df['truth_vep'])
    
    # Create the text to display in each cell
    annot_labels = np.full(concordance_matrix.shape, "", dtype=object)
    for i, idx in enumerate(concordance_matrix.index):
        for j, col in enumerate(concordance_matrix.columns):
            count = concordance_matrix.loc[idx, col]
            if count > 0:
                total_eval = concordance_matrix.loc[idx, :].sum()
                perc = (count / total_eval) * 100 if total_eval > 0 else 0
                annot_labels[i, j] = f"{count}/{total_eval}\n({perc:.1f}%)"

    plt.figure(figsize=(max(12, len(concordance_matrix.columns)), max(10, len(concordance_matrix.index))))
    sns.heatmap(concordance_matrix, annot=annot_labels, fmt="s", cmap="viridis", cbar=True)
    plt.title(f'VEP Concordance for {column}')
    plt.xlabel('Truth VCF Annotation')
    plt.ylabel('Eval VCF Annotation')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    plot_path = os.path.join(output_dir, f"{column}.png")
    plt.savefig(plot_path)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Annotate VCF and/or run benchmarking.")
    parser.add_argument("vcf_unmatched_from_truvari", help="Unmatched VCF from Truvari step.")
    parser.add_argument("closest_bed", help="BED file with closest matches from bedtools.")
    parser.add_argument("prefix", help="Prefix for output files.")

    # Arguments for benchmarking
    parser.add_argument("--create_benchmarks", action="store_true", help="Create benchmark plots and summaries.")
    parser.add_argument("--vcf_truth_snv", help="Original truth VCF with SNVs.")
    parser.add_argument("--vcf_truth_sv", help="Original truth VCF with SVs.")
    parser.add_argument("--exact_matched_vcf", help="VCF of variants matched exactly.")
    parser.add_argument("--truvari_matched_vcf", help="VCF of variants matched by Truvari.")
    parser.add_argument("--contig", help="Contig being processed.")
    args = parser.parse_args()

    # --- Part 1: Annotate Bedtools Matches ---
    bedtools_matches = defaultdict(list)
    try:
        closest_df = pd.read_csv(args.closest_bed, sep='\t')
        for _, row in closest_df.iterrows():
            bedtools_matches[row['query_svid']].append(row['svid_b'])
    except (pd.errors.EmptyDataError, FileNotFoundError):
        pass

    bedtools_matched_path = f"{args.prefix}.bedtools_matched.vcf.gz"
    final_unmatched_path = f"{args.prefix}.final_unmatched.vcf.gz"

    vcf_in_unmatched = pysam.VariantFile(args.vcf_unmatched_from_truvari)
    vcf_in_unmatched.header.info.add('gnomAD_V4_match', '1', 'String', 'Matching status against gnomAD v4.')

    with bgzip.BGZFile(bedtools_matched_path, "w") as matched_out, bgzip.BGZFile(final_unmatched_path, "w") as unmatched_out:
        matched_out.write(str(vcf_in_unmatched.header).encode())
        unmatched_out.write(str(vcf_in_unmatched.header).encode())
        for record in vcf_in_unmatched:
            if record.id in bedtools_matches:
                record.info['gnomAD_V4_match'] = 'BEDTOOLS_CLOSEST'
                # Note: A single eval variant could match multiple truth variants in bedtools.
                # Here we just take the first one reported by the R script.
                record.info['gnomAD_V4_match_ID'] = bedtools_matches[record.id][0]
                matched_out.write(str(record).encode())
            else:
                unmatched_out.write(str(record).encode())
    
    pysam.tabix_index(bedtools_matched_path, preset="vcf", force=True)
    pysam.tabix_index(final_unmatched_path, preset="vcf", force=True)

    # --- Part 2: Combine all parts into a final annotated VCF for the contig ---
    final_vcf_path = f"{args.prefix}.final_annotated.vcf.gz"
    
    # Check if optional VCFs exist and are not empty
    vcf_list = []
    if os.path.exists(args.exact_matched_vcf) and os.path.getsize(args.exact_matched_vcf) > 0:
        vcf_list.append(args.exact_matched_vcf)
    if os.path.exists(args.truvari_matched_vcf) and os.path.getsize(args.truvari_matched_vcf) > 0:
        vcf_list.append(args.truvari_matched_vcf)
    if os.path.exists(bedtools_matched_path) and os.path.getsize(bedtools_matched_path) > 0:
        vcf_list.append(bedtools_matched_path)
    if os.path.exists(final_unmatched_path) and os.path.getsize(final_unmatched_path) > 0:
        vcf_list.append(final_unmatched_path)
    
    if vcf_list:
        subprocess.run(["bcftools", "concat", "-a", "-Oz", "-o", final_vcf_path] + vcf_list, check=True)
        pysam.tabix_index(final_vcf_path, preset="vcf", force=True)
    else:
        # Create an empty gzipped VCF if no variants are present
        with bgzip.BGZFile(final_vcf_path, "w") as f:
            f.write(str(vcf_in_unmatched.header).encode())
        pysam.tabix_index(final_vcf_path, preset="vcf", force=True)


    if not args.create_benchmarks:
        return

    # --- Part 3: Benchmarking ---
    # Setup directories
    output_basedir = f"{args.prefix}_benchmark_results"
    af_plot_dir = os.path.join(output_basedir, "AF_plots", args.contig)
    vep_plot_dir = os.path.join(output_basedir, "VEP_plots", args.contig)
    os.makedirs(af_plot_dir, exist_ok=True)
    os.makedirs(vep_plot_dir, exist_ok=True)

    # Load all truth variants into memory for quick lookup
    truth_variants = {}
    for vcf_path in [args.vcf_truth_snv, args.vcf_truth_sv]:
        with pysam.VariantFile(vcf_path) as vcf_in:
            for record in vcf_in:
                # Key by ID for SVs/Truvari/Bedtools, and by coordinate for exact matches
                truth_variants[record.id] = record.copy()
                exact_key = (record.chrom, record.pos, record.ref, tuple(record.alts))
                truth_variants[exact_key] = record.copy()

    # Process all matched VCFs to build benchmark data
    matched_data = []
    # 1. Exact matches
    with pysam.VariantFile(args.exact_matched_vcf) as vcf_in:
        for record in vcf_in:
            key = (record.chrom, record.pos, record.ref, tuple(record.alts))
            if key in truth_variants:
                matched_data.append({'eval': record.copy(), 'truth': truth_variants[key]})
    
    # 2. Truvari matches
    with pysam.VariantFile(args.truvari_matched_vcf) as vcf_in:
        for record in vcf_in:
            if 'MatchId' in record.info and record.info['MatchId'] in truth_variants:
                 matched_data.append({'eval': record.copy(), 'truth': truth_variants[record.info['MatchId']]})

    # 3. Bedtools matches
    with pysam.VariantFile(bedtools_matched_path) as vcf_in:
        for record in vcf_in:
            if 'gnomAD_V4_match_ID' in record.info and record.info['gnomAD_V4_match_ID'] in truth_variants:
                matched_data.append({'eval': record.copy(), 'truth': truth_variants[record.info['gnomAD_V4_match_ID']]})

    if not matched_data:
        print(f"No matched variants found for contig {args.contig}. Skipping benchmark generation.")
        # Create empty summary and tarball to satisfy outputs
        summary_path = os.path.join(output_basedir, "summary.txt")
        with open(summary_path, "w") as f:
            f.write(f"Contig: {args.contig}\n")
            f.write("Total matched variants: 0\n")
        
        with tarfile.open(f"{args.prefix}.benchmarks.tar.gz", "w:gz") as tar:
            tar.add(output_basedir, arcname=os.path.basename(output_basedir))
        return

    # --- AF Benchmarking ---
    af_data = defaultdict(list)
    for item in matched_data:
        eval_info = item['eval'].info
        truth_info = item['truth'].info
        
        eval_af_fields = {normalize_af_field(k): (k, v) for k, v in eval_info.items() if is_af_field(k)}
        truth_af_fields = {normalize_af_field(k): v for k, v in truth_info.items() if is_af_field(k)}

        for norm_pop, (eval_key, eval_val) in eval_af_fields.items():
            if norm_pop in truth_af_fields:
                try:
                    eval_af = float(eval_val[0] if isinstance(eval_val, tuple) else eval_val)
                    truth_af = float(truth_af_fields[norm_pop][0] if isinstance(truth_af_fields[norm_pop], tuple) else truth_af_fields[norm_pop])
                    if eval_af > 0 and truth_af > 0: # Required for log scale
                        af_data[eval_key].append({'eval_af': eval_af, 'truth_af': truth_af})
                except (ValueError, TypeError):
                    continue
    
    for pop, data_points in af_data.items():
        df = pd.DataFrame(data_points)
        if len(df) > 1:
            plot_af_correlation(df, pop, af_plot_dir)

    # --- VEP Benchmarking ---
    eval_vep_format = get_vep_format(args.exact_matched_vcf) # Assume header is consistent
    truth_vep_format = get_vep_format(args.vcf_truth_snv)
    
    vep_key_eval = 'CSQ' if 'CSQ' in eval_vep_format else 'VEP'
    vep_key_truth = 'vep' if 'vep' in truth_vep_format else 'CSQ'
    
    if eval_vep_format and truth_vep_format:
        vep_df_data = []
        for item in matched_data:
            eval_annos = get_vep_annotations(item['eval'].info, vep_key_eval, eval_vep_format)
            truth_annos = get_vep_annotations(item['truth'].info, vep_key_truth, truth_vep_format)
            
            eval_consequences = {a.get('Consequence') for a in eval_annos if a.get('Consequence')}
            truth_consequences = {a.get('Consequence') for a in truth_annos if a.get('Consequence')}

            # For simplicity, take the first one if multiple exist
            eval_vep = next(iter(eval_consequences), 'N/A')
            truth_vep = next(iter(truth_consequences), 'N/A')
            vep_df_data.append({'eval_vep': eval_vep, 'truth_vep': truth_vep})

        vep_df = pd.DataFrame(vep_df_data)
        if not vep_df.empty:
            plot_vep_heatmap(vep_df, 'Consequence', vep_plot_dir)

    # --- Summary File ---
    total_eval_variants = 0
    with pysam.VariantFile(args.exact_matched_vcf) as f: # just need header
      # This is a hacky way to count original variants. Should be passed in.
      # For now, we'll just report on matched counts.
      pass
    
    summary_path = os.path.join(output_basedir, "summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Contig: {args.contig}\n")
        f.write(f"Total matched variants: {len(matched_data)}\n")
        # Add more summary stats here if needed

    # --- Create final tarball ---
    with tarfile.open(f"{args.prefix}.benchmarks.tar.gz", "w:gz") as tar:
        tar.add(output_basedir, arcname=os.path.basename(output_basedir))

if __name__ == "__main__":
    main() 