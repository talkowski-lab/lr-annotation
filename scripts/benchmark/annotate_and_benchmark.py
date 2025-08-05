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
            if record.key == 'INFO' and record.get('ID'):
                vep_id = record.get('ID').upper()
                if vep_id in ['VEP', 'CSQ']:
                    desc = record.get('Description')
                    if 'Format: ' in desc:
                        return vep_id.lower(), desc.split('Format: ')[-1].strip('"')
    return None, None


def get_vep_annotations(info_dict, vep_key, vep_format):
    """Extract VEP annotations from INFO dictionary."""
    if not vep_format or vep_key not in info_dict:
        return []
    
    vep_value = info_dict.get(vep_key, '')
    if isinstance(vep_value, tuple):
        vep_string = vep_value[0] if vep_value else ''
    else:
        vep_string = vep_value
    
    if not vep_string:
        return []
    
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
    
    r, _ = pearsonr(df['eval_af'], df['truth_af'])
    r_squared = r**2
    plt.text(0.05, 0.95, f'$R^2 = {r_squared:.4f}$', transform=plt.gca().transAxes,
             fontsize=12, verticalalignment='top')
             
    plot_path = os.path.join(output_dir, f"{pop}.png")
    plt.savefig(plot_path)
    plt.close()


def plot_vep_heatmap(eval_values, truth_values, column, output_dir):
    """Generate and save a heatmap for VEP annotation concordance."""
    # Count occurrences
    eval_counts = defaultdict(int)
    matched_counts = defaultdict(int)
    
    # Case insensitive matching
    eval_values_lower = [v.lower() if v else '' for v in eval_values]
    truth_values_lower = [v.lower() if v else '' for v in truth_values]
    
    for eval_val, truth_val in zip(eval_values_lower, truth_values_lower):
        if eval_val:  # Only count non-empty values
            eval_counts[eval_val] += 1
            if eval_val == truth_val:
                matched_counts[eval_val] += 1
    
    if not eval_counts:
        print(f"No data for VEP category {column}, skipping plot")
        return
    
    # Create data for heatmap
    eval_labels = sorted(eval_counts.keys())
    truth_labels = sorted(set(truth_values_lower) - {''})
    
    if not truth_labels:
        print(f"No truth data for VEP category {column}, skipping plot")
        return
    
    # Create matrix
    matrix_data = []
    annotations = []
    
    for truth_val in truth_labels:
        row_data = []
        row_annots = []
        for eval_val in eval_labels:
            if eval_val == truth_val:
                match_count = matched_counts[eval_val]
                total_count = eval_counts[eval_val]
                percentage = (match_count / total_count) * 100 if total_count > 0 else 0
                row_data.append(match_count)
                row_annots.append(f"{match_count}/{total_count}\n({percentage:.1f}%)")
            else:
                row_data.append(0)
                row_annots.append("")
        matrix_data.append(row_data)
        annotations.append(row_annots)
    
    if not matrix_data or not any(any(row) for row in matrix_data):
        print(f"No matching data for VEP category {column}, skipping plot")
        return
    
    # Create plot
    plt.figure(figsize=(max(12, len(eval_labels)), max(10, len(truth_labels))))
    
    sns.heatmap(matrix_data, 
                xticklabels=eval_labels,
                yticklabels=truth_labels,
                annot=annotations,
                fmt="s", 
                cmap="viridis", 
                cbar=True)
    
    plt.title(f'VEP Concordance for {column}')
    plt.xlabel('Eval VCF Annotation')
    plt.ylabel('Truth VCF Annotation')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    plot_path = os.path.join(output_dir, f"{column}.png")
    plt.savefig(plot_path)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Annotate VCF and/or run benchmarking.")
    parser.add_argument("--prefix", help="Prefix for output files.")
    parser.add_argument("--contig", help="Contig being processed.")
    parser.add_argument("--exact_matched_vcf", help="VCF of variants matched exactly.")
    parser.add_argument("--truvari_matched_vcf", help="VCF of variants matched by Truvari.")
    parser.add_argument("--truvari_too_small_vcf", help="VCF of variants too small for Truvari.")
    parser.add_argument("--truvari_unmatched_vcf", help="Unmatched VCF from Truvari step.")
    parser.add_argument("--closest_bed", help="BED file with closest matches from bedtools.")
    parser.add_argument("--vcf_truth_snv", help="Original truth VCF with SNVs.")
    parser.add_argument("--vcf_truth_sv", help="Original truth VCF with SVs.")
    args = parser.parse_args()

    # --- Part 1: Annotate Bedtools Matches and Create Final VCF ---
    bedtools_matches = defaultdict(list)
    with open(args.closest_bed, 'r') as f:
        for line in f:
            if line.startswith('query_svid'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                query_svid = fields[0]
                ref_svid = fields[1]
                if ref_svid != '.':
                    bedtools_matches[query_svid].append(ref_svid)

    # Process unmatched variants from Truvari and annotate bedtools matches
    bedtools_matched_tmp_path = f"{args.prefix}.bedtools_matched.tmp.vcf"
    final_unmatched_tmp_path = f"{args.prefix}.final_unmatched.tmp.vcf"
    bedtools_matched_out_path = f"{args.prefix}.bedtools_matched.vcf.gz"
    final_unmatched_out_path = f"{args.prefix}.final_unmatched.vcf.gz"

    vcf_in_unmatched = pysam.VariantFile(args.truvari_unmatched_vcf)
    if 'gnomAD_V4_match' not in vcf_in_unmatched.header.info:
        vcf_in_unmatched.header.info.add('gnomAD_V4_match', '1', 'String', 'Matching status against gnomAD v4.')
    if 'gnomAD_V4_match_ID' not in vcf_in_unmatched.header.info:
        vcf_in_unmatched.header.info.add('gnomAD_V4_match_ID', '1', 'String', 'Matching variant ID from gnomAD v4.')

    with open(bedtools_matched_tmp_path, "w") as matched_out, open(final_unmatched_tmp_path, "w") as unmatched_out:
        matched_out.write(str(vcf_in_unmatched.header))
        unmatched_out.write(str(vcf_in_unmatched.header))
        for record in vcf_in_unmatched:
            if record.id in bedtools_matches:
                record.info['gnomAD_V4_match'] = 'BEDTOOLS_CLOSEST'
                record.info['gnomAD_V4_match_ID'] = bedtools_matches[record.id][0]
                matched_out.write(str(record))
            else:
                unmatched_out.write(str(record))

    # Convert to compressed VCFs
    subprocess.run(["bcftools", "view", "-Oz", "-o", bedtools_matched_out_path, bedtools_matched_tmp_path], check=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", bedtools_matched_out_path], check=True)
    os.remove(bedtools_matched_tmp_path)

    subprocess.run(["bcftools", "view", "-Oz", "-o", final_unmatched_out_path, final_unmatched_tmp_path], check=True)
    subprocess.run(["tabix", "-p", "vcf", "-f", final_unmatched_out_path], check=True)
    os.remove(final_unmatched_tmp_path)

    # Build final annotated VCF by concatenating all parts
    vcf_files_to_concat = []
    if args.exact_matched_vcf and os.path.exists(args.exact_matched_vcf):
        vcf_files_to_concat.append(args.exact_matched_vcf)
    if args.truvari_matched_vcf and os.path.exists(args.truvari_matched_vcf):
        vcf_files_to_concat.append(args.truvari_matched_vcf)
    if os.path.exists(bedtools_matched_out_path):
        vcf_files_to_concat.append(bedtools_matched_out_path)
    if os.path.exists(final_unmatched_out_path):
        vcf_files_to_concat.append(final_unmatched_out_path)
    if args.truvari_too_small_vcf and os.path.exists(args.truvari_too_small_vcf):
        vcf_files_to_concat.append(args.truvari_too_small_vcf)

    final_vcf_path = f"{args.prefix}.final_annotated.vcf.gz"
    if vcf_files_to_concat:
        concat_cmd = ["bcftools", "concat", "-a", "-Oz", "-o", final_vcf_path] + vcf_files_to_concat
        subprocess.run(concat_cmd, check=True)
        subprocess.run(["tabix", "-p", "vcf", "-f", final_vcf_path], check=True)

    # --- Part 2: Benchmarking using Final VCF ---
    # Setup directories
    output_basedir = f"{args.prefix}_benchmark_results"
    af_plot_dir = os.path.join(output_basedir, "AF_plots", args.contig)
    vep_plot_dir = os.path.join(output_basedir, "VEP_plots", args.contig)
    os.makedirs(af_plot_dir, exist_ok=True)
    os.makedirs(vep_plot_dir, exist_ok=True)

    # Load all truth variants into memory for quick lookup
    truth_variants = {}
    for vcf_path in [args.vcf_truth_snv, args.vcf_truth_sv]:
        if vcf_path and os.path.exists(vcf_path):
            with pysam.VariantFile(vcf_path) as vcf_in:
                for record in vcf_in:
                    truth_variants[record.id] = record.copy()
                    exact_key = (record.chrom, record.pos, record.ref, tuple(record.alts))
                    truth_variants[exact_key] = record.copy()

    # Process final VCF to build benchmark data
    matched_data = []
    with pysam.VariantFile(final_vcf_path) as vcf_in:
        for record in vcf_in:
            if 'gnomAD_V4_match_ID' in record.info:
                match_id = record.info['gnomAD_V4_match_ID']
                if match_id in truth_variants:
                    matched_data.append({'eval': record.copy(), 'truth': truth_variants[match_id]})

    # --- AF Benchmarking ---
    af_data = defaultdict(list)
    for item in matched_data:
        eval_info = item['eval'].info
        truth_info = item['truth'].info
        
        eval_af_fields = {normalize_af_field(k): (k, v) for k, v in eval_info.items() if is_af_field(k)}
        truth_af_fields = {normalize_af_field(k): v for k, v in truth_info.items() if is_af_field(k)}

        for norm_pop, (eval_key, eval_val) in eval_af_fields.items():
            if norm_pop in truth_af_fields:
                eval_af = float(eval_val[0] if isinstance(eval_val, tuple) else eval_val)
                truth_af = float(truth_af_fields[norm_pop][0] if isinstance(truth_af_fields[norm_pop], tuple) else truth_af_fields[norm_pop])
                if eval_af > 0 and truth_af > 0:
                    af_data[eval_key].append({'eval_af': eval_af, 'truth_af': truth_af})
    
    for pop, data_points in af_data.items():
        df = pd.DataFrame(data_points)
        if len(df) > 1:
            plot_af_correlation(df, pop, af_plot_dir)

    # --- VEP Benchmarking ---
    eval_vep_key, eval_vep_format = get_vep_format(final_vcf_path)
    truth_vep_key, truth_vep_format = get_vep_format(args.vcf_truth_snv) if args.vcf_truth_snv else (None, None)
    
    if eval_vep_format and truth_vep_format:
        # Collect all VEP data
        vep_data_by_category = defaultdict(lambda: {'eval': [], 'truth': []})
        
        for item in matched_data:
            eval_annos = get_vep_annotations(item['eval'].info, eval_vep_key, eval_vep_format)
            truth_annos = get_vep_annotations(item['truth'].info, truth_vep_key, truth_vep_format)
            
            # Get all possible categories from both eval and truth formats
            eval_categories = set(eval_vep_format.split('|'))
            truth_categories = set(truth_vep_format.split('|'))
            all_categories = eval_categories.intersection(truth_categories)
            
            # For each category, collect values
            for category in all_categories:
                eval_values = []
                truth_values = []
                
                for anno in eval_annos:
                    val = anno.get(category, '')
                    if val and val.strip():
                        eval_values.append(val.strip())
                
                for anno in truth_annos:
                    val = anno.get(category, '')
                    if val and val.strip():
                        truth_values.append(val.strip())
                
                # Take first non-empty value for simplicity
                eval_val = eval_values[0] if eval_values else ''
                truth_val = truth_values[0] if truth_values else ''
                
                if eval_val or truth_val:  # Only collect if at least one has a value
                    vep_data_by_category[category]['eval'].append(eval_val)
                    vep_data_by_category[category]['truth'].append(truth_val)
        
        # Create plots for each category with data
        for category, data in vep_data_by_category.items():
            if len(data['eval']) > 0:
                plot_vep_heatmap(data['eval'], data['truth'], category, vep_plot_dir)

    # --- Summary File ---
    summary_path = os.path.join(output_basedir, "summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Contig: {args.contig}\n")
        f.write(f"Total matched variants: {len(matched_data)}\n")

    # --- Create final tarball ---
    with tarfile.open(f"{args.prefix}.benchmarks.tar.gz", "w:gz") as tar:
        tar.add(output_basedir, arcname=os.path.basename(output_basedir))

if __name__ == "__main__":
    main() 
