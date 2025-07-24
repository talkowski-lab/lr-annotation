import argparse
import os
import tarfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import re

def parse_info(info_str):
    """Parse INFO string into a dictionary."""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def get_vep_annotations(info_dict, vep_key, vep_format):
    """Extract VEP annotations from INFO dictionary based on the provided format."""
    vep_string = info_dict.get(vep_key, '')
    if not vep_string:
        return []
    
    annotations = []
    vep_fields = vep_format.split('|')
    for annotation_blob in vep_string.split(','):
        values = annotation_blob.split('|')
        annotation_dict = dict(zip(vep_fields, values))
        annotations.append(annotation_dict)
            
    return annotations

def normalize_af_field(field):
    """Normalize AF field to a common format."""
    # Convert to lowercase and replace male/female with xy/xx
    normalized = field.lower().replace('male', 'xy').replace('female', 'xx')
    # Remove 'af' and split by '_'
    parts = [part for part in normalized.split('_') if part != 'af']
    # Return a frozenset of the parts for order-insensitive comparison
    return frozenset(parts)

def benchmark_snvs(eval_tsv_path, truth_tsv_path, eval_vep_format_str, truth_vep_format_str, contig, output_prefix):
    """
    Benchmark SNV annotations between eval and truth VCFs.
    """
    eval_df = pd.read_csv(eval_tsv_path, sep='\t')
    truth_df = pd.read_csv(truth_tsv_path, sep='\t')

    merged_df = pd.merge(eval_df, truth_df, on=['CHROM', 'POS', 'REF', 'ALT'], suffixes=('_eval', '_truth'))

    eval_total = len(eval_df)
    truth_total = len(truth_df)
    matched_count = len(merged_df)
    
    summary_stats = {
        'contig': contig,
        'total_eval_variants': eval_total,
        'total_truth_variants': truth_total,
        'matched_variants': matched_count,
        'percent_eval_in_truth': (matched_count / eval_total) * 100 if eval_total > 0 else 0,
        'percent_truth_in_eval': (matched_count / truth_total) * 100 if truth_total > 0 else 0
    }
    
    # Create output directory
    output_dir = f"{output_prefix}_results"
    os.makedirs(output_dir, exist_ok=True)
    
    af_plot_dir = os.path.join(output_dir, "AF_plots")
    os.makedirs(af_plot_dir, exist_ok=True)

    if not merged_df.empty:
        eval_info_dicts = merged_df['INFO_eval'].apply(parse_info)
        truth_info_dicts = merged_df['INFO_truth'].apply(parse_info)

        eval_af_fields = {key for info in eval_info_dicts for key in info if 'af' in key.lower()}
        truth_af_fields = {key for info in truth_info_dicts for key in info if 'af' in key.lower()}

        normalized_eval_af = {normalize_af_field(f): f for f in eval_af_fields}
        normalized_truth_af = {normalize_af_field(f): f for f in truth_af_fields}

        common_af_keys = set(normalized_eval_af.keys()) & set(normalized_truth_af.keys())

        for key in common_af_keys:
            eval_field = normalized_eval_af[key]
            truth_field = normalized_truth_af[key]
            
            eval_values = eval_info_dicts.apply(lambda x: float(x.get(eval_field, 0)))
            truth_values = truth_info_dicts.apply(lambda x: float(x.get(truth_field, 0)))

            if eval_values.empty or truth_values.empty or eval_values.isnull().all() or truth_values.isnull().all():
                continue

            corr, _ = pearsonr(eval_values, truth_values)

            plt.figure()
            sns.scatterplot(x=eval_values, y=truth_values)
            sns.regplot(x=eval_values, y=truth_values, scatter=False, color='red')
            plt.title(f'{"_".join(sorted(list(key)))} on {contig}')
            plt.xlabel(f'Eval AF ({eval_field})')
            plt.ylabel(f'Truth AF ({truth_field})')
            plt.annotate(f'R = {corr:.2f}', xy=(0.05, 0.95), xycoords='axes fraction')
            plt.savefig(os.path.join(af_plot_dir, f'{"_".join(sorted(list(key)))}_correlation.png'))
            plt.close()

        eval_vep_key, eval_vep_format = eval_vep_format_str.split('|', 1)
        truth_vep_key, truth_vep_format = truth_vep_format_str.split('|', 1)

        eval_vep_annotations = eval_info_dicts.apply(lambda x: get_vep_annotations(x, eval_vep_key, eval_vep_format))
        truth_vep_annotations = truth_info_dicts.apply(lambda x: get_vep_annotations(x, truth_vep_key, truth_vep_format))

        all_vep_fields = set(eval_vep_format.split('|')) | set(truth_vep_format.split('|'))
        vep_concordance = {field: [] for field in all_vep_fields}

        for i in range(len(merged_df)):
            eval_annos = eval_vep_annotations.iloc[i]
            truth_annos = truth_vep_annotations.iloc[i]

            for field in all_vep_fields:
                eval_values = {anno.get(field) for anno in eval_annos if anno.get(field)}
                truth_values = {anno.get(field) for anno in truth_annos if anno.get(field)}
                if eval_values and truth_values:
                    is_match = 1 if eval_values == truth_values else 0
                    vep_concordance[field].append(is_match)

        vep_summary = {}
        for field, matches in vep_concordance.items():
            if matches:
                vep_summary[field] = (sum(matches) / len(matches)) * 100
        
        if vep_summary:
            vep_df = pd.DataFrame.from_dict(vep_summary, orient='index', columns=['match_percentage'])
            vep_df.to_csv(os.path.join(output_dir, 'vep_concordance.tsv'), sep='\t')


    with open(os.path.join(output_dir, "summary.txt"), 'w') as f:
        for key, value in summary_stats.items():
            f.write(f"{key}: {value}\n")

    with tarfile.open(f"{output_prefix}.tar.gz", "w:gz") as tar:
        tar.add(output_dir, arcname=os.path.basename(output_dir))

def main():
    parser = argparse.ArgumentParser(description="Benchmark SNV annotations.")
    parser.add_argument("eval_tsv", help="Path to the eval TSV file.")
    parser.add_argument("truth_tsv", help="Path to the truth TSV file.")
    parser.add_argument("eval_vep_format", help="VEP format string for the eval VCF.")
    parser.add_argument("truth_vep_format", help="VEP format string for the truth VCF.")
    parser.add_argument("contig", help="The contig being analyzed.")
    parser.add_argument("output_prefix", help="Prefix for output files.")
    args = parser.parse_args()

    benchmark_snvs(args.eval_tsv, args.truth_tsv, args.eval_vep_format, args.truth_vep_format, args.contig, args.output_prefix)

if __name__ == "__main__":
    main() 