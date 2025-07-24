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
        return set()
    
    annotations = set()
    vep_fields = vep_format.split('|')
    for annotation_blob in vep_string.split(','):
        values = annotation_blob.split('|')
        annotation_dict = dict(zip(vep_fields, values))
        
        # Create a frozenset of the annotation dictionary items to make it hashable
        # Only include non-empty values
        non_empty_items = frozenset({item for item in annotation_dict.items() if item[1]})
        if non_empty_items:
            annotations.add(non_empty_items)
            
    return annotations

def normalize_af_field(field):
    """Normalize AF field to a common format."""
    field = field.replace('_MALE', '_XY').replace('_FEMALE', '_XX')
    return re.sub(r'^(AF_|_AF)', '', field)

def benchmark_snvs(eval_tsv_path, truth_tsv_path, eval_vep_format_str, truth_vep_format_str, contig, output_prefix):
    """
    Benchmark SNV annotations between eval and truth VCFs.
    """
    eval_df = pd.read_csv(eval_tsv_path, sep='\t')
    truth_df = pd.read_csv(truth_tsv_path, sep='\t')

    merged_df = pd.merge(eval_df, truth_df, on=['CHROM', 'POS', 'REF', 'ALT'], suffixes=('_eval', '_truth'))

    eval_total = len(eval_df)
    matched_count = len(merged_df)
    match_percentage = (matched_count / eval_total) * 100 if eval_total > 0 else 0

    summary_stats = {
        'contig': contig,
        'total_eval_variants': eval_total,
        'matched_variants': matched_count,
        'match_percentage': match_percentage
    }
    
    # Create plot directories
    af_plot_dir = os.path.join(output_prefix, "AF_plots", contig)
    vep_plot_dir = os.path.join(output_prefix, "VEP_plots", contig)
    os.makedirs(af_plot_dir, exist_ok=True)
    os.makedirs(vep_plot_dir, exist_ok=True)

    if not merged_df.empty:
        eval_info_dicts = merged_df['INFO_eval'].apply(parse_info)
        truth_info_dicts = merged_df['INFO_truth'].apply(parse_info)

        eval_af_fields = {key for info in eval_info_dicts for key in info if 'AF' in key}
        truth_af_fields = {key for info in truth_info_dicts for key in info if 'AF' in key}

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
            plt.title(f'{eval_field} vs {truth_field} on {contig}')
            plt.xlabel('Eval AF')
            plt.ylabel('Truth AF')
            plt.annotate(f'R = {corr:.2f}', xy=(0.05, 0.95), xycoords='axes fraction')
            plt.savefig(os.path.join(af_plot_dir, f'{key}_correlation.png'))
            plt.close()

        eval_vep_key, eval_vep_format = eval_vep_format_str.split('|', 1)
        truth_vep_key, truth_vep_format = truth_vep_format_str.split('|', 1)

        eval_vep_sets = eval_info_dicts.apply(lambda x: get_vep_annotations(x, eval_vep_key, eval_vep_format))
        truth_vep_sets = truth_info_dicts.apply(lambda x: get_vep_annotations(x, truth_vep_key, truth_vep_format))

        intersect_percentages_eval_in_truth = []
        intersect_percentages_truth_in_eval = []

        for i in range(len(merged_df)):
            eval_set = eval_vep_sets.iloc[i]
            truth_set = truth_vep_sets.iloc[i]

            intersect_percentages_eval_in_truth.append(len(eval_set.intersection(truth_set)) / len(eval_set) if eval_set else 0)
            intersect_percentages_truth_in_eval.append(len(truth_set.intersection(eval_set)) / len(truth_set) if truth_set else 0)

        if intersect_percentages_eval_in_truth and intersect_percentages_truth_in_eval:
            plt.figure()
            sns.scatterplot(x=intersect_percentages_eval_in_truth, y=intersect_percentages_truth_in_eval)
            plt.title(f'VEP Annotation Concordance on {contig}')
            plt.xlabel('% VEP Eval in Truth')
            plt.ylabel('% VEP Truth in Eval')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.savefig(os.path.join(vep_plot_dir, f'vep_concordance.png'))
            plt.close()

    with open(f"{output_prefix}_summary.txt", 'w') as f:
        for key, value in summary_stats.items():
            f.write(f"{key}: {value}\n")

    with tarfile.open(f"{output_prefix}.tar.gz", "w:gz") as tar:
        tar.add(output_prefix, arcname=os.path.basename(output_prefix))

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