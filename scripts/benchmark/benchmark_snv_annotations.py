import argparse
import os
import tarfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

def parse_info(info_str):
    """Parse INFO string into a dictionary."""
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    return info_dict

def get_vep_annotations(info_dict):
    """Extract VEP annotations from INFO dictionary."""
    vep_string = info_dict.get('vep', '')
    if not vep_string:
        return set()
    
    annotations = set()
    for annotation_blob in vep_string.split(','):
        fields = annotation_blob.split('|')
        # Assuming the consequence is the second field
        if len(fields) > 1 and fields[1]:
             annotations.add(fields[1])
    return annotations

def benchmark_snvs(eval_tsv_path, truth_tsv_path, contig, output_prefix):
    """
    Benchmark SNV annotations between eval and truth VCFs.
    """
    eval_df = pd.read_csv(eval_tsv_path, sep='\t')
    truth_df = pd.read_csv(truth_tsv_path, sep='\t')

    # Merge based on CHROM, POS, REF, ALT
    merged_df = pd.merge(eval_df, truth_df, on=['CHROM', 'POS', 'REF', 'ALT'], suffixes=('_eval', '_truth'))

    # 1. Count and percentage of matched variants
    eval_total = len(eval_df)
    matched_count = len(merged_df)
    match_percentage = (matched_count / eval_total) * 100 if eval_total > 0 else 0

    summary_stats = {
        'contig': contig,
        'total_eval_variants': eval_total,
        'matched_variants': matched_count,
        'match_percentage': match_percentage
    }

    # Create directories for plots
    af_plot_dir = f"{output_prefix}_plots/AF_plots/{contig}"
    vep_plot_dir = f"{output_prefix}_plots/VEP_plots/{contig}"
    os.makedirs(af_plot_dir, exist_ok=True)
    os.makedirs(vep_plot_dir, exist_ok=True)

    # 2. AF correlation plots
    eval_info_dicts = merged_df['INFO_eval'].apply(parse_info)
    truth_info_dicts = merged_df['INFO_truth'].apply(parse_info)
    
    af_fields = sorted({key for info in eval_info_dicts for key in info if key.startswith('AF_')})

    for af_field in af_fields:
        truth_af_field = af_field.replace('_MALE', '_XY').replace('_FEMALE', '_XX')
        
        eval_values = eval_info_dicts.apply(lambda x: float(x.get(af_field, 0)))
        truth_values = truth_info_dicts.apply(lambda x: float(x.get(truth_af_field, 0)))

        if eval_values.empty or truth_values.empty:
            continue

        corr, _ = pearsonr(eval_values, truth_values)

        plt.figure()
        sns.scatterplot(x=eval_values, y=truth_values)
        sns.regplot(x=eval_values, y=truth_values, scatter=False, color='red')
        plt.title(f'{af_field} vs {truth_af_field} on {contig}')
        plt.xlabel('Eval AF')
        plt.ylabel('Truth AF')
        plt.annotate(f'R = {corr:.2f}', xy=(0.05, 0.95), xycoords='axes fraction')
        plt.savefig(os.path.join(af_plot_dir, f'{contig}_{af_field}_correlation.png'))
        plt.close()

    # 3. VEP annotation comparison
    eval_vep_sets = eval_info_dicts.apply(get_vep_annotations)
    truth_vep_sets = truth_info_dicts.apply(get_vep_annotations)

    intersect_percentages_eval_in_truth = []
    intersect_percentages_truth_in_eval = []

    for i in range(len(merged_df)):
        eval_set = eval_vep_sets.iloc[i]
        truth_set = truth_vep_sets.iloc[i]

        if eval_set:
            intersect_percentages_eval_in_truth.append(len(eval_set.intersection(truth_set)) / len(eval_set))
        else:
            intersect_percentages_eval_in_truth.append(0)
        
        if truth_set:
            intersect_percentages_truth_in_eval.append(len(truth_set.intersection(eval_set)) / len(truth_set))
        else:
             intersect_percentages_truth_in_eval.append(0)

    plt.figure()
    sns.scatterplot(x=intersect_percentages_eval_in_truth, y=intersect_percentages_truth_in_eval)
    plt.title(f'VEP Annotation Concordance on {contig}')
    plt.xlabel('% VEP Eval in Truth')
    plt.ylabel('% VEP Truth in Eval')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.savefig(os.path.join(vep_plot_dir, f'{contig}_vep_concordance.png'))
    plt.close()

    # Write summary stats to file
    with open(f"{output_prefix}_summary.txt", 'w') as f:
        for key, value in summary_stats.items():
            f.write(f"{key}: {value}\n")

    # Create tarball of plots
    with tarfile.open(f"{output_prefix}_plots.tar.gz", "w:gz") as tar:
        tar.add(f"{output_prefix}_plots", arcname=os.path.basename(f"{output_prefix}_plots"))
        
def main():
    parser = argparse.ArgumentParser(description="Benchmark SNV annotations.")
    parser.add_argument("eval_tsv", help="Path to the eval TSV file.")
    parser.add_argument("truth_tsv", help="Path to the truth TSV file.")
    parser.add_argument("contig", help="The contig being analyzed.")
    parser.add_argument("output_prefix", help="Prefix for output files.")
    args = parser.parse_args()

    benchmark_snvs(args.eval_tsv, args.truth_tsv, args.contig, args.output_prefix)

if __name__ == "__main__":
    main() 