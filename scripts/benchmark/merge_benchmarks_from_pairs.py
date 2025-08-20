#!/usr/bin/env python3

import argparse
import gzip
import os
import tarfile
from typing import List, Tuple

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import numpy as np
from collections import Counter


def parse_vep_header_line(header_line: str) -> Tuple[str, List[str]]:
    line = header_line.strip()
    vep_id = 'CSQ' if 'ID=CSQ' in line else 'VEP'
    fmt_part = line.split('Format:')[-1].strip().strip('"').lower()
    fields = [f.strip() for f in fmt_part.split('|')]
    return vep_id, fields


def plot_af_correlation(df: pd.DataFrame, pop: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    plt.figure(figsize=(10, 10))
    sns.scatterplot(data=df, x='eval_af', y='truth_af')
    plt.xlabel('Eval VCF AF')
    plt.ylabel('Truth VCF AF')
    plt.title(pop)
    plt.grid(True, which='both', ls='--')
    if len(df) > 1 and df['eval_af'].nunique() > 1 and df['truth_af'].nunique() > 1:
        r, _ = pearsonr(df['eval_af'], df['truth_af'])
        plt.text(0.05, 0.95, f'$R^2 = {r**2:.4f}$', transform=plt.gca().transAxes,
                 fontsize=12, verticalalignment='top')
    plt.savefig(os.path.join(output_dir, f"{pop}.png"))
    plt.close()


def write_vep_table(eval_values, truth_values, column, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    eval_values = [('N/A' if pd.isna(v) else str(v)) for v in eval_values]
    truth_values = [('N/A' if pd.isna(v) else str(v)) for v in truth_values]
    eval_labels = sorted(list(set(eval_values)))
    truth_labels = sorted(list(set(truth_values)))
    if not eval_labels:
        return
    concordance = pd.DataFrame(0, index=truth_labels, columns=eval_labels)
    for e, t in zip(eval_values, truth_values):
        if e in eval_labels and t in truth_labels:
            concordance.loc[t, e] += 1
    concordance.to_csv(os.path.join(output_dir, f"{column}.tsv"), sep='\t')


def plot_vep_heatmap(eval_values, truth_values, column, output_dir):
    # Match reference style: top-N + 'Other', percent annotations, white zero cells
    eval_values = [('N/A' if pd.isna(v) else str(v)) for v in eval_values]
    truth_values = [('N/A' if pd.isna(v) else str(v)) for v in truth_values]

    eval_value_counts = Counter(eval_values)
    if not eval_value_counts:
        return
    if len(eval_value_counts) == 1 and 'N/A' in eval_value_counts:
        return

    top_n = 10
    if len(eval_value_counts) > top_n:
        top_labels = [val for val, count in eval_value_counts.most_common(top_n)]
        col_labels = top_labels + ['Other']
        row_labels = top_labels + ['Other']
    else:
        all_unique_values = set(eval_values) | set(truth_values)
        sorted_labels = sorted(list(all_unique_values), key=lambda x: eval_value_counts.get(x, 0), reverse=True)
        col_labels = sorted_labels
        row_labels = sorted_labels

    concordance = pd.DataFrame(0, index=row_labels, columns=col_labels)
    for e, t in zip(eval_values, truth_values):
        col = e if e in col_labels and e != 'Other' else 'Other'
        row = t if t in row_labels and t != 'Other' else 'Other'
        if col in concordance.columns and row in concordance.index:
            concordance.loc[row, col] += 1

    annot_matrix = pd.DataFrame('', index=row_labels, columns=col_labels)
    percent_matrix = pd.DataFrame(0.0, index=row_labels, columns=col_labels)
    zero_count_mask = pd.DataFrame(False, index=row_labels, columns=col_labels)
    column_totals = concordance.sum(axis=0)
    for c in col_labels:
        total = column_totals[c]
        if total > 0:
            for r in row_labels:
                count = concordance.loc[r, c]
                percentage = (count / total) * 100
                percent_matrix.loc[r, c] = percentage
                annot_matrix.loc[r, c] = f"{count}/{total}\n({percentage:.1f}%)"
                if count == 0:
                    zero_count_mask.loc[r, c] = True
        else:
            for r in row_labels:
                annot_matrix.loc[r, c] = "0/0\n(0.0%)"
                zero_count_mask.loc[r, c] = True

    fig_height = max(10, len(row_labels) * 0.6)
    fig_width = max(12, len(col_labels) * 0.8)
    plt.figure(figsize=(fig_width, fig_height))

    masked_percent_matrix = np.ma.masked_where(zero_count_mask.values, percent_matrix.values)

    sns.heatmap(
        masked_percent_matrix,
        annot=annot_matrix,
        fmt='s',
        cmap='viridis',
        cbar=True,
        cbar_kws={'label': 'Concordance (%)'},
        linewidths=0.5,
        linecolor='black'
    )

    ax = plt.gca()
    for i, r in enumerate(row_labels):
        for j, c in enumerate(col_labels):
            if zero_count_mask.loc[r, c]:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, facecolor='white',
                                           linewidth=0.5, edgecolor='black'))
                ax.text(j + 0.5, i + 0.5, annot_matrix.loc[r, c],
                        ha='center', va='center', fontsize=10)

    plt.title(f'{column}')
    plt.xlabel('Eval VCF Annotation')
    plt.ylabel('Truth VCF Annotation')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout(pad=2.0)

    plot_path = os.path.join(output_dir, f"{column}.png")
    plt.savefig(plot_path, bbox_inches='tight')
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prefix', required=True)
    ap.add_argument('--contig', required=True)
    ap.add_argument('--af_pair_tsvs', required=True, help='Comma-separated list of gz TSVs')
    ap.add_argument('--vep_pair_tsvs', required=True, help='Comma-separated list of gz TSVs')
    ap.add_argument('--truth_vep_header', required=True)
    args = ap.parse_args()

    af_pair_paths = [p for p in args.af_pair_tsvs.split(',') if p]
    vep_pair_paths = [p for p in args.vep_pair_tsvs.split(',') if p]

    out_base = f"{args.prefix}_benchmark_results"
    af_dir = os.path.join(out_base, 'AF_plots', args.contig)
    vep_plot_dir = os.path.join(out_base, 'VEP_plots', args.contig)

    af_groups = {}
    for p in af_pair_paths:
        with gzip.open(p, 'rt') as f:
            df = pd.read_csv(f, sep='\t', keep_default_na=False)
        if df.empty or 'af_key' not in df.columns:
            continue
        for key, sub in df.groupby('af_key'):
            if key not in af_groups:
                af_groups[key] = []
            af_groups[key].append(sub[['eval_af', 'truth_af']])
    for key, pieces in af_groups.items():
        df_all = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=['eval_af','truth_af'])
        plot_af_correlation(df_all, key, af_dir)

    with open(args.truth_vep_header, 'r') as f:
        _, truth_vep_fields = parse_vep_header_line(f.readline())

    vep_map = {}
    for p in vep_pair_paths:
        with gzip.open(p, 'rt') as f:
            df = pd.read_csv(f, sep='\t', keep_default_na=False)
        if df.empty or 'category' not in df.columns:
            continue
        for cat, sub in df.groupby('category'):
            if cat not in vep_map:
                vep_map[cat] = {'eval': [], 'truth': []}
            vep_map[cat]['eval'].extend(sub['eval'].tolist())
            vep_map[cat]['truth'].extend(sub['truth'].tolist())

    os.makedirs(vep_plot_dir, exist_ok=True)
    vep_table_dir = os.path.join(out_base, 'VEP_tables', args.contig)
    os.makedirs(vep_table_dir, exist_ok=True)

    for cat, data in vep_map.items():
        safe_cat = cat.replace('/', '_')
        plot_vep_heatmap(data['eval'], data['truth'], safe_cat, vep_plot_dir)
        write_vep_table(data['eval'], data['truth'], safe_cat, vep_table_dir)

    with tarfile.open(f"{args.prefix}.benchmarks.tar.gz", "w:gz") as tar:
        tar.add(out_base, arcname=os.path.basename(out_base))


if __name__ == '__main__':
    main() 