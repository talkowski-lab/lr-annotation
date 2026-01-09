#!/usr/bin/env python

import pandas as pd
import argparse
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

def load_stats_files(stats_files_list):
    all_data = []
    with open(stats_files_list, 'r') as f:
        stats_files = [line.strip() for line in f]
    
    for stats_file in stats_files:
        df = pd.read_csv(stats_file, sep='\t')
        all_data.append(df)
    
    return pd.concat(all_data, ignore_index=True)

def write_aggregated_match_data(output_path, df):
    df.to_csv(output_path, sep='\t', index=False)

def write_concordance_matrix(output_path, df):
    matrix = np.zeros((3, 3), dtype=int)
    for _, row in df.iterrows():
        t = int(row['trgt_non_ref_count'])
        v = int(row['vamos_non_ref_count'])
        if 0 <= t <= 2 and 0 <= v <= 2:
            matrix[t, v] += 1
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['NonRef_Haps', 'Vamos_0', 'Vamos_1', 'Vamos_2'])
        for i in range(3):
            writer.writerow([f"TRGT_{i}", matrix[i, 0], matrix[i, 1], matrix[i, 2]])

def plot_similarity_across_callers(output_path, df):
    if len(df) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(12, 6))
    bins = np.linspace(0, 1, 101)
    
    divergent = df[df['category'] == 'DIVERGENT']['similarity'].values
    similar = df[df['category'] == 'SIMILAR']['similarity'].values
    exact = df[df['category'] == 'EXACT']['similarity'].values
    
    total = len(df)
    div_pct = len(divergent) / total * 100
    sim_pct = len(similar) / total * 100
    exact_pct = len(exact) / total * 100
    
    ax.hist([divergent, similar, exact], bins=bins, stacked=True, 
            label=[f'< 90% Similar ({div_pct:.1f}%)', 
                   f'≥ 90% Similar ({sim_pct:.1f}%)', 
                   f'Identical ({exact_pct:.1f}%)'],
            color=['#e74c3c', '#3498db', '#2ecc71'], 
            edgecolor='black', linewidth=0.3)
    
    ax.set_xlabel('Similarity Score (1 - Edit Distance / Locus Size)', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Similarity Across Callers', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_edit_distance_across_callers(output_path, df):
    if len(df) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    divergent = df[df['category'] == 'DIVERGENT']['edit_dist'].values
    similar = df[df['category'] == 'SIMILAR']['edit_dist'].values
    exact = df[df['category'] == 'EXACT']['edit_dist'].values
    
    all_dists = np.concatenate([divergent, similar, exact])
    bins = np.arange(0, min(max(all_dists) + 2, 102))
    
    total = len(df)
    div_pct = len(divergent) / total * 100
    sim_pct = len(similar) / total * 100
    exact_pct = len(exact) / total * 100
    
    ax.hist([divergent, similar, exact], bins=bins, stacked=True,
            label=[f'< 90% Similar ({div_pct:.1f}%)', 
                   f'≥ 90% Similar ({sim_pct:.1f}%)', 
                   f'Identical ({exact_pct:.1f}%)'],
            color=['#e74c3c', '#3498db', '#2ecc71'], 
            edgecolor='black', linewidth=0.3)
    
    ax.set_xlabel('Summed Edit Distance', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('Edit Distance Across Callers', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_length_difference_across_callers(output_path, df):
    if len(df) == 0:
        return
    
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    size_ranges = [(1, 6, '1-6bp'), (7, 20, '7-20bp'), (21, float('inf'), '≥21bp')]
    
    def get_size_bin(size):
        for i, (lo, hi, _) in enumerate(size_ranges):
            if lo <= size <= hi:
                return i
        return -1
    
    trgt_by_bin = [[] for _ in size_ranges]
    vamos_by_bin = [[] for _ in size_ranges]
    
    for _, row in df.iterrows():
        trgt_bin = get_size_bin(row['trgt_motif_size'])
        vamos_bin = get_size_bin(row['vamos_motif_size'])
        if trgt_bin >= 0:
            trgt_by_bin[trgt_bin].append(row['length_diff'])
        if vamos_bin >= 0:
            vamos_by_bin[vamos_bin].append(row['length_diff'])
    
    all_diffs = df['length_diff'].values
    
    def plot_length_hist(ax, diffs, title):
        if len(diffs) > 0:
            max_abs = max(abs(min(diffs)), abs(max(diffs)))
            bins = np.arange(-min(max_abs + 2, 52), min(max_abs + 3, 53))
            neg = [d for d in diffs if d < 0]
            zero = [d for d in diffs if d == 0]
            pos = [d for d in diffs if d > 0]
            ax.hist([neg, zero, pos], bins=bins, stacked=True, 
                   color=['#3498db', '#95a5a6', '#e67e22'],
                   label=[f'Vamos Shorter (n={len(neg)})', 
                          f'Same (n={len(zero)})', 
                          f'Vamos Longer (n={len(pos)})'],
                   edgecolor='black', linewidth=0.3)
            ax.axvline(np.mean(diffs), color='darkred', linestyle='--', 
                      linewidth=2, label=f'Mean: {np.mean(diffs):.2f}')
            ax.legend(fontsize=7, loc='upper right')
        ax.set_xlabel('Length Difference (Vamos - TRGT)', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.grid(axis='y', alpha=0.3)
    
    for col, (diffs, label) in enumerate([(all_diffs, 'All Motifs')] + 
                                         [(trgt_by_bin[i], size_ranges[i][2]) for i in range(3)]):
        plot_length_hist(axes[0, col], diffs, f'TRGT Motif - {label} (n={len(diffs)})')
    
    for col, (diffs, label) in enumerate([(all_diffs, 'All Motifs')] + 
                                         [(vamos_by_bin[i], size_ranges[i][2]) for i in range(3)]):
        plot_length_hist(axes[1, col], diffs, f'Vamos Motif - {label} (n={len(diffs)})')
    
    fig.suptitle('Length Difference Across Callers', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_length_difference_vs_locus_size(output_path, df):
    if len(df) == 0:
        return
    
    fig, ax = plt.subplots(figsize=(12, 6))
    diffs = df['length_diff'].values
    sizes = df['locus_size'].values
    colors = ['#e67e22' if d > 0 else '#3498db' if d < 0 else '#95a5a6' for d in diffs]
    
    ax.scatter(sizes, diffs, c=colors, alpha=0.5, s=10)
    ax.axhline(0, color='black', linestyle='-', linewidth=1)
    ax.set_xlabel('Locus Size (bp)', fontsize=12)
    ax.set_ylabel('Length Difference (Vamos - TRGT)', fontsize=12)
    ax.set_title('Length Difference Across Callers vs Locus Size', fontsize=14)
    ax.grid(alpha=0.3)
    
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor='#e67e22', label='Vamos longer'),
        Patch(facecolor='#3498db', label='Vamos shorter'),
        Patch(facecolor='#95a5a6', label='Same length')
    ], fontsize=10)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_edit_distance_to_reference(output_path, df):
    if len(df) == 0:
        return
    
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    size_ranges = [(1, 6, '1-6bp'), (7, 20, '7-20bp'), (21, float('inf'), '≥21bp')]
    
    def get_size_bin(size):
        for i, (lo, hi, _) in enumerate(size_ranges):
            if lo <= size <= hi:
                return i
        return -1
    
    trgt_by_bin = [[] for _ in size_ranges]
    vamos_by_bin = [[] for _ in size_ranges]
    
    for _, row in df.iterrows():
        trgt_bin = get_size_bin(row['trgt_motif_size'])
        vamos_bin = get_size_bin(row['vamos_motif_size'])
        if trgt_bin >= 0:
            trgt_by_bin[trgt_bin].append(row['trgt_ref_dist'])
        if vamos_bin >= 0:
            vamos_by_bin[vamos_bin].append(row['vamos_ref_dist'])
    
    trgt_all = df['trgt_ref_dist'].values
    vamos_all = df['vamos_ref_dist'].values
    
    for col, (dists, label) in enumerate([(trgt_all, 'All Motifs')] + 
                                         [(trgt_by_bin[i], size_ranges[i][2]) for i in range(3)]):
        ax = axes[0, col]
        if len(dists) > 0:
            max_dist = max(dists)
            bins = np.arange(0, min(max_dist + 2, 102))
            ax.hist(dists, bins=bins, color='#3498db', edgecolor='black', 
                   linewidth=0.3, alpha=0.8)
            ax.axvline(np.mean(dists), color='red', linestyle='--', 
                      label=f'Mean: {np.mean(dists):.1f}')
            ax.axvline(np.median(dists), color='orange', linestyle='--', 
                      label=f'Median: {np.median(dists):.1f}')
            ax.legend(fontsize=9)
        ax.set_xlabel('Summed Edit Distance', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(f'TRGT - {label} (n={len(dists)})', fontsize=12)
        ax.grid(axis='y', alpha=0.3)
    
    for col, (dists, label) in enumerate([(vamos_all, 'All Motifs')] + 
                                         [(vamos_by_bin[i], size_ranges[i][2]) for i in range(3)]):
        ax = axes[1, col]
        if len(dists) > 0:
            max_dist = max(dists)
            bins = np.arange(0, min(max_dist + 2, 102))
            ax.hist(dists, bins=bins, color='#e74c3c', edgecolor='black', 
                   linewidth=0.3, alpha=0.8)
            ax.axvline(np.mean(dists), color='red', linestyle='--', 
                      label=f'Mean: {np.mean(dists):.1f}')
            ax.axvline(np.median(dists), color='orange', linestyle='--', 
                      label=f'Median: {np.median(dists):.1f}')
            ax.legend(fontsize=9)
        ax.set_xlabel('Summed Edit Distance', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(f'Vamos - {label} (n={len(dists)})', fontsize=12)
        ax.grid(axis='y', alpha=0.3)
    
    fig.suptitle('Edit Distance to Reference', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_length_difference_to_reference(output_path, df):
    if len(df) == 0:
        return
    
    fig, axes = plt.subplots(2, 4, figsize=(22, 10))
    size_ranges = [(1, 6, '1-6bp'), (7, 20, '7-20bp'), (21, float('inf'), '≥21bp')]
    
    def get_size_bin(size):
        for i, (lo, hi, _) in enumerate(size_ranges):
            if lo <= size <= hi:
                return i
        return -1
    
    vamos_ref_diffs = []
    trgt_ref_diffs = []
    
    for _, row in df.iterrows():
        ref_len = row['locus_size']
        vamos_total_len = row['vamos_ref_dist']
        trgt_total_len = row['trgt_ref_dist']
        vamos_ref_diffs.append(vamos_total_len - ref_len)
        trgt_ref_diffs.append(trgt_total_len - ref_len)
    
    vamos_by_bin = [[] for _ in size_ranges]
    trgt_by_bin = [[] for _ in size_ranges]
    
    for i, row in df.iterrows():
        vamos_bin = get_size_bin(row['vamos_motif_size'])
        trgt_bin = get_size_bin(row['trgt_motif_size'])
        if vamos_bin >= 0:
            vamos_by_bin[vamos_bin].append(vamos_ref_diffs[i])
        if trgt_bin >= 0:
            trgt_by_bin[trgt_bin].append(trgt_ref_diffs[i])
    
    def plot_ref_diff_hist(ax, diffs, title):
        if len(diffs) > 0:
            max_abs = max(abs(min(diffs)), abs(max(diffs)))
            bins = np.arange(-min(max_abs + 2, 102), min(max_abs + 3, 103))
            neg = [d for d in diffs if d < 0]
            zero = [d for d in diffs if d == 0]
            pos = [d for d in diffs if d > 0]
            ax.hist([neg, zero, pos], bins=bins, stacked=True, 
                   color=['#3498db', '#95a5a6', '#e67e22'],
                   label=[f'Shorter (n={len(neg)})', 
                          f'Same (n={len(zero)})', 
                          f'Longer (n={len(pos)})'],
                   edgecolor='black', linewidth=0.3)
            ax.axvline(np.mean(diffs), color='darkred', linestyle='--', 
                      linewidth=2, label=f'Mean: {np.mean(diffs):.2f}')
            ax.legend(fontsize=7, loc='upper right')
        ax.set_xlabel('Length Difference (Caller - Reference)', fontsize=11)
        ax.set_ylabel('Frequency', fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.grid(axis='y', alpha=0.3)
    
    for col, (diffs, label) in enumerate([(vamos_ref_diffs, 'All Motifs')] + 
                                         [(vamos_by_bin[i], size_ranges[i][2]) for i in range(3)]):
        plot_ref_diff_hist(axes[0, col], diffs, f'Vamos - {label} (n={len(diffs)})')
    
    for col, (diffs, label) in enumerate([(trgt_ref_diffs, 'All Motifs')] + 
                                         [(trgt_by_bin[i], size_ranges[i][2]) for i in range(3)]):
        plot_ref_diff_hist(axes[1, col], diffs, f'TRGT - {label} (n={len(diffs)})')
    
    fig.suptitle('Length Difference to Reference', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def generate_outputs(output_prefix, df, suffix, include_ref_plots):
    write_aggregated_match_data(f"{output_prefix}.{suffix}.match_data.tsv", df)
    write_concordance_matrix(f"{output_prefix}.{suffix}.genotype_concordance_matrix.tsv", df)
    plot_similarity_across_callers(f"{output_prefix}.{suffix}.similarity_across_callers.png", df)
    plot_edit_distance_across_callers(f"{output_prefix}.{suffix}.edit_distance_across_callers.png", df)
    plot_length_difference_across_callers(f"{output_prefix}.{suffix}.length_difference_across_callers.png", df)
    plot_length_difference_vs_locus_size(f"{output_prefix}.{suffix}.length_difference_vs_locus_size.png", df)
    
    if include_ref_plots:
        plot_edit_distance_to_reference(f"{output_prefix}.{suffix}.edit_distance_to_reference.png", df)
        plot_length_difference_to_reference(f"{output_prefix}.{suffix}.length_difference_to_reference.png", df)

def main():
    parser = argparse.ArgumentParser(description='Aggregate STR benchmark statistics and generate plots')
    parser.add_argument('--stats-files', required=True, help='File containing list of per-sample stats TSVs')
    parser.add_argument('--output-prefix', required=True, help='Output prefix for all files')
    parser.add_argument('--include-all-regions', required=True, choices=['true', 'false'],
                       help='Whether to generate outputs for all regions (in addition to non-ref)')
    args = parser.parse_args()
    
    include_all = args.include_all_regions == 'true'
    
    df_all = load_stats_files(args.stats_files)
    
    if include_all:
        generate_outputs(args.output_prefix, df_all, 'all_regions', include_ref_plots=True)
    
    df_non_ref = df_all[(df_all['trgt_non_ref_count'] > 0) | (df_all['vamos_non_ref_count'] > 0)]
    generate_outputs(args.output_prefix, df_non_ref, 'non_ref_regions', include_ref_plots=False)
    
    print(f"\nProcessed {len(df_all)} total regions across {df_all['sample_id'].nunique()} samples")
    print(f"Non-ref regions: {len(df_non_ref)} ({len(df_non_ref)/len(df_all)*100:.1f}%)")
    
    for category in ['EXACT', 'SIMILAR', 'DIVERGENT']:
        count = len(df_non_ref[df_non_ref['category'] == category])
        pct = count / len(df_non_ref) * 100 if len(df_non_ref) > 0 else 0
        print(f"  {category}: {count} ({pct:.1f}%)")

if __name__ == "__main__":
    main()
