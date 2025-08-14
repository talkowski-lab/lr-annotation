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
    if len(df) > 1:
        r, _ = pearsonr(df['eval_af'], df['truth_af'])
        plt.text(0.05, 0.95, f'$R^2 = {r**2:.4f}$', transform=plt.gca().transAxes,
                 fontsize=12, verticalalignment='top')
    plt.savefig(os.path.join(output_dir, f"{pop}.png"))
    plt.close()


def write_vep_table(eval_values, truth_values, column, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    eval_labels = sorted(list(set(eval_values)))
    truth_labels = sorted(list(set(truth_values)))
    if not eval_labels:
        return
    concordance = pd.DataFrame(0, index=truth_labels, columns=eval_labels)
    for e, t in zip(eval_values, truth_values):
        if e in eval_labels and t in truth_labels:
            concordance.loc[t, e] += 1
    concordance.to_csv(os.path.join(output_dir, f"{column}.tsv"), sep='\t')


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

    # AF aggregation
    af_groups = {}
    for p in af_pair_paths:
        with gzip.open(p, 'rt') as f:
            df = pd.read_csv(f, sep='\t')
        if df.empty or 'af_key' not in df.columns:
            continue
        for key, sub in df.groupby('af_key'):
            if key not in af_groups:
                af_groups[key] = []
            af_groups[key].append(sub[['eval_af', 'truth_af']])
    for key, pieces in af_groups.items():
        df_all = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=['eval_af','truth_af'])
        plot_af_correlation(df_all, key, af_dir)

    # VEP aggregation
    with open(args.truth_vep_header, 'r') as f:
        _, truth_vep_fields = parse_vep_header_line(f.readline())

    vep_map = {}
    for p in vep_pair_paths:
        with gzip.open(p, 'rt') as f:
            df = pd.read_csv(f, sep='\t')
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

    # reuse simpler heatmap: generate tables only here; plots can be built later from tables if needed
    for cat, data in vep_map.items():
        write_vep_table(data['eval'], data['truth'], cat.replace('/', '_'), vep_table_dir)

    # package
    with tarfile.open(f"{args.prefix}.benchmarks.tar.gz", "w:gz") as tar:
        tar.add(out_base, arcname=os.path.basename(out_base))


if __name__ == '__main__':
    main() 