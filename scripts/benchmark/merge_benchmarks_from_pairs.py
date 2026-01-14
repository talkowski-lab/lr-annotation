#!/usr/bin/env python3

import argparse
import gzip
import os
import tarfile
import gc
from typing import List, Tuple, Dict

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
import numpy as np
from collections import defaultdict


def parse_vep_header_line(header_line: str) -> Tuple[str, List[str]]:
    line = header_line.strip()
    vep_id = "CSQ" if "ID=CSQ" in line else "VEP"
    fmt_part = line.split("Format:")[-1].strip().strip('"').lower()
    fields = [f.strip() for f in fmt_part.split("|")]
    return vep_id, fields


def plot_af_correlation(df: pd.DataFrame, pop: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    plt.figure(figsize=(10, 10))
    sns.scatterplot(data=df, x="eval_af", y="truth_af")
    plt.xlabel("Eval VCF AF")
    plt.ylabel("Truth VCF AF")
    plt.title(pop)
    plt.grid(True, which="both", ls="--")
    if len(df) > 1 and df["eval_af"].nunique() > 1 and df["truth_af"].nunique() > 1:
        r, _ = pearsonr(df["eval_af"], df["truth_af"])
        plt.text(
            0.05,
            0.95,
            f"$R^2 = {r**2:.4f}$",
            transform=plt.gca().transAxes,
            fontsize=12,
            verticalalignment="top",
        )
    plt.savefig(os.path.join(output_dir, f"{pop}.png"))
    plt.close()


def write_vep_table_agg(pairs_df: pd.DataFrame, column: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    eval_labels = sorted(pairs_df["eval"].unique().tolist())
    truth_labels = sorted(pairs_df["truth"].unique().tolist())
    table = pd.DataFrame(0, index=truth_labels, columns=eval_labels)
    for _, row in pairs_df.iterrows():
        table.loc[row["truth"], row["eval"]] += int(row["count"])
    table.to_csv(os.path.join(output_dir, f"{column}.tsv"), sep="\t")


def plot_vep_heatmap_agg(pairs_df: pd.DataFrame, column: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    top_n = 10
    cnts = pairs_df.groupby("eval")["count"].sum().sort_values(ascending=False)
    if len(cnts) > top_n:
        top_labels = cnts.index.tolist()[:top_n]
        col_labels = top_labels + ["Other"]
        row_labels = col_labels
    else:
        all_vals = sorted(
            list(set(pairs_df["eval"]).union(set(pairs_df["truth"]))),
            key=lambda x: cnts.get(x, 0),
            reverse=True,
        )
        col_labels = all_vals
        row_labels = all_vals

    concordance = pd.DataFrame(0, index=row_labels, columns=col_labels)
    for _, row in pairs_df.iterrows():
        e = (
            row["eval"]
            if row["eval"] in col_labels and row["eval"] != "Other"
            else "Other"
        )
        t = (
            row["truth"]
            if row["truth"] in row_labels and row["truth"] != "Other"
            else "Other"
        )
        if e in concordance.columns and t in concordance.index:
            concordance.loc[t, e] += int(row["count"])

    annot_matrix = pd.DataFrame("", index=row_labels, columns=col_labels)
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

    masked_percent_matrix = np.ma.masked_where(
        zero_count_mask.values, percent_matrix.values
    )

    sns.heatmap(
        masked_percent_matrix,
        annot=annot_matrix,
        fmt="s",
        cmap="viridis",
        cbar=True,
        cbar_kws={"label": "Concordance (%)"},
        linewidths=0.5,
        linecolor="black",
        xticklabels=col_labels,
        yticklabels=row_labels,
    )

    ax = plt.gca()
    for i, r in enumerate(row_labels):
        for j, c in enumerate(col_labels):
            if zero_count_mask.loc[r, c]:
                ax.add_patch(
                    plt.Rectangle(
                        (j, i),
                        1,
                        1,
                        fill=True,
                        facecolor="white",
                        linewidth=0.5,
                        edgecolor="black",
                    )
                )
                ax.text(
                    j + 0.5,
                    i + 0.5,
                    annot_matrix.loc[r, c],
                    ha="center",
                    va="center",
                    fontsize=10,
                )

    plt.title(f"{column}")
    plt.xlabel("Eval VCF Annotation")
    plt.ylabel("Truth VCF Annotation")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout(pad=2.0)

    plot_path = os.path.join(output_dir, f"{column}.png")
    plt.savefig(plot_path, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument(
        "--af_pair_tsvs", required=True, help="Comma-separated list of gz TSVs"
    )
    ap.add_argument(
        "--vep_pair_tsvs", required=True, help="Comma-separated list of gz TSVs"
    )
    ap.add_argument("--truth_vep_header", required=True)
    ap.add_argument(
        "--skip_vep_categories", help="Comma-separated list of VEP categories to skip"
    )
    args = ap.parse_args()

    skip_categories = set()
    if args.skip_vep_categories:
        skip_categories = {
            cat.strip().lower() for cat in args.skip_vep_categories.split(",")
        }

    af_pair_paths = [p for p in args.af_pair_tsvs.split(",") if p]
    vep_pair_paths = [p for p in args.vep_pair_tsvs.split(",") if p]
    out_base = f"{args.prefix}_benchmark_results"
    af_dir = os.path.join(out_base, "AF_plots", args.contig)
    vep_plot_dir = os.path.join(out_base, "VEP_plots", args.contig)
    af_groups: Dict[str, List[pd.DataFrame]] = defaultdict(list)
    for i, p in enumerate(af_pair_paths):
        with gzip.open(p, "rt") as fh:
            chunk_count = 0
            reader = pd.read_csv(fh, sep="\t", chunksize=200000, iterator=True)
            for chunk in reader:
                chunk_count += 1
                for key, sub in chunk.groupby("af_key"):
                    af_groups[key].append(sub[["eval_af", "truth_af"]].copy())

        gc.collect()

    for key, pieces in af_groups.items():
        if not pieces:
            continue
        df_all = pd.concat(pieces, ignore_index=True)
        plot_af_correlation(df_all, key, af_dir)

        del df_all
        pieces.clear()

    af_groups.clear()
    gc.collect()

    vep_lists: Dict[str, List[pd.DataFrame]] = defaultdict(list)
    for _, p in enumerate(vep_pair_paths):
        with gzip.open(p, "rt") as fh:
            chunk_count = 0
            reader = pd.read_csv(
                fh, sep="\t", keep_default_na=False, chunksize=200000, iterator=True
            )
            for chunk in reader:
                chunk_count += 1
                for cat, sub in chunk.groupby("category"):
                    if cat.lower() in skip_categories:
                        continue
                    vep_lists[cat].append(sub[["eval", "truth", "count"]].copy())

        gc.collect()

    for cat, frames in vep_lists.items():
        if cat.lower() in skip_categories:
            continue

        if not frames:
            continue

        df_cat = pd.concat(frames, ignore_index=True)
        df_grouped = df_cat.groupby(["eval", "truth"], as_index=False)["count"].sum()

        del df_cat
        gc.collect()

        plot_vep_heatmap_agg(df_grouped, cat.replace("/", "_"), vep_plot_dir)
        write_vep_table_agg(
            df_grouped,
            cat.replace("/", "_"),
            os.path.join(out_base, "VEP_tables", args.contig),
        )

        del df_grouped
        frames.clear()
        gc.collect()

    os.makedirs(out_base, exist_ok=True)
    with tarfile.open(f"{args.prefix}.benchmarks.tar.gz", "w:gz") as tar:
        tar.add(out_base, arcname=os.path.basename(out_base))


if __name__ == "__main__":
    main()
