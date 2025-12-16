#!/usr/bin/env python3

import argparse
import os
import subprocess
import matplotlib.pyplot as plt

from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate genotype distribution plots from TRGT VCF files."
    )
    parser.add_argument(
        "--vcf",
        action="append",
        required=False,
        help="Path to VCF file (can be specified multiple times)"
    )
    parser.add_argument(
        "--vcf-list",
        type=str,
        required=False,
        help="File containing line-separated VCF paths (local or gs://)"
    )
    parser.add_argument(
        "--chromosomes",
        type=str,
        help="Comma-separated list of chromosomes to include (default: all)"
    )
    parser.add_argument(
        "--ancestry",
        required=False,
        help="Path to ancestry file (tab-separated: sample_name, ancestry)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for plots"
    )
    return parser.parse_args()


def load_ancestry_map(ancestry_path):
    ancestry_map = {}
    with open(ancestry_path, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 2:
                ancestry_map[fields[0]] = fields[1]
    return ancestry_map


def get_sample_name(vcf_path):
    if vcf_path.startswith("gs://"):
        cmd = [
            "bcftools", "query", "-l", f"{vcf_path}"
        ]
    else:
        cmd = ["bcftools", "query", "-l", vcf_path]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def count_genotypes_for_sample(vcf_path, chromosomes):
    gt_counts = defaultdict(int)
    if vcf_path.startswith("gs://"):
        # Stream from GCS using bcftools
        cmd = ["bcftools", "query", "-f", "[%GT]\n"]
        if chromosomes:
            cmd.extend(["-r", ",".join(chromosomes)])
        cmd.append(vcf_path)
    else:
        cmd = ["bcftools", "query", "-f", "[%GT]\n", vcf_path]
        if chromosomes:
            cmd.extend(["-r", ",".join(chromosomes)])
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in result.stdout.strip().split("\n"):
        if line:
            gt_counts[line] += 1
    return dict(gt_counts)


def get_ancestry_colors():
    return {
        "afr": "#E41A1C",
        "amr": "#377EB8",
        "asj": "#4DAF4A",
        "eas": "#984EA3",
        "fin": "#FF7F00",
        "mid": "#FFFF33",
        "nfe": "#A65628",
        "sas": "#F781BF",
        "oth": "#999999",
    }


def plot_genotype_distribution(gt_data, ancestry_map, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # Create subfolders for plots and histograms
    plots_dir = os.path.join(output_dir, "plots")
    hists_dir = os.path.join(output_dir, "histograms")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(hists_dir, exist_ok=True)

    # Write stats.txt with (sample, genotype, count) per row
    stats_path = os.path.join(output_dir, "stats.txt")
    with open(stats_path, "w") as stats_file:
        stats_file.write("sample\tgenotype\tcount\n")
        for sample_name, gts in gt_data.items():
            for gt, count in gts.items():
                stats_file.write(f"{sample_name}\t{gt}\t{count}\n")

    all_genotypes = set()
    for sample_gts in gt_data.values():
        all_genotypes.update(sample_gts.keys())

    ancestry_colors = get_ancestry_colors() if ancestry_map else None

    for gt in sorted(all_genotypes):
        samples_with_gt = []
        for sample_name, gts in gt_data.items():
            if gt in gts:
                ancestry = ancestry_map.get(sample_name, "oth") if ancestry_map else None
                samples_with_gt.append({
                    "sample": sample_name,
                    "count": gts[gt],
                    "ancestry": ancestry
                })

        if ancestry_map:
            samples_with_gt.sort(key=lambda x: (x["ancestry"], x["sample"]))
        else:
            samples_with_gt.sort(key=lambda x: x["sample"])

        # Bar plot (per-sample genotype distribution)
        fig, ax = plt.subplots(figsize=(max(10, len(samples_with_gt) * 0.5), 6))
        x_positions = range(len(samples_with_gt))
        if ancestry_map:
            colors = [ancestry_colors.get(s["ancestry"], "#999999") for s in samples_with_gt]
        else:
            colors = ["#377EB8"] * len(samples_with_gt)
        counts = [s["count"] for s in samples_with_gt]
        labels = [s["sample"] for s in samples_with_gt]
        bars = ax.bar(x_positions, counts, color=colors)
        for bar, count in zip(bars, counts):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                str(count),
                ha="center",
                va="bottom",
                fontsize=7
            )
        ax.set_xticks(x_positions)
        ax.set_xticklabels(labels, rotation=90, fontsize=8)
        ax.set_xlabel("Sample")
        ax.set_ylabel("Count")
        ax.set_title(f"Genotype {gt} Distribution")
        if ancestry_map:
            unique_ancestries = sorted(set(s["ancestry"] for s in samples_with_gt))
            legend_handles = [
                plt.Rectangle((0, 0), 1, 1, color=ancestry_colors.get(anc, "#999999"))
                for anc in unique_ancestries
            ]
            ax.legend(legend_handles, unique_ancestries, title="Ancestry", loc="upper right")
        plt.tight_layout()
        gt_filename = gt.replace("/", "_").replace(".", "missing")
        plt.savefig(os.path.join(plots_dir, f"genotype_{gt_filename}.png"), dpi=150)
        plt.close()

        # Histogram (distribution of counts per genotype)
        if counts:
            fig2, ax2 = plt.subplots(figsize=(8, 6))
            ax2.hist(counts, bins='auto', color="#377EB8", edgecolor='black')
            ax2.set_xlabel("Count")
            ax2.set_ylabel("Number of Samples")
            ax2.set_title(f"Histogram of Counts for Genotype {gt}")
            plt.tight_layout()
            plt.savefig(os.path.join(hists_dir, f"histogram_{gt_filename}.png"), dpi=150)
            plt.close()


def main():
    args = parse_args()

    chromosomes = None
    if args.chromosomes:
        chromosomes = [c.strip() for c in args.chromosomes.split(",")]

    ancestry_map = load_ancestry_map(args.ancestry) if args.ancestry else None

    # Collect all VCFs from --vcf and --vcf-list
    vcf_paths = []
    if args.vcf:
        vcf_paths.extend(args.vcf)
    if args.vcf_list:
        with open(args.vcf_list, "r") as f:
            vcf_paths.extend([line.strip() for line in f if line.strip()])
    if not vcf_paths:
        raise ValueError("No VCF files provided. Use --vcf and/or --vcf-list.")

    gt_data = {}
    for vcf_path in vcf_paths:
        sample_name = get_sample_name(vcf_path)
        gt_counts = count_genotypes_for_sample(vcf_path, chromosomes)
        gt_data[sample_name] = gt_counts

        # Delete .vcf.gz.tbi index file if it exists and was created locally
        if vcf_path.endswith(".vcf.gz"):
            tbi_path = vcf_path + ".tbi"
            if os.path.exists(tbi_path):
                try:
                    os.remove(tbi_path)
                except Exception:
                    pass

    plot_genotype_distribution(gt_data, ancestry_map, args.output)


if __name__ == "__main__":
    main()
