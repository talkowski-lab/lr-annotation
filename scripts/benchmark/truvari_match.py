#!/usr/bin/env python3

import argparse
import pysam
import subprocess
import os

def run_truvari(vcf_eval, vcf_truth, pctseq, prefix):
    output_dir = f"{prefix}_truvari_{pctseq}"
    os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        "truvari", "bench",
        "-b", vcf_truth,
        "-c", vcf_eval,
        "-o", output_dir,
        "--pctseq", str(pctseq),
        "--passonly",
        "-r", "1000",
        "-s", "10",
    ]
    subprocess.run(cmd, check=True)
    
    return f"{output_dir}/tp-comp.vcf.gz", f"{output_dir}/fn.vcf.gz"

def main():
    parser = argparse.ArgumentParser(description="Iteratively run Truvari to find matches.")
    parser.add_argument("vcf_eval_unmatched", help="Unmatched evaluation VCF file from exact matching.")
    parser.add_argument("vcf_truth", help="Truth VCF file.")
    parser.add_argument("prefix", help="Prefix for output files.")
    args = parser.parse_args()

    # Filter truth VCF to exclude SNVs using bcftools
    truth_non_snv_path = f"{args.prefix}.truth.non_snv.vcf.gz"
    subprocess.run([
        "bcftools", "view",
        "-e", 'INFO/variant_type="snv"',
        "-Oz",
        "-o", truth_non_snv_path,
        args.vcf_truth
    ], check=True)
    subprocess.run(["tabix", "-p", "vcf", truth_non_snv_path], check=True)

    # Filter eval VCF to include only variants with abs(SVLEN) >= 10 using bcftools
    eval_svlen_ge10_path = f"{args.prefix}.eval.svlen_ge10.vcf.gz"
    subprocess.run([
        "bcftools", "view",
        "-i", 'ABS(INFO/SVLEN)>=10',
        "-Oz",
        "-o", eval_svlen_ge10_path,
        args.vcf_eval_unmatched
    ], check=True)
    subprocess.run(["tabix", "-p", "vcf", eval_svlen_ge10_path], check=True)
    
    pctseq_passes = [0.9, 0.7, 0.5]
    remaining_eval_vcf = eval_svlen_ge10_path
    all_matched_vcfs = []

    for pctseq in pctseq_passes:
        matched_vcf, remaining_eval_vcf = run_truvari(remaining_eval_vcf, truth_non_snv_path, pctseq, f"{args.prefix}_{pctseq}")
        
        # Annotate the matched VCF with the match type
        annotated_matched_tmp_path = f"{args.prefix}.truvari_matched.{pctseq}.tmp.vcf"
        annotated_matched_path = f"{args.prefix}.truvari_matched.{pctseq}.vcf.gz"
        with pysam.VariantFile(matched_vcf) as vcf_in, open(annotated_matched_tmp_path, "w") as vcf_out:
             vcf_in.header.info.add('gnomAD_V4_match', '1', 'String', 'Matching status against gnomAD v4.')
             vcf_out.write(str(vcf_in.header))
             for record in vcf_in:
                 record.info['gnomAD_V4_match'] = f'TRUVARI_{pctseq}'
                 vcf_out.write(str(record))
        subprocess.run(["bcftools", "view", "-Oz", "-o", annotated_matched_path, annotated_matched_tmp_path], check=True)
        subprocess.run(["tabix", "-p", "vcf", annotated_matched_path], check=True)
        os.remove(annotated_matched_tmp_path)
        all_matched_vcfs.append(annotated_matched_path)

    # `remaining_eval_vcf` from the last truvari run is the final unmatched set
    final_unmatched_path = f"{args.prefix}.truvari_unmatched.vcf.gz"
    os.rename(remaining_eval_vcf, final_unmatched_path)
    subprocess.run(["tabix", "-p", "vcf", final_unmatched_path], check=True)

    # Combine all truvari-matched VCFs
    if all_matched_vcfs:
        combined_matched_path = f"{args.prefix}.truvari_matched.combined.vcf.gz"
        subprocess.run(["bcftools", "concat", "-a", "-Oz", "-o", combined_matched_path] + all_matched_vcfs, check=True)
        subprocess.run(["tabix", "-p", "vcf", combined_matched_path], check=True)

if __name__ == "__main__":
    main() 