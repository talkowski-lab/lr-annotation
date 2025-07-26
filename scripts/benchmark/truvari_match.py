#!/usr/bin/env python3

import argparse
import pysam
import subprocess
import os
import shutil

def run_truvari(vcf_eval, vcf_truth, ref_fasta, pctseq, prefix):
    output_dir = f"{prefix}_truvari_{pctseq}"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    cmd = [
        "truvari", "bench",
        "-b", vcf_truth,
        "-c", vcf_eval,
        "-o", output_dir,
        "--reference", ref_fasta,
        "--pctseq", str(pctseq),
        "--sizemin", "5",
        "--sizefilt", "10"
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Truvari for pctseq={pctseq}:")
        print(e.stderr)
        raise
    
    return f"{output_dir}/tp-comp.vcf.gz", f"{output_dir}/fn.vcf.gz"

def main():
    parser = argparse.ArgumentParser(description="Iteratively run Truvari to find matches.")
    parser.add_argument("vcf_eval_unmatched", help="Unmatched evaluation VCF file from exact matching.")
    parser.add_argument("vcf_truth", help="Truth VCF file.")
    parser.add_argument("ref_fasta", help="Reference FASTA file.")
    parser.add_argument("prefix", help="Prefix for output files.")
    args = parser.parse_args()

    # Filter truth VCF to exclude SNVs
    truth_non_snv_path = f"{args.prefix}.truth.non_snv.vcf.gz"
    subprocess.run([
        "bcftools", "view",
        "-e", 'INFO/variant_type="snv"',
        "-Oz",
        "-o", truth_non_snv_path,
        args.vcf_truth
    ], check=True)
    subprocess.run(["tabix", "-p", "vcf", truth_non_snv_path], check=True)
    
    pctseq_passes = [0.9, 0.7, 0.5]
    remaining_eval_vcf = args.vcf_eval_unmatched
    all_matched_vcfs = []

    for pctseq in pctseq_passes:
        # The 'fn.vcf.gz' is the set of unmatched variants
        matched_vcf, remaining_eval_vcf = run_truvari(remaining_eval_vcf, truth_non_snv_path, args.ref_fasta, pctseq, f"{args.prefix}_{pctseq}")
        
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
    shutil.move(remaining_eval_vcf, final_unmatched_path)
    subprocess.run(["tabix", "-p", "vcf", final_unmatched_path], check=True)

    # Combine all truvari-matched VCFs
    if all_matched_vcfs:
        combined_matched_path = f"{args.prefix}.truvari_matched.combined.vcf.gz"
        subprocess.run(["bcftools", "concat", "-a", "-Oz", "-o", combined_matched_path] + all_matched_vcfs, check=True)
        subprocess.run(["tabix", "-p", "vcf", combined_matched_path], check=True)

if __name__ == "__main__":
    main() 