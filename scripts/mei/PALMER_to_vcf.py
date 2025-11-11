#!/usr/bin/env python3

import argparse
import sys
import numpy


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert PALMER MEI calls to VCF format."
    )
    parser.add_argument("--palmer_calls", required=True, help="PALMER calls file")
    parser.add_argument("--mei_type", required=True, help="Mobile element type (e.g., ALU, LINE1, SVA)")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--ref_fai", required=True, help="Reference genome FAI index file")
    parser.add_argument("--haplotype", required=True, help="Haplotype genotype (e.g., 1|0 or 0|1)")
    return parser.parse_args()


def write_vcf_header(ref_fai_path, sample):
    print("##fileformat=VCFv4.2")
    with open(ref_fai_path, "r") as ref_fai:
        for line in ref_fai:
            contig, length = line.strip().split("\t")[:2]
            print(f"##contig=<ID={contig},length={length}>")

    print('##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##INFO=<ID=CONF_READS,Number=.,Type=Integer,Description="Number of confident supporting reads">')
    print('##INFO=<ID=ORI,Number=.,Type=String,Description="Orientation of insertion">')
    print('##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">')
    print('##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Type of SV">')
    print('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Length of SV">')
    print('##INFO=<ID=POLYA_LEN,Number=.,Type=Integer,Description="Length of polyA tail">')
    print('##INFO=<ID=TSD_5PRIME_LEN,Number=.,Type=Integer,Description="Length of 5prime TSD">')
    print('##INFO=<ID=TSD_3PRIME_LEN,Number=.,Type=Integer,Description="Length of 3prime TSD">')
    print('##INFO=<ID=TRANSD_LEN,Number=.,Type=Integer,Description="Length of transduction">')
    print('##INFO=<ID=INVERSION_5PRIME,Number=0,Type=Flag,Description="Whether the MEI has a 5 prime inversion">')
    print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")


def parse_palmer_calls(callfile_path, min_conf=1):
    calls = {}
    seen_records = set()

    with open(callfile_path, "r") as callfile:
        for line in callfile:
            if line.startswith("cluster_id\t"):
                continue

            fields = line.strip().split("\t")
            
            # Filter by confidence threshold (column 10, 0-indexed)
            conf = int(fields[10])
            if conf < min_conf:
                continue

            # Create unique key from columns 2-10 and 14 (0-indexed: 2-10 and 14)
            unique_key = tuple(fields[2:11] + [fields[14]])
            if unique_key in seen_records:
                continue
            seen_records.add(unique_key)

            # Parse call information
            cid = fields[0]
            call = {
                "chrom": fields[1],
                "pos": round(numpy.median([int(fields[2]), int(fields[3]), int(fields[4]), int(fields[5])])),
                "conf": fields[10],
                "ori": fields[13],
                "polyAlen": fields[14],
                "5TSDlen": fields[15],
                "3TSDlen": fields[16],
                "transDlen": fields[17],
                "5inv": int(fields[18]),
            }

            # Calculate MEI length
            tmp_end = round(numpy.median([int(fields[8]), int(fields[9])]))
            tmp_start = round(numpy.median([int(fields[6]), int(fields[7])]))
            if call["5inv"]:
                tmp_5invend = int(fields[19])
                tmp_5invstart = int(fields[20])
                call["length"] = tmp_end - tmp_start + int(call["transDlen"]) + tmp_5invend - tmp_5invstart
            else:
                call["length"] = tmp_end - tmp_start + int(call["transDlen"])

            calls[cid] = call

    return calls


def write_vcf_records(calls, mei_type, haplotype):
    for cid, call in calls.items():
        info_fields = (
            f"SVTYPE=INS;ME_TYPE={mei_type};SVLEN={call['length']};"
            f"CONF_READS={call['conf']};ORI={call['ori']};POLYA_LEN={call['polyAlen']};"
            f"TSD_5PRIME_LEN={call['5TSDlen']};TSD_3PRIME_LEN={call['3TSDlen']};"
            f"TRANSD_LEN={call['transDlen']}"
        )
        if call["5inv"]:
            info_fields += ";INVERSION_5PRIME"
        
        print(f"{call['chrom']}\t{call['pos']}\t{cid}\tN\tN\t60\t.\t{info_fields}\tGT\t{haplotype}")


def main():
    args = parse_args()
    
    write_vcf_header(args.ref_fai, args.sample)
    calls = parse_palmer_calls(args.palmer_calls)
    write_vcf_records(calls, args.mei_type, args.haplotype)


if __name__ == "__main__":
    main()
