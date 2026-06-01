version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TruvariMatch {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        String prefix

        Int min_sv_length_eval
        Int min_sv_length_truth
        String length_field_eval = "SVLEN"

        File ref_fa
        File ref_fai

        String utils_docker
        
        RuntimeAttr? runtime_attr_subset_eval
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_run_truvari

        Int truvari_chunksize = 10000
    }

    call Helpers.SubsetVcfByLength as SubsetEval {
        input:
            vcf = vcf_eval,
            vcf_idx = vcf_eval_idx,
            length_field = length_field_eval,
            min_length = min_sv_length_eval,
            prefix = "~{prefix}.subset_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_eval
    }

    call Helpers.SubsetVcfByArgs as SubsetTruth {
        input:
            vcf = vcf_truth,
            vcf_idx = vcf_truth_idx,
            include_args = 'abs(ILEN) >= ~{min_sv_length_truth}',
            prefix = "~{prefix}.subset_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_truth
    }

    call RunTruvari as RunTruvari_07 {
        input:
            vcf_eval = SubsetEval.subset_vcf,
            vcf_eval_idx = SubsetEval.subset_vcf_idx,
            vcf_truth = SubsetTruth.subset_vcf,
            vcf_truth_idx = SubsetTruth.subset_vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.7,
            sizemin = 0,
            sizefilt = 0,
            chunksize = truvari_chunksize,
            tag_value = "TRUVARI_0.7",
            prefix = "~{prefix}.0.7",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari
    }

    output {
        File annotation_tsv = RunTruvari_07.annotation_tsv
        File unmatched_vcf = RunTruvari_07.unmatched_vcf
        File unmatched_vcf_idx = RunTruvari_07.unmatched_vcf_idx
    }
}

task RunTruvari {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        File ref_fa
        File ref_fai
        Float pctseq
        Int sizemin
        Int sizefilt
        Int chunksize
        String tag_value
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ -d "~{prefix}_truvari" ]; then
            rm -r "~{prefix}_truvari"
        fi

        truvari bench \
            -b ~{vcf_truth} \
            -c ~{vcf_eval} \
            -o "~{prefix}_truvari" \
            --reference ~{ref_fa} \
            --pctseq ~{pctseq} \
            --sizemin ~{sizemin} \
            --sizefilt ~{sizefilt} \
            --chunksize ~{chunksize} \
            --dup-to-ins

        
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/MatchId\n' "~{prefix}_truvari/tp-comp.vcf.gz" \
            | awk 'BEGIN{FS=OFS="\t"} {split($6,a,","); print $1 "|" a[1],$1,$2,$3,$4,$5}' \
            | LC_ALL=C sort -k1,1 > comp.mid.tsv

        bcftools query -f '%CHROM\t%ID\t%INFO/MatchId\t%INFO/AF\n' "~{prefix}_truvari/tp-base.vcf.gz" \
            | awk 'BEGIN{FS=OFS="\t"} {split($3,a,","); print $1 "|" a[1],$2,$4}' \
            | LC_ALL=C sort -k1,1 > base.mid2id.tsv

        LC_ALL=C join -t $'\t' -1 1 -2 1 comp.mid.tsv base.mid2id.tsv \
            | awk -F'\t' -v tag="~{tag_value}" 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6,tag,$7,"SNV_indel",$8}' \
            > ~{prefix}.annotation.tsv
    >>>

    output {
        File annotation_tsv = "~{prefix}.annotation.tsv"
        File unmatched_vcf = "~{prefix}_truvari/fp.vcf.gz"
        File unmatched_vcf_idx = "~{prefix}_truvari/fp.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 4 * ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}



 
