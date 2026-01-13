version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl" as Helpers

workflow TruvariMatch {
    input {
        File vcf_eval
        File vcf_eval_index
        File vcf_truth
        File vcf_truth_index
        File ref_fa
        File ref_fai
        String prefix
        Int min_sv_length
        String utils_docker
        
        RuntimeAttr? runtime_attr_filter_eval_vcf
        RuntimeAttr? runtime_attr_filter_truth_vcf
        RuntimeAttr? runtime_attr_run_truvari_09
        RuntimeAttr? runtime_attr_annotate_matched_09
        RuntimeAttr? runtime_attr_run_truvari_07
        RuntimeAttr? runtime_attr_annotate_matched_07
        RuntimeAttr? runtime_attr_run_truvari_05
        RuntimeAttr? runtime_attr_annotate_matched_05
        RuntimeAttr? runtime_attr_concat_matched
    }

    call FilterEvalVcf {
        input:
            vcf_eval = vcf_eval,
            vcf_eval_index = vcf_eval_index,
            min_sv_length = min_sv_length,
            prefix = "~{prefix}.filtered_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_filter_eval_vcf
    }

    call FilterTruthVcf {
        input:
            vcf_truth = vcf_truth,
            vcf_truth_index = vcf_truth_index,
            prefix = "~{prefix}.filtered_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_filter_truth_vcf
    }

    call RunTruvari as RunTruvari_09 {
        input:
            vcf_eval = FilterEvalVcf.retained_vcf,
            vcf_eval_index = FilterEvalVcf.retained_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.retained_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.retained_vcf_index,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.9,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.9",
            prefix = "~{prefix}.0.9",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_09
    }

    call RunTruvari as RunTruvari_07 {
        input:
            vcf_eval = RunTruvari_09.unmatched_vcf,
            vcf_eval_index = RunTruvari_09.unmatched_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.retained_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.retained_vcf_index,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.7,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.7",
            prefix = "~{prefix}.0.7",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_07
    }

    call RunTruvari as RunTruvari_05 {
        input:
            vcf_eval = RunTruvari_07.unmatched_vcf,
            vcf_eval_index = RunTruvari_07.unmatched_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.retained_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.retained_vcf_index,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.5,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.5",
            prefix = "~{prefix}.0.5",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_05
    }

    call Helpers.ConcatTsvs as ConcatAnnotationTsvs {
        input:
            tsvs = [RunTruvari_09.annotation_tsv, RunTruvari_07.annotation_tsv, RunTruvari_05.annotation_tsv],
            prefix = "~{prefix}.truvari_combined",
            skip_sort = true,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_matched
    }

    output {
        File annotation_tsv = ConcatAnnotationTsvs.concatenated_tsv
        File unmatched_vcf = RunTruvari_05.unmatched_vcf
        File unmatched_vcf_index = RunTruvari_05.unmatched_vcf_index
        File dropped_vcf = FilterEvalVcf.dropped_vcf
        File dropped_vcf_index = FilterEvalVcf.dropped_vcf_index
    }
}

task FilterEvalVcf {
    input {
        File vcf_eval
        File vcf_eval_index
        Int min_sv_length
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        bcftools view -i "ABS(INFO/SVLEN)>=~{min_sv_length}" ~{vcf_eval} -Oz -o ~{prefix}.retained.vcf.gz
        tabix -p vcf -f ~{prefix}.retained.vcf.gz

        bcftools view -e 'ABS(INFO/SVLEN)>=~{min_sv_length}' ~{vcf_eval} -Oz -o ~{prefix}.dropped.vcf.gz
        tabix -p vcf -f ~{prefix}.dropped.vcf.gz
    >>>

    output {
        File retained_vcf = "~{prefix}.retained.vcf.gz"
        File retained_vcf_index = "~{prefix}.retained.vcf.gz.tbi"
        File dropped_vcf = "~{prefix}.dropped.vcf.gz"
        File dropped_vcf_index = "~{prefix}.dropped.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf_eval, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
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

task FilterTruthVcf {
    input {
        File vcf_truth
        File vcf_truth_index
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -e 'INFO/variant_type="snv"' ~{vcf_truth} \
            | bcftools view -i 'ABS(ILEN)>=5' -Oz -o ~{prefix}.vcf.gz
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File retained_vcf = "~{prefix}.vcf.gz"
        File retained_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf_truth, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
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

task RunTruvari {
    input {
        File vcf_eval
        File vcf_eval_index
        File vcf_truth_filtered
        File vcf_truth_filtered_index
        File ref_fa
        File ref_fai
        Float pctseq
        Int sizemin
        Int sizefilt
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
            -b ~{vcf_truth_filtered} \
            -c ~{vcf_eval} \
            -o "~{prefix}_truvari" \
            --reference ~{ref_fa} \
            --pctseq ~{pctseq} \
            --sizemin ~{sizemin} \
            --sizefilt ~{sizefilt}
        
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/MatchId\n' "~{prefix}_truvari/tp-comp.vcf.gz" \
            | awk 'BEGIN{FS=OFS="\t"} {split($6,a,","); print $1,$2,$3,$4,$5,a[1]}' \
            | LC_ALL=C sort -k6,6 > comp.mid.tsv

        bcftools query -f '%ID\t%INFO/MatchId\n' "~{prefix}_truvari/tp-base.vcf.gz" \
            | awk 'BEGIN{FS=OFS="\t"} {split($2,a,","); print a[1],$1}' \
            | LC_ALL=C sort -k1,1 > base.mid2id.tsv

        LC_ALL=C join -t $'\t' -1 6 -2 1 comp.mid.tsv base.mid2id.tsv \
            | awk -F'\t' -v tag="~{tag_value}" 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6,tag,$7}' \
            > ~{prefix}.annotation.tsv
    >>>

    output {
        File annotation_tsv = "~{prefix}.annotation.tsv"
        File unmatched_vcf = "~{prefix}_truvari/fp.vcf.gz"
        File unmatched_vcf_index = "~{prefix}_truvari/fp.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf_eval, "GB") + size(vcf_truth_filtered, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
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



 
