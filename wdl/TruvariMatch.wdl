version 1.0

import "Structs.wdl"
import "Helpers.wdl" as Helpers

workflow TruvariMatch {
    input {
        File vcf_eval_unmatched
        File vcf_eval_unmatched_index
        File vcf_truth
        File vcf_truth_index
        File ref_fasta
        File ref_fasta_fai
        String prefix
        String pipeline_docker
        String truvari_docker
        RuntimeAttr? runtime_attr_override
    }

    call SplitEvalVcfBySize {
        input:
            vcf = vcf_eval_unmatched,
            vcf_index = vcf_eval_unmatched_index,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    call FilterTruthVcf {
        input:
            vcf = vcf_truth,
            vcf_index = vcf_truth_index,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    # Pass 1: pctseq=0.9
    call RunTruvari as RunTruvari_09 {
        input:
            vcf_eval = SplitEvalVcfBySize.passed_vcf,
            vcf_eval_index = SplitEvalVcfBySize.passed_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.filtered_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.filtered_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            pctseq = 0.9,
            sizemin = 0,
            sizefilt = 0,
            prefix = "~{prefix}.0.9",
            truvari_docker = truvari_docker,
            runtime_attr_override = runtime_attr_override
    }
    call AnnotateVcf as Annotate_09 {
        input:
            vcf_in = RunTruvari_09.matched_vcf,
            vcf_in_index = RunTruvari_09.matched_vcf_index,
            prefix = "~{prefix}.0.9.annotated",
            tag_value = "TRUVARI_0.9",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    # Pass 2: pctseq=0.7
    call RunTruvari as RunTruvari_07 {
        input:
            vcf_eval = RunTruvari_09.unmatched_vcf,
            vcf_eval_index = RunTruvari_09.unmatched_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.filtered_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.filtered_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            pctseq = 0.7,
            sizemin = 0,
            sizefilt = 0,
            prefix = "~{prefix}.0.7",
            truvari_docker = truvari_docker,
            runtime_attr_override = runtime_attr_override
    }
    call AnnotateVcf as Annotate_07 {
        input:
            vcf_in = RunTruvari_07.matched_vcf,
            vcf_in_index = RunTruvari_07.matched_vcf_index,
            prefix = "~{prefix}.0.7.annotated",
            tag_value = "TRUVARI_0.7",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    # Pass 3: pctseq=0.5
    call RunTruvari as RunTruvari_05 {
        input:
            vcf_eval = RunTruvari_07.unmatched_vcf,
            vcf_eval_index = RunTruvari_07.unmatched_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.filtered_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.filtered_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            pctseq = 0.5,
            sizemin = 0,
            sizefilt = 0,
            prefix = "~{prefix}.0.5",
            truvari_docker = truvari_docker,
            runtime_attr_override = runtime_attr_override
    }
    call AnnotateVcf as Annotate_05 {
        input:
            vcf_in = RunTruvari_05.matched_vcf,
            vcf_in_index = RunTruvari_05.matched_vcf_index,
            prefix = "~{prefix}.0.5.annotated",
            tag_value = "TRUVARI_0.5",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    call ConcatTruvariResults as ConcatMatched {
        input:
            vcfs = [Annotate_09.vcf_out, Annotate_07.vcf_out, Annotate_05.vcf_out],
            vcfs_idx = [Annotate_09.vcf_out_index, Annotate_07.vcf_out_index, Annotate_05.vcf_out_index],
            final_unmatched_vcf = RunTruvari_05.unmatched_vcf,
            prefix = "~{prefix}.truvari_matched",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File matched_vcf = ConcatMatched.concat_vcf
        File matched_vcf_index = ConcatMatched.concat_vcf_idx
        File unmatched_vcf = RunTruvari_05.unmatched_vcf
        File unmatched_vcf_index = RunTruvari_05.unmatched_vcf_index
        File unmatched_too_small_vcf = SplitEvalVcfBySize.failed_vcf
        File unmatched_too_small_vcf_index = SplitEvalVcfBySize.failed_vcf_index
    }
}

task FilterTruthVcf {
    input {
        File vcf
        File vcf_index
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -euxo pipefail
        bcftools view -e 'INFO/variant_type="snv"' ~{vcf} \
            | bcftools view -i 'ABS(ILEN)>=5' -Oz -o ~{prefix}.filtered.vcf.gz
        
        tabix -p vcf -f ~{prefix}.filtered.vcf.gz
    >>>
    output {
        File filtered_vcf = "~{prefix}.filtered.vcf.gz"
        File filtered_vcf_index = "~{prefix}.filtered.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1, mem_gb: 4, disk_gb: ceil(size(vcf, "GB")) * 2 + 10,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SplitEvalVcfBySize {
    input {
        File vcf
        File vcf_index
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -euxo pipefail
        
        bcftools view -i 'ABS(INFO/SVLEN)>=10' ~{vcf} -Oz -o ~{prefix}.passed_size_filt.vcf.gz
        tabix -p vcf -f ~{prefix}.passed_size_filt.vcf.gz

        bcftools view -e 'ABS(INFO/SVLEN)>=10' ~{vcf} -Oz -o ~{prefix}.failed_size_filt.vcf.gz
        tabix -p vcf -f ~{prefix}.failed_size_filt.vcf.gz
    >>>
    output {
        File passed_vcf = "~{prefix}.passed_size_filt.vcf.gz"
        File passed_vcf_index = "~{prefix}.passed_size_filt.vcf.gz.tbi"
        File failed_vcf = "~{prefix}.failed_size_filt.vcf.gz"
        File failed_vcf_index = "~{prefix}.failed_size_filt.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, mem_gb: 4, disk_gb: ceil(size(vcf, "GB")) * 3 + 10,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
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
        File ref_fasta
        File ref_fasta_fai
        Float pctseq
        Int sizemin
        Int sizefilt
        String prefix
        String truvari_docker
        RuntimeAttr? runtime_attr_override
    }
    command <<<
        set -euxo pipefail

        if [ -d "~{prefix}_truvari" ]; then
            rm -r "~{prefix}_truvari"
        fi
        
        truvari bench \
            -b ~{vcf_truth_filtered} \
            -c ~{vcf_eval} \
            -o "~{prefix}_truvari" \
            --reference ~{ref_fasta} \
            --pctseq ~{pctseq} \
            --sizemin ~{sizemin} \
            --sizefilt ~{sizefilt}
        
        tabix -p vcf -f ~{prefix}_truvari/tp-comp.vcf.gz
        tabix -p vcf -f ~{prefix}_truvari/fp.vcf.gz
    >>>
    output {
        File matched_vcf = "~{prefix}_truvari/tp-comp.vcf.gz"
        File matched_vcf_index = "~{prefix}_truvari/tp-comp.vcf.gz.tbi"
        File unmatched_vcf = "~{prefix}_truvari/fp.vcf.gz"
        File unmatched_vcf_index = "~{prefix}_truvari/fp.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1, mem_gb: 8, disk_gb: ceil(size(vcf_eval, "GB") + size(vcf_truth_filtered, "GB")) * 5 + 20,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: truvari_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateVcf {
    input {
        File vcf_in
        File vcf_in_index
        String prefix
        String tag_value
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        
        bcftools annotate -a ~{vcf_in} -c INFO/gnomAD_V4_match -m "+'~{tag_value}'" ~{vcf_in} -Oz -o ~{prefix}.vcf.gz
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_out = "~{prefix}.vcf.gz"
        File vcf_out_index = "~{prefix}.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1, mem_gb: 4, disk_gb: ceil(size(vcf_in, "GB")) * 2 + 10,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatTruvariResults {
    input {
        Array[File] vcfs
        Array[File] vcfs_idx
        File final_unmatched_vcf # Used to get the header for reheadering
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        bcftools view -h ~{final_unmatched_vcf} > header.txt
        
        for vcf in ~{sep=' ' vcfs}; do
            bcftools reheader -h header.txt ${vcf} | bcftools view -Oz >> reheadered_vcfs.list
        done

        if [ -s reheadered_vcfs.list ]; then
            bcftools concat -a -f reheadered_vcfs.list -Oz -o ~{prefix}.combined.vcf.gz
            tabix -p vcf -f ~{prefix}.combined.vcf.gz
        else
            touch ~{prefix}.combined.vcf.gz
            touch ~{prefix}.combined.vcf.gz.tbi
        fi
    >>>
    
    output {
        File concat_vcf = "~{prefix}.combined.vcf.gz"
        File concat_vcf_idx = "~{prefix}.combined.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, mem_gb: 8, disk_gb: ceil(size(vcfs, "GB")) * 2 + 10,
        boot_disk_gb: 10, preemptible_tries: 2, max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
} 