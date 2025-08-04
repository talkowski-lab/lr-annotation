version 1.0

import "Structs.wdl"

workflow TruvariMatch {
    input {
        File vcf_eval
        File vcf_eval_index
        File vcf_truth
        File vcf_truth_index
        File ref_fasta
        File ref_fasta_fai
        String prefix
        String pipeline_docker
        String truvari_docker
        RuntimeAttr? runtime_attr_override
    }

    call FilterEvalVcf {
        input:
            vcf_eval = vcf_eval,
            vcf_eval_index = vcf_eval_index,
            prefix = "~{prefix}.filtered_eval",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    call FilterTruthVcf {
        input:
            vcf_truth = vcf_truth,
            vcf_truth_index = vcf_truth_index,
            prefix = "~{prefix}.filtered_truth",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    # Pass 1: pctseq = 0.9
    call RunTruvari as RunTruvari_09 {
        input:
            vcf_eval = FilterEvalVcf.retained_vcf,
            vcf_eval_index = FilterEvalVcf.retained_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.retained_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.retained_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            pctseq = 0.9,
            sizemin = 0,
            sizefilt = 0,
            prefix = "~{prefix}.0.9",
            truvari_docker = truvari_docker,
            runtime_attr_override = runtime_attr_override
    }

    call AnnotateVcf as AnnotateMatched_09 {
        input:
            vcf = RunTruvari_09.matched_vcf,
            vcf_index = RunTruvari_09.matched_vcf_index,
            tag_name = "gnomAD_V4_match",
            tag_value = "TRUVARI_0.9",
            prefix = "~{prefix}.0.9.annotated",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    # Pass 2: pctseq = 0.7
    call RunTruvari as RunTruvari_07 {
        input:
            vcf_eval = RunTruvari_09.unmatched_vcf,
            vcf_eval_index = RunTruvari_09.unmatched_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.retained_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.retained_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            pctseq = 0.7,
            sizemin = 0,
            sizefilt = 0,
            prefix = "~{prefix}.0.7",
            truvari_docker = truvari_docker,
            runtime_attr_override = runtime_attr_override
    }

    call AnnotateVcf as AnnotateMatched_07 {
        input:
            vcf = RunTruvari_07.matched_vcf,
            vcf_index = RunTruvari_07.matched_vcf_index,
            tag_name = "gnomAD_V4_match",
            tag_value = "TRUVARI_0.7",
            prefix = "~{prefix}.0.7.annotated",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    # Pass 3: pctseq = 0.5
    call RunTruvari as RunTruvari_05 {
        input:
            vcf_eval = RunTruvari_07.unmatched_vcf,
            vcf_eval_index = RunTruvari_07.unmatched_vcf_index,
            vcf_truth_filtered = FilterTruthVcf.retained_vcf,
            vcf_truth_filtered_index = FilterTruthVcf.retained_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            pctseq = 0.5,
            sizemin = 0,
            sizefilt = 0,
            prefix = "~{prefix}.0.5",
            truvari_docker = truvari_docker,
            runtime_attr_override = runtime_attr_override
    }

    call AnnotateVcf as AnnotateMatched_05 {
        input:
            vcf = RunTruvari_05.matched_vcf,
            vcf_index = RunTruvari_05.matched_vcf_index,
            tag_name = "gnomAD_V4_match",
            tag_value = "TRUVARI_0.5",
            prefix = "~{prefix}.0.5.annotated",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    call ConcatTruvariResults as ConcatMatched {
        input:
            vcfs = [AnnotateMatched_09.vcf_out, AnnotateMatched_07.vcf_out, AnnotateMatched_05.vcf_out],
            vcfs_idx = [AnnotateMatched_09.vcf_out_index, AnnotateMatched_07.vcf_out_index, AnnotateMatched_05.vcf_out_index],
            prefix = "~{prefix}.truvari_combined",
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File matched_vcf = ConcatMatched.concat_vcf
        File matched_vcf_index = ConcatMatched.concat_vcf_idx
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
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        bcftools view -i "ABS(INFO/SVLEN)>=10" ~{vcf_eval} -Oz -o ~{prefix}.retained.vcf.gz
        tabix -p vcf -f ~{prefix}.retained.vcf.gz

        bcftools view -e 'ABS(INFO/SVLEN)>=10' ~{vcf_eval} -Oz -o ~{prefix}.dropped.vcf.gz
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
        mem_gb: 8, 
        disk_gb: ceil(size(vcf_eval, "GB")) * 3 + 5,
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
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task FilterTruthVcf {
    input {
        File vcf_truth
        File vcf_truth_index
        String prefix
        String pipeline_docker
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
        mem_gb: 8, 
        disk_gb: ceil(size(vcf_truth, "GB")) * 3 + 5,
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
        set -euo pipefail

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
    >>>

    output {
        File matched_vcf = "~{prefix}_truvari/tp-comp.vcf.gz"
        File matched_vcf_index = "~{prefix}_truvari/tp-comp.vcf.gz.tbi"
        File unmatched_vcf = "~{prefix}_truvari/fp.vcf.gz"
        File unmatched_vcf_index = "~{prefix}_truvari/fp.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 8, 
        disk_gb: ceil(size(vcf_eval, "GB") + size(vcf_truth_filtered, "GB")) * 5 + 20,
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
        docker: truvari_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateVcf {
    input {
        File vcf
        File vcf_index
        String tag_name
        String tag_value
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        echo '##INFO=<ID=~{tag_name},Number=1,Type=String,Description="Matching status against gnomAD v4.">' > header.hdr
        echo '##INFO=<ID=gnomAD_V4_match_ID,Number=1,Type=String,Description="Matching variant ID from gnomAD v4.">' >> header.hdr

        # Extract MatchId values and create annotation file
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t~{tag_value}\t%INFO/MatchId\n' ~{vcf} | \
            awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, ($6 != "." && $6 != "") ? $6 : ""}' | \
            bgzip -c > annots.tab.gz
        tabix -s 1 -b 2 -e 2 annots.tab.gz

        bcftools annotate \
            -a annots.tab.gz \
            -h header.hdr \
            -c CHROM,POS,REF,ALT,~{tag_name},gnomAD_V4_match_ID \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}
        
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_out = "~{prefix}.vcf.gz"
        File vcf_out_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 4, 
        disk_gb: ceil(size(vcf, "GB")) * 2 + 5,
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
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatTruvariResults {
    input {
        Array[File] vcfs
        Array[File] vcfs_idx
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat -a -Oz -o ~{prefix}.vcf.gz ~{sep=' ' vcfs}
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File concat_vcf = "~{prefix}.vcf.gz"
        File concat_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

     RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 4, 
        disk_gb: ceil(size(vcfs, "GB")) * 2 + 5,
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
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
} 