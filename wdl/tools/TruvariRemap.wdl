version 1.0

import "../utils/Structs.wdl"

workflow CallInsRemap {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        File ref_fa
        File ref_fai
        Array[File] ref_bwa_indices

        Int? minlength
        Int? maxlength
        Int? mm2_threshold
        Float? cov_threshold

        RuntimeAttr? runtime_attr_ins_remap
        RuntimeAttr? runtime_attr_merge
    }

    scatter (contig in contigs) {
        call InsRemap {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                ref_bwa_indices = ref_bwa_indices,
                minlength = minlength,
                maxlength = maxlength,
                mm2_threshold = mm2_threshold,
                cov_threshold = cov_threshold,
                prefix = "~{prefix}.~{contig}.remapped",
                runtime_attr_override = runtime_attr_ins_remap
        }
    }

    call MergePerChrCalls {
        input:
            vcfs = InsRemap.remapped_vcf,
            prefix = "~{prefix}.remapped",
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File remapped_vcf = MergePerChrCalls.vcf
        File remapped_vcf_idx = MergePerChrCalls.tbi
    }
}

task InsRemap {
    input {
        File vcf
        File vcf_idx
        String contig
        File ref_fa
        File ref_fai
        Array[File] ref_bwa_indices
        Int minlength = 50
        Int maxlength = 10000000
        Int mm2_threshold = 5000
        Float cov_threshold = 0.8
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        bcftools view \
            -Oz \
            ~{vcf} \
            ~{contig} \
            > ~{contig}.vcf.gz
        
        tabix ~{contig}.vcf.gz

        mkdir ref_files
        mv ~{ref_fa} ref_files/
        mv ~{ref_fai} ref_files/
        for x in ~{sep=' ' ref_bwa_indices}; do
            mv $x ref_files/
        done

        ref_fa_basename=$(basename ~{ref_fa})

        truvari anno remap \
            -r ref_files/$ref_fa_basename \
            -o ~{prefix}.vcf.gz \
            --min-length ~{minlength} \
            --max-length ~{maxlength} \
            --mm2-threshold ~{mm2_threshold} \
            --threads ${N_THREADS} \
            --cov-threshold ~{cov_threshold} \
            ~{contig}.vcf.gz
    >>>
    
    output {
        File remapped_vcf = "~{prefix}.vcf.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 16,
        disk_gb: ceil(size(vcf, "GB") + size(ref_fa, "GB")) * 3 + 50,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-remap:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergePerChrCalls {
    input {
        Array[File] vcfs
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        VCF_WITH_HEADER=~{vcfs[0]}

        bcftools view -h $VCF_WITH_HEADER > tmp.vcf

        for vcf in ~{sep=" " vcfs}; do
            bcftools view -H $vcf >> tmp.vcf
        done

        bcftools sort -Oz tmp.vcf > ~{prefix}.vcf.gz
        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File unsorted_vcf = "tmp.vcf"
        File vcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 24,
        disk_gb: 5 * ceil(size(vcfs, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}