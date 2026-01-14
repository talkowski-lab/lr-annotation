version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow TRGTMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs

        String prefix
        Array[String] contigs

        File ref_fa
        File ref_fai

        String trgt_docker
        String utils_docker

        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call MergeVCFsByContig {
            input:
                vcfs = vcfs,
                vcf_idxs = vcf_idxs,
                prefix = prefix,
                contig = contig,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                docker = trgt_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeVCFsByContig.merged_vcf,
            vcfs_idx = MergeVCFsByContig.merged_vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File trgt_merged_vcf = ConcatVcfs.concat_vcf
        File trgt_merged_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task MergeVCFsByContig {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String contig
        File ref_fa
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    
    command <<<
        set -euo pipefail

        echo "~{sep='\n' vcfs}" > vcf_list.txt

        trgt merge \
            --vcf-list vcf_list.txt \
            --contig ~{contig} \
            --genome ~{ref_fa} \
            --output-type z \
            --output ~{prefix}.~{contig}.merged.trgt.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.~{contig}.merged.trgt.vcf.gz"
        File merged_vcf_idx = "~{prefix}.~{contig}.merged.trgt.vcf.gz.tbi"
    }

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
