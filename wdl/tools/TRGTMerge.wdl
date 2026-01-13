version 1.0

import "../utils/Structs.wdl"
import "../utils/MergeSplitVCF.wdl" as MergeSplitVCF

workflow TRGTMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs

        String prefix
        Array[String] contigs

        String trgt_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_override_merge_vcfs
        RuntimeAttr? runtime_attr_override_concat_vcfs
    }

    scatter (contig in contigs) {
        call MergeVCFsByContig {
            input:
                vcfs = vcfs,
                vcf_idxs = vcf_idxs,
                prefix = prefix,
                contig = contig,
                docker = trgt_docker,
                runtime_attr_override = runtime_attr_override_merge_vcfs
        }
    }

    call MergeSplitVCF.CombineVCFs {
        input:
            vcf_files = MergeVCFsByContig.merged_vcf,
            vcf_indices = MergeVCFsByContig.merged_vcf_idx,
            naive = true,
            allow_overlaps = false,
            sv_base_mini_docker = sv_base_mini_docker,
            cohort_prefix = prefix + ".final",
            sort_after_merge = false,
            runtime_attr_override = runtime_attr_override_concat_vcfs
    }

    output {
        File trgt_merged_vcf = CombineVCFs.combined_vcf
        File trgt_merged_vcf_idx = CombineVCFs.combined_vcf_idx
    }
}

task MergeVCFsByContig {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String contig
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcfs, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    
    Int threads = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euo pipefail

        # Create VCF list file
        echo "~{sep='\n' vcfs}" > vcf_list.txt

        trgt merge \
            --vcf-list vcf_list.txt \
            --contig ~{contig} \
            --threads ~{threads} \
            --write-index \
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
