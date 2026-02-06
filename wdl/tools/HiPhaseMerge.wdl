version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "TRGTMerge.wdl"

workflow HiPhaseMerge {
    input {
        Array[File] phased_vcfs
        Array[File] phased_vcf_idxs
        Array[String] contigs
        String prefix

        String merge_args = "--merge id"

        File ref_fa
        File ref_fai

        String utils_docker
        String trgt_docker

        RuntimeAttr? runtime_attr_subset_trgt
        RuntimeAttr? runtime_attr_subset_integrated
        RuntimeAttr? runtime_attr_fix_al_header
        RuntimeAttr? runtime_attr_merge_trgt
        RuntimeAttr? runtime_attr_merge_integrated
        RuntimeAttr? runtime_attr_concat_trgt
        RuntimeAttr? runtime_attr_concat_integrated
    }

    scatter (i in range(length(phased_vcfs))) {
        call Helpers.SubsetVcfByArgs as SubsetTRGT {
            input:
                vcf = phased_vcfs[i],
                vcf_idx = phased_vcf_idxs[i],
                include_args = 'INFO/TRID != "."',
                prefix = "~{prefix}.~{i}.trgt",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_trgt
        }

        call FixALHeader {
            input:
                vcf = SubsetTRGT.subset_vcf,
                vcf_idx = SubsetTRGT.subset_vcf_idx,
                prefix = "~{prefix}.~{i}.trgt.fixed_al",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_fix_al_header
        }

        call Helpers.SubsetVcfByArgs as SubsetIntegrated {
            input:
                vcf = phased_vcfs[i],
                vcf_idx = phased_vcf_idxs[i],
                include_args = 'INFO/TRID = "."',
                prefix = "~{prefix}.~{i}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_integrated
        }
    }

    scatter (contig in contigs) {
        call TRGTMerge.TRGTMergeContig as MergeTRGTVcfs {
            input:
                vcfs = FixALHeader.fixed_vcf,
                vcf_idxs = FixALHeader.fixed_vcf_idx,
                prefix = "~{prefix}.~{contig}.trgt",
                contig = contig,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                docker = trgt_docker,
                runtime_attr_override = runtime_attr_merge_trgt
        }

        call Helpers.MergeVcfs as MergeIntegratedVcfs {
            input:
                vcfs = SubsetIntegrated.subset_vcf,
                vcf_idxs = SubsetIntegrated.subset_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.integrated",
                extra_args = merge_args,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_integrated
        }
    }

    call Helpers.ConcatVcfs as ConcatTRGTVcfs {
        input:
            vcfs = MergeTRGTVcfs.merged_vcf,
            vcf_idxs = MergeTRGTVcfs.merged_vcf_idx,
            prefix = "~{prefix}.trgt",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_trgt
    }

    call Helpers.ConcatVcfs as ConcatIntegratedVcfs {
        input:
            vcfs = MergeIntegratedVcfs.merged_vcf,
            vcf_idxs = MergeIntegratedVcfs.merged_vcf_idx,
            prefix = "~{prefix}.integrated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_integrated
    }

    output {
        File hiphase_merged_trgt_vcf = ConcatTRGTVcfs.concat_vcf
        File hiphase_merged_trgt_vcf_idx = ConcatTRGTVcfs.concat_vcf_idx
        File hiphase_merged_integrated_vcf = ConcatIntegratedVcfs.concat_vcf
        File hiphase_merged_integrated_vcf_idx = ConcatIntegratedVcfs.concat_vcf_idx
    }
}

task FixALHeader {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -h ~{vcf} \
            | grep -v '^##FORMAT=<ID=AL,' \
            > new_header.txt
        echo '##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">' >> new_header.txt

        bcftools reheader -h new_header.txt ~{vcf} \
            | bcftools view -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File fixed_vcf = "~{prefix}.vcf.gz"
        File fixed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
