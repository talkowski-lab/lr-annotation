version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CombineVcfs {
    input {
        File a_vcf
        File a_vcf_idx
        File b_vcf
        File b_vcf_idx
        String contig
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_a
        RuntimeAttr? runtime_attr_subset_b
        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
    }

    call Helpers.SubsetVcfToContig as SubsetA {
        input:
            vcf = a_vcf,
            vcf_idx = a_vcf_idx,
            contig = contig,
            prefix = "~{prefix}.a.~{contig}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_a
    }

    call Helpers.SubsetVcfToContig as SubsetB {
        input:
            vcf = b_vcf,
            vcf_idx = b_vcf_idx,
            contig = contig,
            prefix = "~{prefix}.b.~{contig}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_b
    }

    call CheckSampleMatch {
        input:
            vcf_a = SubsetA.subset_vcf,
            vcf_a_idx = SubsetA.subset_vcf_idx,
            vcf_b = SubsetB.subset_vcf,
            vcf_b_idx = SubsetB.subset_vcf_idx,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    call Helpers.MergeVcfs {
        input:
            vcfs = [SubsetA.subset_vcf, SubsetB.subset_vcf],
            vcf_idxs = [SubsetA.subset_vcf_idx, SubsetB.subset_vcf_idx],
            prefix = "~{prefix}.~{contig}.merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = [SubsetA.subset_vcf, SubsetB.subset_vcf],
            vcf_idxs = [SubsetA.subset_vcf_idx, SubsetB.subset_vcf_idx],
            prefix = "~{prefix}.~{contig}.concatenated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File merged_vcf = MergeVcfs.merged_vcf
        File merged_vcf_idx = MergeVcfs.merged_vcf_idx
        File concatenated_vcf = ConcatVcfs.concat_vcf
        File concatenated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task CheckSampleMatch {
    input {
        File vcf_a
        File vcf_a_idx
        File vcf_b
        File vcf_b_idx
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf_a} | sort > samples_a.txt
        bcftools query -l ~{vcf_b} | sort > samples_b.txt

        if ! diff -q samples_a.txt samples_b.txt > /dev/null; then
            echo "ERROR: Sample mismatch between VCFs"
            echo "Samples only in VCF A:"
            comm -23 samples_a.txt samples_b.txt
            echo "Samples only in VCF B:"
            comm -13 samples_a.txt samples_b.txt
            exit 1
        fi

        echo "SUCCESS: Samples match between VCFs"
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf_a, vcf_b], "GB")) + 10,
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
