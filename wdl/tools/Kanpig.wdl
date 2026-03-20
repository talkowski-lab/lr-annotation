version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow Kanpig {
    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[String] sample_ids
        Array[File] bams
        Array[File] bais
        Array[String] sexes
        String prefix

        File ref_fa
        File ref_fai
        File ploidy_bed_male
        File ploidy_bed_female

        String merge_args = "--merge id"
        String kanpig_params = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"

        String kanpig_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_to_sample
        RuntimeAttr? runtime_attr_run_kanpig
        RuntimeAttr? runtime_attr_merge_vcfs
    }

    scatter (i in range(length(sample_ids))) {
        File ploidy_bed = if sexes[i] == "M" then ploidy_bed_male else ploidy_bed_female

        call Helpers.SubsetVcfToSamples {
            input:
                vcf = cohort_vcf,
                vcf_idx = cohort_vcf_idx,
                samples = [sample_ids[i]],
                filter_to_sample = false,
                prefix = sample_ids[i] + ".subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_to_sample
        }

        call RunKanpig {
            input:
                input_vcf = SubsetVcfToSamples.subset_vcf,
                input_vcf_idx = SubsetVcfToSamples.subset_vcf_idx,
                bam = bams[i],
                bai = bais[i],
                sample_id = sample_ids[i],
                ploidy_bed = ploidy_bed,
                kanpig_params = kanpig_params,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = sample_ids[i],
                docker = kanpig_docker,
                runtime_attr_override = runtime_attr_run_kanpig
        }
    }

    call Helpers.MergeVcfs {
        input:
            vcfs = RunKanpig.regenotyped_vcf,
            vcf_idxs = RunKanpig.regenotyped_vcf_idx,
            prefix = prefix,
            extra_args = merge_args,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    output {
        File kanpig_regenotyped_vcf = MergeVcfs.merged_vcf
        File kanpig_regenotyped_vcf_idx = MergeVcfs.merged_vcf_idx
    }
}

task RunKanpig {
    input {
        File input_vcf
        File input_vcf_idx
        File bam
        File bai
        String sample_id
        File ploidy_bed
        String kanpig_params
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        nproc=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        kanpig gt \
            --threads "${nproc}" \
            --sizemin 10 \
            ~{kanpig_params} \
            --reference ~{ref_fa} \
            --ploidy-bed ~{ploidy_bed} \
            --input ~{input_vcf} \
            --reads ~{bam} \
            --sample ~{sample_id} \
            --out ~{prefix}.kanpig.vcf

        bcftools sort \
            -Oz -o ~{prefix}.kanpig.vcf.gz \
            ~{prefix}.kanpig.vcf
        
        rm ~{prefix}.kanpig.vcf

        tabix -p vcf ~{prefix}.kanpig.vcf.gz
    >>>

    output {
        File regenotyped_vcf = "~{prefix}.kanpig.vcf.gz"
        File regenotyped_vcf_idx = "~{prefix}.kanpig.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 16,
        disk_gb: ceil(size(bam, "GB")) * 2 + ceil(size(input_vcf, "GB")) + 20,
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
