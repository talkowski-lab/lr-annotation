version 1.0

import "Structs.wdl"

workflow AnnotateSTRs {
    input {
        Array[String] sample_ids
        File vcf
        File vcf_index
        File high_confidence_bed
        File ref_fasta
        File ref_fasta_fai
        String docker_image

        Boolean exclude_homopolymers = false
        Boolean only_pure_repeats = false
        Boolean keep_loci_that_have_overlapping_variants = false
        Int min_str_length = 9
        Int min_str_repeats = 3
        String output_suffix = "_str_variants"

        RuntimeAttr? runtime_attr_override
    }

    scatter (sample_id in sample_ids) {
        call PreprocessVcfSingleSample {
            input:
                sample_id = sample_id,
                vcf = vcf,
                vcf_index = vcf_index,
                high_confidence_bed = high_confidence_bed,
                docker_image = docker_image,
                runtime_attr_override = runtime_attr_override
        }

        call FilterVcfToSTR {
            input:
                sample_id = sample_id,
                vcf = PreprocessVcfSingleSample.high_confidence_vcf,
                vcf_index = PreprocessVcfSingleSample.high_confidence_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                docker_image = docker_image,
                exclude_homopolymers = exclude_homopolymers,
                only_pure_repeats = only_pure_repeats,
                keep_loci_that_have_overlapping_variants = keep_loci_that_have_overlapping_variants,
                min_str_length = min_str_length,
                min_str_repeats = min_str_repeats,
                output_suffix = output_suffix,
                runtime_attr_override = runtime_attr_override
        }

        call RestoreVariantID {
            input:
                vcf = FilterVcfToSTR.filtered_vcf,
                vcf_index = FilterVcfToSTR.filtered_vcf_index,
                sample_id = sample_id,
                output_suffix = output_suffix,
                docker_image = docker_image,
                runtime_attr_override = runtime_attr_override
        }
    }

    call MergeSTRVcfs {
        input:
            vcfs = RestoreVariantID.restored_id_vcf,
            vcf_indices = RestoreVariantID.restored_id_vcf_index,
            prefix = sample_ids[0],
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_override
    }

    call CombineWithOriginalVcf {
        input:
            original_vcf = vcf,
            original_vcf_index = vcf_index,
            str_vcf = MergeSTRVcfs.merged_vcf,
            str_vcf_index = MergeSTRVcfs.merged_vcf_index,
            prefix = sample_ids[0],
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File str_filtered_vcf = CombineWithOriginalVcf.final_vcf
        File str_filtered_vcf_index = CombineWithOriginalVcf.final_vcf_index
    }
}

task PreprocessVcfSingleSample {
    input {
        String sample_id
        File vcf
        File vcf_index
        File high_confidence_bed
        String docker_image

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: ceil(size(vcf, "GB") + size(high_confidence_bed, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -exuo pipefail

        if [ ! -s "~{high_confidence_bed}" ]; then
            echo "High confidence BED file is empty. Exiting."
            exit 1
        fi

        bcftools view -s ~{sample_id} -c 1 -Oz -o ~{sample_id}.subset.vcf.gz ~{vcf}
        tabix ~{sample_id}.subset.vcf.gz

        bedtools intersect -header -f 1 -wa -u \
            -a ~{sample_id}.subset.vcf.gz \
            -b ~{high_confidence_bed} \
            | bgzip > ~{sample_id}.high_confidence_regions.vcf.gz
        tabix -f ~{sample_id}.high_confidence_regions.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File high_confidence_vcf = "~{sample_id}.high_confidence_regions.vcf.gz"
        File high_confidence_vcf_index = "~{sample_id}.high_confidence_regions.vcf.gz.tbi"
    }
}

task FilterVcfToSTR {
    input {
        String sample_id
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fasta_fai
        String docker_image

        Boolean exclude_homopolymers
        Boolean only_pure_repeats
        Boolean keep_loci_that_have_overlapping_variants
        Int min_str_length
        Int min_str_repeats
        String output_suffix

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: ceil(size(vcf, "GB") + size(ref_fasta, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String min_repeat_unit_length = if exclude_homopolymers then "2" else "1"
    String allow_interruptions = if only_pure_repeats then "no" else "only-if-pure-repeats-not-found"
    String keep_loci_arg = if keep_loci_that_have_overlapping_variants then "--keep-loci-that-have-overlapping-variants" else ""

    command <<<
        set -exuo pipefail

        if [[ $(zgrep -vc "^#" ~{vcf}) -eq 0 ]]; then
            touch ~{sample_id}~{output_suffix}.vcf.gz
            touch ~{sample_id}~{output_suffix}.vcf.gz.tbi
            touch ~{sample_id}~{output_suffix}.variants.tsv.gz
            touch ~{sample_id}~{output_suffix}.alleles.tsv.gz
            touch ~{sample_id}~{output_suffix}.variants.bed.gz
            touch ~{sample_id}~{output_suffix}.variants.bed.gz.tbi
            touch ~{sample_id}~{output_suffix}.filter_vcf.log
            exit 0
        fi

        python3 -u -m str_analysis.filter_vcf_to_STR_variants \
            -R ~{ref_fasta} \
            --allow-interruptions ~{allow_interruptions} \
            --write-bed-file \
            --min-str-length ~{min_str_length} \
            --min-str-repeats ~{min_str_repeats} \
            --min-repeat-unit-length ~{min_repeat_unit_length} \
            ~{keep_loci_arg} \
            --output-prefix ~{sample_id}~{output_suffix} \
            --verbose \
            ~{vcf} |& tee ~{sample_id}~{output_suffix}.filter_vcf.log
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File filtered_vcf = "~{sample_id}~{output_suffix}.vcf.gz"
        File filtered_vcf_index = "~{sample_id}~{output_suffix}.vcf.gz.tbi"
        File variants_tsv = "~{sample_id}~{output_suffix}.variants.tsv.gz"
        File alleles_tsv = "~{sample_id}~{output_suffix}.alleles.tsv.gz"
        File variants_bed = "~{sample_id}~{output_suffix}.variants.bed.gz"
        File variants_bed_index = "~{sample_id}~{output_suffix}.variants.bed.gz.tbi"
        File filter_log = "~{sample_id}~{output_suffix}.filter_vcf.log"
    }
}

task RestoreVariantID {
    input {
        String sample_id
        File vcf
        File vcf_index
        String output_suffix
        String docker_image

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -exuo pipefail

        if [[ $(zgrep -vc "^#" ~{vcf}) -eq 0 ]]; then
            touch ~{sample_id}~{output_suffix}.restored_id.vcf.gz
            touch ~{sample_id}~{output_suffix}.restored_id.vcf.gz.tbi
        else
            bcftools annotate --set-id '%INFO/ID' ~{vcf} -Oz -o ~{sample_id}~{output_suffix}.restored_id.vcf.gz
            tabix ~{sample_id}~{output_suffix}.restored_id.vcf.gz
        fi
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File restored_id_vcf = "~{sample_id}~{output_suffix}.restored_id.vcf.gz"
        File restored_id_vcf_index = "~{sample_id}~{output_suffix}.restored_id.vcf.gz.tbi"
    }
}

task MergeSTRVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        String prefix
        String docker_image

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: ceil(size(vcfs, "GB") * 1.5) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -exuo pipefail

        ls ~{sep=' ' vcfs} > vcf_list.txt
        bcftools merge --file-list vcf_list.txt -Oz -o ~{prefix}.merged.vcf.gz
        tabix ~{prefix}.merged.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File merged_vcf = "~{prefix}.merged.vcf.gz"
        File merged_vcf_index = "~{prefix}.merged.vcf.gz.tbi"
    }
}

task CombineWithOriginalVcf {
    input {
        File original_vcf
        File original_vcf_index
        File str_vcf
        File str_vcf_index
        String prefix
        String docker_image

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: ceil(size(original_vcf, "GB") + size(str_vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -exuo pipefail

        bcftools query -f '%ID\n' ~{str_vcf} > str_ids.txt

        bcftools view ~{original_vcf} -e 'ID=@str_ids.txt' -Oz -o non_str_variants.vcf.gz

        bcftools concat -a -D -Oz -o ~{prefix}.final.vcf.gz non_str_variants.vcf.gz ~{str_vcf}
        tabix ~{prefix}.final.vcf.gz
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }

    output {
        File final_vcf = "~{prefix}.final.vcf.gz"
        File final_vcf_index = "~{prefix}.final.vcf.gz.tbi"
    }
} 