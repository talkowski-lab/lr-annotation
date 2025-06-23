version 1.0

import "Structs.wdl"

workflow AnnotateSTRs {
    input {
        String sample_id
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

    call FilterVcfToSTR {
        input:
            sample_id = sample_id,
            vcf = vcf,
            vcf_index = vcf_index,
            high_confidence_bed = high_confidence_bed,
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

    output {
        File str_filtered_vcf = FilterVcfToSTR.filtered_vcf
        File str_filtered_vcf_index = FilterVcfToSTR.filtered_vcf_index
        File str_variants_tsv = FilterVcfToSTR.variants_tsv
        File str_alleles_tsv = FilterVcfToSTR.alleles_tsv
        File str_variants_bed = FilterVcfToSTR.variants_bed
        File str_variants_bed_index = FilterVcfToSTR.variants_bed_index
        File str_filter_log = FilterVcfToSTR.filter_log
        File str_high_confidence_vcf = FilterVcfToSTR.high_confidence_vcf
        File str_high_confidence_vcf_index = FilterVcfToSTR.high_confidence_vcf_index
    }
}

task FilterVcfToSTR {
    input {
        String sample_id
        File vcf
        File vcf_index
        File high_confidence_bed
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
        disk_gb: ceil(size(vcf, "GB") + size(ref_fasta, "GB") + size(high_confidence_bed, "GB")) + 20,
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

        if [ ! -s "~{high_confidence_bed}" ]; then
            echo "High confidence BED file is empty. Exiting."
            exit 1
        fi

        bedtools intersect -header -f 1 -wa -u \
            -a ~{vcf} \
            -b ~{high_confidence_bed} \
            | bgzip > ~{sample_id}.high_confidence_regions.vcf.gz

        tabix -f ~{sample_id}.high_confidence_regions.vcf.gz

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
            ~{sample_id}.high_confidence_regions.vcf.gz |& tee ~{sample_id}~{output_suffix}.filter_vcf.log
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
        File high_confidence_vcf = "~{sample_id}.high_confidence_regions.vcf.gz"
        File high_confidence_vcf_index = "~{sample_id}.high_confidence_regions.vcf.gz.tbi"
    }
} 