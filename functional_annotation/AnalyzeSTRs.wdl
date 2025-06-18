version 1.0

workflow AnalyzeSTRs {
    input {
        String sample_id
        File vcf
        File vcf_index
        File high_confidence_bed
        File ref_fasta
        File ref_fasta_fai

        String? docker_image
        Boolean exclude_homopolymers = false
        Boolean only_pure_repeats = false
        Boolean keep_loci_that_have_overlapping_variants = false
        Int min_str_length = 9
        Int min_str_repeats = 3
        String output_suffix = "_str_variants"

        Int preemptible_tries = 3
        Int machine_mem_gb = 16
        Int machine_cpu = 2
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
            preemptible_tries = preemptible_tries,
            machine_mem_gb = machine_mem_gb,
            machine_cpu = machine_cpu
    }

    output {
        File filtered_vcf = FilterVcfToSTR.filtered_vcf
        File filtered_vcf_index = FilterVcfToSTR.filtered_vcf_index
        File variants_tsv = FilterVcfToSTR.variants_tsv
        File alleles_tsv = FilterVcfToSTR.alleles_tsv
        File variants_bed = FilterVcfToSTR.variants_bed
        File variants_bed_index = FilterVcfToSTR.variants_bed_index
        File filter_log = FilterVcfToSTR.filter_log
        File high_confidence_vcf = FilterVcfToSTR.high_confidence_vcf
        File high_confidence_vcf_index = FilterVcfToSTR.high_confidence_vcf_index
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

        String? docker_image
        Boolean exclude_homopolymers
        Boolean only_pure_repeats
        Boolean keep_loci_that_have_overlapping_variants
        Int min_str_length
        Int min_str_repeats
        String output_suffix

        Int preemptible_tries
        Int machine_mem_gb
        Int machine_cpu
    }

    Int machine_disk_gb = ceil(size(vcf, "GB") + size(ref_fasta, "GB") + size(high_confidence_bed, "GB")) + 20

    String min_repeat_unit_length = if exclude_homopolymers then "2" else "1"
    String allow_interruptions = if only_pure_repeats then "no" else "only-if-pure-repeats-not-found"
    String keep_loci_arg = if keep_loci_that_have_overlapping_variants then "--keep-loci-that-have-overlapping-variants" else ""
    String docker = select_first([docker_image, "broadinstitute/gatk:latest"])

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
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu: machine_cpu
        disks: "local-disk " + machine_disk_gb + " HDD"
        preemptible: preemptible_tries
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