version 1.0

import "utils/Structs.wdl"
import "utils/Helpers.wdl" as Helpers

workflow AnnotateSTRs {
    input {
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fasta_fai
        String pipeline_docker

        Boolean exclude_homopolymers = false
        Boolean only_pure_repeats = false
        Boolean keep_loci_that_have_overlapping_variants = false
        Int min_str_length = 9
        Int min_str_repeats = 3
        String output_suffix = "_str_variants"
        Int? variants_per_shard

        RuntimeAttr? runtime_attr_preprocess
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_postprocess
        RuntimeAttr? runtime_attr_annotate
    }

    String prefix = basename(vcf, ".vcf.gz")

    Int variants_per_shard_eff = select_first([variants_per_shard, 1000000000])

    call Helpers.SplitVcfIntoShards {
        input:
            input_vcf = vcf,
            input_vcf_index = vcf_index,
            variants_per_shard = variants_per_shard_eff,
            output_prefix = prefix,
            docker = pipeline_docker
    }

    scatter (shard in zip(SplitVcfIntoShards.split_vcfs, SplitVcfIntoShards.split_vcf_indexes)) {
        String shard_prefix = basename(shard.left, ".vcf.gz")
        
        call PreprocessVcf {
            input:
                vcf = shard.left,
                vcf_index = shard.right,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_preprocess
        }

        call FilterVcfToSTR {
            input:
                prefix = shard_prefix,
                vcf = PreprocessVcf.dummy_vcf,
                vcf_index = PreprocessVcf.dummy_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                pipeline_docker = pipeline_docker,
                exclude_homopolymers = exclude_homopolymers,
                only_pure_repeats = only_pure_repeats,
                keep_loci_that_have_overlapping_variants = keep_loci_that_have_overlapping_variants,
                min_str_length = min_str_length,
                min_str_repeats = min_str_repeats,
                output_suffix = output_suffix,
                runtime_attr_override = runtime_attr_filter
        }

        call PostprocessVcf {
            input:
                vcf = FilterVcfToSTR.filtered_vcf,
                vcf_index = FilterVcfToSTR.filtered_vcf_index,
                prefix = shard_prefix,
                output_suffix = output_suffix,
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_postprocess
        }
    }

    call Helpers.ConcatVcfs as ConcatSTRs {
        input:
            vcfs = PostprocessVcf.restored_id_vcf,
            vcfs_idx = PostprocessVcf.restored_id_vcf_index,
            outfile_prefix = prefix,
            docker = pipeline_docker
    }

    call AnnotateOriginalVcf {
        input:
            original_vcf = vcf,
            original_vcf_index = vcf_index,
            str_vcf = ConcatSTRs.concat_vcf,
            str_vcf_index = ConcatSTRs.concat_vcf_idx,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_annotate
    }

    output {
        File str_annotated_vcf = AnnotateOriginalVcf.final_vcf
        File str_annotated_vcf_index = AnnotateOriginalVcf.final_vcf_index
    }
}

task PreprocessVcf {
    input {
        File vcf
        File vcf_index
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -exuo pipefail

        SAMPLE=$(bcftools query -l ~{vcf} | head -n 1)

        bcftools view -s "${SAMPLE}" -O v --no-version ~{vcf} | \
            bcftools reheader --samples <(echo "input_sample") | \
            bcftools +setGT -- -t a -n c:'0/1' | \
            bgzip > "dummy_sample.vcf.gz"

        tabix "dummy_sample.vcf.gz"
    >>>

    output {
        File dummy_vcf = "dummy_sample.vcf.gz"
        File dummy_vcf_index = "dummy_sample.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
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

task FilterVcfToSTR {
    input {
        String prefix
        File vcf
        File vcf_index
        File ref_fasta
        File ref_fasta_fai
        String pipeline_docker

        Boolean exclude_homopolymers
        Boolean only_pure_repeats
        Boolean keep_loci_that_have_overlapping_variants
        Int min_str_length
        Int min_str_repeats
        String output_suffix

        RuntimeAttr? runtime_attr_override
    }

    String min_repeat_unit_length = if exclude_homopolymers then "2" else "1"
    String allow_interruptions = if only_pure_repeats then "no" else "only-if-pure-repeats-not-found"
    String keep_loci_arg = if keep_loci_that_have_overlapping_variants then "--keep-loci-that-have-overlapping-variants" else ""

    command <<<
        set -exuo pipefail

        FINAL_PREFIX="~{prefix}~{output_suffix}"

        python3 -u -m str_analysis.filter_vcf_to_STR_variants \
            -R ~{ref_fasta} \
            --allow-interruptions ~{allow_interruptions} \
            --write-bed-file \
            --min-str-length ~{min_str_length} \
            --min-str-repeats ~{min_str_repeats} \
            --min-repeat-unit-length ~{min_repeat_unit_length} \
            ~{keep_loci_arg} \
            --output-prefix ${FINAL_PREFIX} \
            --verbose \
            ~{vcf} |& tee ${FINAL_PREFIX}.filter_vcf.log

        cat <<EOT > new_header_fields.txt
##INFO=<ID=LocusId,Number=1,Type=String,Description="Locus ID">
##INFO=<ID=Locus,Number=1,Type=String,Description="Locus coordinates">
##INFO=<ID=Motif,Number=1,Type=String,Description="Repeat motif">
##INFO=<ID=NumRepeatsShortAllele,Number=1,Type=Integer,Description="Number of repeats in the shorter allele">
##INFO=<ID=NumRepeatsLongAllele,Number=A,Type=Integer,Description="Number of repeats in the longer allele">
##INFO=<ID=NumRepeatsInReference,Number=1,Type=Float,Description="Number of repeats in the reference">
##INFO=<ID=IsPureRepeat,Number=1,Type=String,Description="Whether the repeat is pure">
##INFO=<ID=MotifInterruptionIndex,Number=A,Type=Integer,Description="Index of motif interruption">
EOT
        
        bcftools view -h ${FINAL_PREFIX}.vcf.gz > current_header.txt
        head -n -1 current_header.txt > new_full_header.txt
        cat new_header_fields.txt >> new_full_header.txt
        tail -n 1 current_header.txt >> new_full_header.txt
        
        bcftools reheader -h new_full_header.txt ${FINAL_PREFIX}.vcf.gz -o ${FINAL_PREFIX}.reheadered.vcf.gz
        mv ${FINAL_PREFIX}.reheadered.vcf.gz ${FINAL_PREFIX}.vcf.gz

        tabix -f ${FINAL_PREFIX}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}~{output_suffix}.vcf.gz"
        File filtered_vcf_index = "~{prefix}~{output_suffix}.vcf.gz.tbi"
        File variants_tsv = "~{prefix}~{output_suffix}.variants.tsv.gz"
        File alleles_tsv = "~{prefix}~{output_suffix}.alleles.tsv.gz"
        File variants_bed = "~{prefix}~{output_suffix}.variants.bed.gz"
        File variants_bed_index = "~{prefix}~{output_suffix}.variants.bed.gz.tbi"
        File filter_log = "~{prefix}~{output_suffix}.filter_vcf.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(size(vcf, "GB")) + 5,
        disk_gb: ceil(size(vcf, "GB") + size(ref_fasta, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
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

task PostprocessVcf {
    input {
        String prefix
        File vcf
        File vcf_index
        String output_suffix
        String pipeline_docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -exuo pipefail

        bcftools annotate --set-id '%INFO/ID' ~{vcf} -Oz -o ~{prefix}~{output_suffix}.restored_id.vcf.gz
        tabix ~{prefix}~{output_suffix}.restored_id.vcf.gz
    >>>

    output {
        File restored_id_vcf = "~{prefix}~{output_suffix}.restored_id.vcf.gz"
        File restored_id_vcf_index = "~{prefix}~{output_suffix}.restored_id.vcf.gz.tbi"
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

task AnnotateOriginalVcf {
    input {
        File original_vcf
        File original_vcf_index
        File str_vcf
        File str_vcf_index
        String prefix
        String pipeline_docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -exuo pipefail

        bcftools annotate \
            -a ~{str_vcf} \
            -c "INFO/LocusId,INFO/Locus,INFO/Motif,INFO/NumRepeatsShortAllele,INFO/NumRepeatsLongAllele,INFO/NumRepeatsInReference,INFO/IsPureRepeat,INFO/MotifInterruptionIndex" \
            -Oz -o "~{prefix}.final.vcf.gz" \
            ~{original_vcf}

        tabix "~{prefix}.final.vcf.gz"
    >>>

    output {
        File final_vcf = "~{prefix}.final.vcf.gz"
        File final_vcf_index = "~{prefix}.final.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(original_vcf, "GB") + size(str_vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
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