version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotatePALMER {
    input {
        File vcf
        File vcf_idx
        File PALMER_vcf
        File PALMER_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] mei_types
        Int min_length
        File rm_out
        Int rm_buffer
        File ref_fai
        Float ins_reciprocal_overlap_ALU = 0.9
        Float ins_reciprocal_overlap_SVA = 0.9
        Float ins_reciprocal_overlap_LINE = 0.9
        Float ins_reciprocal_overlap_HERVK = 0.9
        Float ins_size_similarity_ALU = 0.9
        Float ins_size_similarity_SVA = 0.9
        Float ins_size_similarity_LINE = 0.9
        Float ins_size_similarity_HERVK = 0.9
        Float ins_sequence_similarity_ALU = 0.9
        Float ins_sequence_similarity_SVA = 0.9
        Float ins_sequence_similarity_LINE = 0.9
        Float ins_sequence_similarity_HERVK = 0.9
        Int ins_breakpoint_window_ALU = 500
        Int ins_breakpoint_window_SVA = 500
        Int ins_breakpoint_window_LINE = 500
        Int ins_breakpoint_window_HERVK = 500
        Int ins_min_shared_samples_ALU = 0
        Int ins_min_shared_samples_SVA = 0
        Int ins_min_shared_samples_LINE = 0
        Int ins_min_shared_samples_HERVK = 0

        String annotate_palmer_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_filter_palmer
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByLength {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                min_length = min_length,
                extra_args = if single_contig then "" else "--regions ~{contig}",
                prefix = "~{prefix}.~{contig}.filtered",
                docker = annotate_palmer_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call FilterPALMER {
            input:
                vcf = SubsetVcfByLength.subset_vcf,
                vcf_idx = SubsetVcfByLength.subset_vcf_idx,
                PALMER_vcf = PALMER_vcf,
                PALMER_vcf_idx = PALMER_vcf_idx,
                mei_types = mei_types,
                rm_out = rm_out,
                rm_buffer = rm_buffer,
                ref_fai = ref_fai,
                ins_reciprocal_overlap_ALU = ins_reciprocal_overlap_ALU,
                ins_reciprocal_overlap_LINE = ins_reciprocal_overlap_LINE,
                ins_reciprocal_overlap_SVA = ins_reciprocal_overlap_SVA,
                ins_reciprocal_overlap_HERVK = ins_reciprocal_overlap_HERVK,
                ins_size_similarity_ALU = ins_size_similarity_ALU,
                ins_size_similarity_LINE = ins_size_similarity_LINE,
                ins_size_similarity_SVA = ins_size_similarity_SVA,
                ins_size_similarity_HERVK = ins_size_similarity_HERVK,
                ins_sequence_similarity_ALU = ins_sequence_similarity_ALU,
                ins_sequence_similarity_LINE = ins_sequence_similarity_LINE,
                ins_sequence_similarity_SVA = ins_sequence_similarity_SVA,
                ins_sequence_similarity_HERVK = ins_sequence_similarity_HERVK,
                ins_breakpoint_window_ALU = ins_breakpoint_window_ALU,
                ins_breakpoint_window_LINE = ins_breakpoint_window_LINE,
                ins_breakpoint_window_SVA = ins_breakpoint_window_SVA,
                ins_breakpoint_window_HERVK = ins_breakpoint_window_HERVK,
                ins_min_shared_samples_ALU = ins_min_shared_samples_ALU,
                ins_min_shared_samples_LINE = ins_min_shared_samples_LINE,
                ins_min_shared_samples_SVA = ins_min_shared_samples_SVA,
                ins_min_shared_samples_HERVK = ins_min_shared_samples_HERVK,
                prefix = "~{prefix}.~{contig}.filtered",
                docker = annotate_palmer_docker,
                runtime_attr_override = runtime_attr_filter_palmer
        }
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotations {
            input:
                tsvs = FilterPALMER.annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.palmer_annotations",
                docker = annotate_palmer_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_palmer = select_first([MergeAnnotations.concatenated_tsv, FilterPALMER.annotations_tsv[0]])
    }
}

task FilterPALMER {
    input {
        File vcf
        File vcf_idx
        File PALMER_vcf
        File PALMER_vcf_idx
        Array[String] mei_types
        File rm_out
        Int rm_buffer
        File ref_fai
        Float ins_reciprocal_overlap_ALU
        Float ins_reciprocal_overlap_LINE
        Float ins_reciprocal_overlap_SVA
        Float ins_reciprocal_overlap_HERVK
        Float ins_size_similarity_ALU
        Float ins_size_similarity_LINE
        Float ins_size_similarity_SVA
        Float ins_size_similarity_HERVK
        Float ins_sequence_similarity_ALU
        Float ins_sequence_similarity_LINE
        Float ins_sequence_similarity_SVA
        Float ins_sequence_similarity_HERVK
        Int ins_breakpoint_window_ALU
        Int ins_breakpoint_window_LINE
        Int ins_breakpoint_window_SVA
        Int ins_breakpoint_window_HERVK
        Int ins_min_shared_samples_ALU
        Int ins_min_shared_samples_LINE
        Int ins_min_shared_samples_SVA
        Int ins_min_shared_samples_HERVK
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cut -f1,2 ~{ref_fai} > genome_file

        mei_types=(~{sep=' ' mei_types})

        for ME_type in "${mei_types[@]}"; do
            bcftools view \
                -i "INFO/ME_TYPE='${ME_type}'" \
                -Oz \
                -o ${ME_type}_subset.vcf.gz \
                ~{PALMER_vcf}
            tabix ${ME_type}_subset.vcf.gz

            MEfilter=""
            if [[ "${ME_type}" == "LINE" ]]; then MEfilter="LINE/L1";
            elif [[ "${ME_type}" == "ALU" ]]; then MEfilter="SINE/Alu";
            elif [[ "${ME_type}" == "SVA" ]]; then MEfilter="Retroposon/SVA";
            elif [[ "${ME_type}" == "HERVK" ]]; then MEfilter="LTR/ERVK";
            fi

            reciprocal_overlap_threshold=0
            size_similarity_threshold=0
            sequence_similarity_threshold=0
            breakpoint_window=0
            min_shared_samples=0
            if [[ "${ME_type}" == "ALU" ]]; then
                reciprocal_overlap_threshold=~{ins_reciprocal_overlap_ALU}
                size_similarity_threshold=~{ins_size_similarity_ALU}
                sequence_similarity_threshold=~{ins_sequence_similarity_ALU}
                breakpoint_window=~{ins_breakpoint_window_ALU}
                min_shared_samples=~{ins_min_shared_samples_ALU}
            elif [[ "${ME_type}" == "LINE" ]]; then
                reciprocal_overlap_threshold=~{ins_reciprocal_overlap_LINE}
                size_similarity_threshold=~{ins_size_similarity_LINE}
                sequence_similarity_threshold=~{ins_sequence_similarity_LINE}
                breakpoint_window=~{ins_breakpoint_window_LINE}
                min_shared_samples=~{ins_min_shared_samples_LINE}
            elif [[ "${ME_type}" == "SVA" ]]; then
                reciprocal_overlap_threshold=~{ins_reciprocal_overlap_SVA}
                size_similarity_threshold=~{ins_size_similarity_SVA}
                sequence_similarity_threshold=~{ins_sequence_similarity_SVA}
                breakpoint_window=~{ins_breakpoint_window_SVA}
                min_shared_samples=~{ins_min_shared_samples_SVA}
            elif [[ "${ME_type}" == "HERVK" ]]; then
                reciprocal_overlap_threshold=~{ins_reciprocal_overlap_HERVK}
                size_similarity_threshold=~{ins_size_similarity_HERVK}
                sequence_similarity_threshold=~{ins_sequence_similarity_HERVK}
                breakpoint_window=~{ins_breakpoint_window_HERVK}
                min_shared_samples=~{ins_min_shared_samples_HERVK}
            fi

            awk '$11==FILTER' FILTER="${MEfilter}" ~{rm_out} | \
                awk 'OFS="\t" {print $5,$7,$8}'| sed 's/(//'|sed 's/)//'|awk 'OFS="\t" {print $1,$2+$3}' | \
                sed 's/:/\t/'|sed 's/;/\t/' | awk 'OFS="\t" {print $1,$2-1,$2,$2,$3,$4}'| sort -k1,1 -k2,2n | uniq | \
                bedtools slop -g genome_file -b ~{rm_buffer} | \
                bedtools merge -c 4,5,6 -o collapse > RM_filtered.bed

            bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/allele_length\t[%SAMPLE,]\t%REF\t%ALT\n' -i 'GT=="alt"' ${ME_type}_subset.vcf.gz \
                | awk 'OFS="\t" {print $1,$2-1,$2,$3,$4,$5,$6,$7}' \
                | sort -k1,1 -k2,2n \
                > PALMER_calls.bed

            bedtools intersect \
                -wo \
                -a RM_filtered.bed \
                -b PALMER_calls.bed \
                > intersection.bed

            python /opt/gnomad-lr/scripts/mei/PALMER_transfer_annotations.py \
                --intersection intersection.bed \
                --target-vcf ~{vcf} \
                --me-type ${ME_type} \
                --output ${ME_type}_annotations.tsv \
                --reciprocal-overlap "${reciprocal_overlap_threshold}" \
                --size-similarity "${size_similarity_threshold}" \
                --sequence-similarity "${sequence_similarity_threshold}" \
                --breakpoint-window "${breakpoint_window}" \
                --min-shared-samples "${min_shared_samples}"
        done

        annotation_files=()
        for ME_type in "${mei_types[@]}"; do
            annotation_files+=("${ME_type}_annotations.tsv")
        done

        cat "${annotation_files[@]}" | sort -k1,1 -k2,2n | uniq > ~{prefix}.palmer_annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.palmer_annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3*ceil(size(vcf, "GB")+size(PALMER_vcf, "GB")+size(rm_out, "GB"))+20,
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
