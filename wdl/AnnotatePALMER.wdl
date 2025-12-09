version 1.0

import "general/Structs.wdl"

workflow AnnotatePALMER {
    input {
        File vcf
        File vcf_idx

        String prefix
        File PALMER_vcf
        File PALMER_vcf_idx

        File rm_out
        Int rm_buffer = 50
        File ref_fai

        Float reciprocal_overlap_ALU = 0.0
        Float reciprocal_overlap_SVA = 0.0
        Float reciprocal_overlap_LINE = 0.0
        Float reciprocal_overlap_HERVK = 0.0

        Float size_similarity_ALU = 0.9
        Float size_similarity_SVA = 0.9
        Float size_similarity_LINE = 0.9
        Float size_similarity_HERVK = 0.9

        Float sequence_similarity_ALU = 0.9
        Float sequence_similarity_SVA = 0.9
        Float sequence_similarity_LINE = 0.9
        Float sequence_similarity_HERVK = 0.9

        Int breakpoint_window_ALU = 100000
        Int breakpoint_window_SVA = 100000
        Int breakpoint_window_LINE = 100000
        Int breakpoint_window_HERVK = 100000

        Int min_shared_samples_ALU = 0
        Int min_shared_samples_SVA = 0
        Int min_shared_samples_LINE = 0
        Int min_shared_samples_HERVK = 0

        String annotate_palmer_docker

        RuntimeAttr? runtime_attr_filter_palmer
    }

    call FilterPALMER {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix,
            PALMER_vcf = PALMER_vcf,
            PALMER_vcf_idx = PALMER_vcf_idx,
            rm_out = rm_out,
            rm_buffer = rm_buffer,
            ref_fai = ref_fai,
            docker = annotate_palmer_docker,
            runtime_attr_override = runtime_attr_filter_palmer,
            reciprocal_overlap_ALU = reciprocal_overlap_ALU,
            reciprocal_overlap_LINE = reciprocal_overlap_LINE,
            reciprocal_overlap_SVA = reciprocal_overlap_SVA,
            reciprocal_overlap_HERVK = reciprocal_overlap_HERVK,
            size_similarity_ALU = size_similarity_ALU,
            size_similarity_LINE = size_similarity_LINE,
            size_similarity_SVA = size_similarity_SVA,
            size_similarity_HERVK = size_similarity_HERVK,
            sequence_similarity_ALU = sequence_similarity_ALU,
            sequence_similarity_LINE = sequence_similarity_LINE,
            sequence_similarity_SVA = sequence_similarity_SVA,
            sequence_similarity_HERVK = sequence_similarity_HERVK,
            breakpoint_window_ALU = breakpoint_window_ALU,
            breakpoint_window_LINE = breakpoint_window_LINE,
            breakpoint_window_SVA = breakpoint_window_SVA,
            breakpoint_window_HERVK = breakpoint_window_HERVK,
            min_shared_samples_ALU = min_shared_samples_ALU,
            min_shared_samples_LINE = min_shared_samples_LINE,
            min_shared_samples_SVA = min_shared_samples_SVA,
            min_shared_samples_HERVK = min_shared_samples_HERVK
    }

    output {
        File palmer_annotated_vcf = FilterPALMER.annotated_vcf
        File palmer_annotated_vcf_idx = FilterPALMER.annotated_vcf_idx
    }
}

task FilterPALMER {
    input {
        File vcf
        File vcf_idx
        String prefix
        File PALMER_vcf
        File PALMER_vcf_idx
        File rm_out
        Int rm_buffer
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override

        Float reciprocal_overlap_ALU
        Float reciprocal_overlap_LINE
        Float reciprocal_overlap_SVA
        Float reciprocal_overlap_HERVK

        Float size_similarity_ALU
        Float size_similarity_LINE
        Float size_similarity_SVA
        Float size_similarity_HERVK

        Float sequence_similarity_ALU
        Float sequence_similarity_LINE
        Float sequence_similarity_SVA
        Float sequence_similarity_HERVK

        Int breakpoint_window_ALU
        Int breakpoint_window_LINE
        Int breakpoint_window_SVA
        Int breakpoint_window_HERVK

        Int min_shared_samples_ALU
        Int min_shared_samples_LINE
        Int min_shared_samples_SVA
        Int min_shared_samples_HERVK
    }

    command <<<
        set -euo pipefail

        cur_vcf=~{vcf}
        cut -f1,2 ~{ref_fai} > genome_file

        types=("ALU" "LINE" "SVA" "HERVK")

        for ME_type in "${types[@]}"; do
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
                reciprocal_overlap_threshold=~{reciprocal_overlap_ALU}
                size_similarity_threshold=~{size_similarity_ALU}
                sequence_similarity_threshold=~{sequence_similarity_ALU}
                breakpoint_window=~{breakpoint_window_ALU}
                min_shared_samples=~{min_shared_samples_ALU}
            elif [[ "${ME_type}" == "LINE" ]]; then
                reciprocal_overlap_threshold=~{reciprocal_overlap_LINE}
                size_similarity_threshold=~{size_similarity_LINE}
                sequence_similarity_threshold=~{sequence_similarity_LINE}
                breakpoint_window=~{breakpoint_window_LINE}
                min_shared_samples=~{min_shared_samples_LINE}
            elif [[ "${ME_type}" == "SVA" ]]; then
                reciprocal_overlap_threshold=~{reciprocal_overlap_SVA}
                size_similarity_threshold=~{size_similarity_SVA}
                sequence_similarity_threshold=~{sequence_similarity_SVA}
                breakpoint_window=~{breakpoint_window_SVA}
                min_shared_samples=~{min_shared_samples_SVA}
            elif [[ "${ME_type}" == "HERVK" ]]; then
                reciprocal_overlap_threshold=~{reciprocal_overlap_HERVK}
                size_similarity_threshold=~{size_similarity_HERVK}
                sequence_similarity_threshold=~{sequence_similarity_HERVK}
                breakpoint_window=~{breakpoint_window_HERVK}
                min_shared_samples=~{min_shared_samples_HERVK}
            fi

            awk '$11==FILTER' FILTER="${MEfilter}" ~{rm_out} | \
                awk 'OFS="\t" {print $5,$7,$8}'| sed 's/(//'|sed 's/)//'|awk 'OFS="\t" {print $1,$2+$3}' | \
                sed 's/:/\t/'|sed 's/;/\t/' | awk 'OFS="\t" {print $1,$2-1,$2,$2,$3,$4}'| sort -k1,1 -k2,2n | uniq | \
                bedtools slop -g genome_file -b ~{rm_buffer} | \
                bedtools merge -c 4,5,6 -o collapse > RM_filtered.bed

            bcftools query -f '%CHROM\t%POS\t%ID\t%SVLEN\t[%SAMPLE,]\t%REF\t%ALT\n' -i 'GT=="alt"' ${ME_type}_subset.vcf.gz \
                | awk 'OFS="\t" {print $1,$2-1,$2,$3,$4,$5,$6,$7}' \
                | sort -k1,1 -k2,2n \
                > PALMER_calls.bed

            bedtools intersect -wo -a RM_filtered.bed -b PALMER_calls.bed > intersection

            python /opt/gnomad-lr/scripts/mei/PALMER_transfer_annotations.py \
                --intersection intersection \
                --target-vcf ${cur_vcf} \
                --me-type ${ME_type} \
                --output ${ME_type}_annotations.raw.vcf.gz \
                --reciprocal-overlap "${reciprocal_overlap_threshold}" \
                --size-similarity "${size_similarity_threshold}" \
                --sequence-similarity "${sequence_similarity_threshold}" \
                --breakpoint-window "${breakpoint_window}" \
                --min-shared-samples "${min_shared_samples}"

            gunzip -c ${ME_type}_annotations.raw.vcf.gz > ${ME_type}_annotations.raw.vcf
            grep '^#' ${ME_type}_annotations.raw.vcf > annotations_to_transfer.vcf
            (grep -v '^#' ${ME_type}_annotations.raw.vcf || true) | sort -k1,1 -k2,2n -u >> annotations_to_transfer.vcf
            bgzip -f annotations_to_transfer.vcf
            tabix -f -p vcf annotations_to_transfer.vcf.gz
            rm -f ${ME_type}_annotations.raw.vcf ${ME_type}_annotations.raw.vcf.gz

            echo '##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">' > line.header
            echo '##INFO=<ID=ME_LEN,Number=.,Type=Integer,Description="Length of SV">' >> line.header
            echo '##INFO=<ID=ME_ID,Number=.,Type=String,Description="ID of the PALMER variant">' >> line.header

            bcftools annotate \
                -Oz \
                -a annotations_to_transfer.vcf.gz \
                -c CHROM,POS,REF,ALT,+ME_TYPE,+ME_LEN,+ME_ID \
                -h line.header \
                ${cur_vcf} \
                > tmp.vcf.gz
            cp tmp.vcf.gz ~{prefix}.anno.vcf.gz
            tabix -f ~{prefix}.anno.vcf.gz
            cur_vcf=~{prefix}.anno.vcf.gz
        done
    >>>

    output {
        File annotated_vcf = "~{prefix}.anno.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.anno.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3*ceil(size(vcf, "GB")+size(PALMER_vcf, "GB")+size(rm_out, "GB"))+20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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
