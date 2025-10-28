version 1.0

import "general/Structs.wdl"

workflow AnnotatePALMER {
    input {
        File vcf
        File vcf_idx

        String prefix
        Array[File] PALMER_vcfs
        Array[String] ME_types

        File rm_out
        File rm_fa
        File ref_fai

        String palmer_docker

        RuntimeAttr? runtime_attr_filter_palmer
    }

    call FilterPALMER {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix,
            PALMER_vcfs = PALMER_vcfs,
            ME_types = ME_types,
            rm_out = rm_out,
            rm_fa = rm_fa,
            ref_fai = ref_fai,
            docker = palmer_docker,
            runtime_attr_override = runtime_attr_filter_palmer
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
        Array[File] PALMER_vcfs
        Array[String] ME_types
        File rm_out
        File rm_fa
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cur_vcf=~{vcf}
        cut -f1,2 ~{ref_fai} > genome_file
        
        call_files=(~{sep=" " PALMER_vcfs})
        types=(~{sep=" " ME_types})

        for i in "${!types[@]}"; do
            ME_type=${types[i]}
            PALMER_vcf=${call_files[i]}
            echo "Processing ${ME_type} calls from ${PALMER_vcf}"

            MEfilter=""
            if [[ "${ME_type}" == "LINE" ]]; then MEfilter="LINE/L1";
            elif [[ "${ME_type}" == "ALU" ]]; then MEfilter="SINE/Alu";
            elif [[ "${ME_type}" == "SVA" ]]; then MEfilter="Retroposon/SVA";
            elif [[ "${ME_type}" == "HERVK" ]]; then MEfilter="LTR/ERVK";
            fi

            if [ -z "$MEfilter" ]; then echo "ERROR: filter type is unset (check that MEI type is one of permitted options)"; exit 1; fi

            awk '$11==FILTER' FILTER="${MEfilter}" ~{rm_out} | \
                awk 'OFS="\t" {print $5,$7,$8}'| sed 's/(//'|sed 's/)//'|awk 'OFS="\t" {print $1,$2+$3}' | \
                sed 's/:/\t/'|sed 's/;/\t/' | awk 'OFS="\t" {print $1,$2-1,$2,$2,$3,$4}'| sort -k1,1 -k2,2n | uniq | \
                bedtools slop -g genome_file -b 50 | \
                bedtools merge -c 4,5,6 -o collapse > RMfilter.50bpbuffer.bed

            bcftools query -f '%CHROM\t%POS\t%SVLEN\t[%SAMPLE,]\n' -i 'GT=="alt"' ${PALMER_vcf} \
                | awk 'OFS="\t" {print $1,$2-1,$2,$3,$4}' \
                | sort -k1,1 -k2,2n \
                > PALMER_calls.bed

            bedtools intersect -wo -a RMfilter.50bpbuffer.bed -b PALMER_calls.bed > intersection

            python /opt/gnomad-lr/scripts/palmer/PALMER_filter_and_transfer_annotations.py \
                intersection \
                ~{vcf} \
                ~{rm_fa} \
                ${ME_type} \
                | sort -k1,1 -k2,2n \
                | uniq \
                | bgzip \
                > annotations_to_transfer.tsv.gz
            tabix -f -b 2 -e 2 annotations_to_transfer.tsv.gz

            echo '##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">' > line.header
            echo '##INFO=<ID=ME_LEN,Number=.,Type=Integer,Description="Length of SV">' >> line.header

            bcftools annotate \
                -Oz \
                -a annotations_to_transfer.tsv.gz \
                -c CHROM,POS,REF,ALT,+ME_TYPE,+ME_LEN \
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
        disk_gb: 5*ceil(size(PALMER_vcfs, "GB")+size(rm_out, "GB")+size(rm_fa, "GB")+size(ref_fai, "GB")+size(vcf, "GB"))+20,
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
