version 1.0

import "Structs.wdl"

workflow AnnotatePALMER {
    input {
        File SV_vcf
        File SV_vcf_tbi

        Array[File] PALMER_vcfs
        Array[String] ME_types

        File RM_out
        File RM_fa
        File ref_fai

        RuntimeAttr? runtime_attr_override_filter_palmer
    }

    call FilterPALMER { 
        input: 
            PALMER_vcfs=PALMER_vcfs, 
            RM_out=RM_out, 
            RM_fa=RM_fa, 
            ref_fai=ref_fai, 
            ME_types=ME_types, 
            SV_vcf=SV_vcf, 
            SV_vcf_tbi=SV_vcf_tbi,
            runtime_attr_override=runtime_attr_override_filter_palmer
    }

    output {
        File palmer_annotated_vcf = FilterPALMER.vcf_anno
        File palmer_annotated_vcf_index = FilterPALMER.vcf_anno_tbi
    }
}

task FilterPALMER {
    input {
        Array[File] PALMER_vcfs
        File RM_out
        File RM_fa
        File ref_fai
        File SV_vcf
        File SV_vcf_tbi
        Array[String] ME_types

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 5*ceil(size(PALMER_vcfs, "GB")+size(RM_out, "GB")+size(RM_fa, "GB")+size(ref_fai, "GB")+size(SV_vcf, "GB"))+20

    String vcfbase = basename(SV_vcf, ".vcf.gz")

    command <<<
        set -euxo pipefail

        cur_vcf=~{SV_vcf}
        cut -f1,2 ~{ref_fai} > genome_file
        
        call_files=(~{sep=" " PALMER_vcfs})
        types=(~{sep=" " ME_types})

        for i in "${!types[@]}"; do
            ME_type=${types[i]}
            PALMER_vcf=${call_files[i]}
            echo "Processing ${ME_type} calls from ${PALMER_vcf}"

            # process repeatmasker outputs
            MEfilter=""
            if [[ "${ME_type}" == "LINE" ]]; then MEfilter="LINE/L1";
            elif [[ "${ME_type}" == "ALU" ]]; then MEfilter="SINE/Alu";
            elif [[ "${ME_type}" == "SVA" ]]; then MEfilter="Retroposon/SVA";
            elif [[ "${ME_type}" == "HERVK" ]]; then MEfilter="LTR/ERVK";
            fi

            # check that filter is set - throw error if not
            if [ -z "$MEfilter" ]; then echo "ERROR: filter type is unset (check that MEI type is one of permitted options)"; exit 1; fi

            # create bed file with 50bp buffer around INS positions that match the ME type per repeatmasker
            awk '$11==FILTER' FILTER="${MEfilter}"  ~{RM_out} | \
                awk 'OFS="\t" {print $5,$7,$8}'| sed 's/(//'|sed 's/)//'|awk 'OFS="\t" {print $1,$2+$3}' | \
                sed 's/:/\t/'|sed 's/;/\t/' | awk 'OFS="\t" {print $1,$2-1,$2,$2,$3,$4}'| sort -k1,1 -k2,2n | uniq | \
                bedtools slop -g genome_file -b 50 | \
                bedtools merge -c 4,5,6 -o collapse > RMfilter.50bpbuffer.bed

            bcftools query -f '%CHROM\t%POS\t%SVLEN\t[%SAMPLE,]\n' -i 'GT=="alt"' ${PALMER_vcf} |awk 'OFS="\t" {print $1,$2-1,$2,$3,$4}' |sort -k1,1 -k2,2n > PALMER_calls.bed

            bedtools intersect -wo -a RMfilter.50bpbuffer.bed -b PALMER_calls.bed > intersection

            python /PALMER_filter_and_transfer_annotations.py intersection ~{SV_vcf} ~{RM_fa} ${ME_type} | sort -k1,1 -k2,2n | uniq | bgzip > annotations_to_transfer.tsv.gz
            tabix -f -b 2 -e 2 annotations_to_transfer.tsv.gz

            echo '##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">' > line.header
            echo '##INFO=<ID=ME_LEN,Number=.,Type=Integer,Description="Length of SV">' >> line.header

            bcftools annotate -Oz -a annotations_to_transfer.tsv.gz -c CHROM,POS,REF,ALT,+ME_TYPE,+ME_LEN -h line.header ${cur_vcf} > tmp.vcf.gz
            cp tmp.vcf.gz ~{vcfbase}.anno.vcf.gz
            tabix -f ~{vcfbase}.anno.vcf.gz
            cur_vcf=~{vcfbase}.anno.vcf.gz
        done
    >>>

    output {
        File vcf_anno = "~{vcfbase}.anno.vcf.gz"
        File vcf_anno_tbi = "~{vcfbase}.anno.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-palmer-filter:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
