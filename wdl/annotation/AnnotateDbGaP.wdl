version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateDbGaP {
    input {
        File vcf
        File vcf_idx
        File dbgap_vcf
        File dbgap_vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_dbgap_vcf
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetVcf {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.SubsetVcfToContig as SubsetDbGaPVcf {
            input:
                vcf = dbgap_vcf,
                vcf_idx = dbgap_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.dbgap",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_dbgap_vcf
        }

        call AnnotateDbGaPIds {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                dbgap_vcf = SubsetDbGaPVcf.subset_vcf,
                dbgap_vcf_idx = SubsetDbGaPVcf.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateDbGaPIds.annotated_vcf,
            vcf_idxs = AnnotateDbGaPIds.annotated_vcf_idx,
            prefix = prefix + ".dbgap_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File dbgap_annotated_vcf = ConcatVcfs.concat_vcf
        File dbgap_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AnnotateDbGaPIds {
    input {
        File vcf
        File vcf_idx
        File dbgap_vcf
        File dbgap_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch new_headers.txt
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=dbGaP_ID'; then
            echo '##INFO=<ID=dbGaP_ID,Number=1,Type=String,Description="dbGaP variant ID for matching SNV">' >> new_headers.txt
        fi

        bcftools annotate \
            -h new_headers.txt \
            -Oz -o temp.vcf.gz \
            ~{vcf}
        
        tabix -p vcf temp.vcf.gz

        # Extract SNVs from input VCF for variants with allele_type='snv'
        bcftools query \
            -i 'INFO/allele_type="snv"' \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            temp.vcf.gz \
            > snvs_from_input.txt

        # Extract all variants from dbGaP VCF
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            ~{dbgap_vcf} \
            > dbgap_variants.txt

        # Match SNVs with dbGaP variants on CHROM, POS, REF and ALT
        awk 'BEGIN{OFS="\t"} 
        NR==FNR {
            key = $1"\t"$2"\t"$3"\t"$4
            dbgap[key] = $5
            next
        } 
        {
            key = $1"\t"$2"\t"$3"\t"$4
            if (key in dbgap) {
                print $1, $2, $3, $4, $5, dbgap[key]
            }
        }' dbgap_variants.txt snvs_from_input.txt \
            | bgzip -c > annotations.txt.gz
        
        tabix -s1 -b2 -e2 annotations.txt.gz

        # Annotate VCF with dbGaP_ID
        bcftools annotate \
            -a annotations.txt.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/dbGaP_ID \
            -Oz -o ~{prefix}.vcf.gz \
            temp.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, dbgap_vcf], "GB")) + 10,
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
