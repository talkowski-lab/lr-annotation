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

    call RenameDbGaPContigs {
        input:
            vcf = dbgap_vcf,
            vcf_idx = dbgap_vcf_idx,
            prefix = prefix + ".dbgap_renamed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_dbgap_vcf
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
                vcf = RenameDbGaPContigs.renamed_vcf,
                vcf_idx = RenameDbGaPContigs.renamed_vcf_idx,
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

task RenameDbGaPContigs {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > chr_name_mapping.txt <<EOF
NC_000001.11 chr1
NC_000002.12 chr2
NC_000003.12 chr3
NC_000004.12 chr4
NC_000005.10 chr5
NC_000006.12 chr6
NC_000007.14 chr7
NC_000008.11 chr8
NC_000009.12 chr9
NC_000010.11 chr10
NC_000011.10 chr11
NC_000012.12 chr12
NC_000013.11 chr13
NC_000014.9 chr14
NC_000015.10 chr15
NC_000016.10 chr16
NC_000017.11 chr17
NC_000018.10 chr18
NC_000019.10 chr19
NC_000020.11 chr20
NC_000021.9 chr21
NC_000022.11 chr22
NC_000023.11 chrX
NC_000024.10 chrY
EOF

        bcftools annotate \
            --rename-chrs chr_name_mapping.txt \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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
