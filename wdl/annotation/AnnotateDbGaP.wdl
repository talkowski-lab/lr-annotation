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
            prefix = "~{prefix}.dbgap_renamed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_dbgap_vcf
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
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
        }

        File contig_vcf = select_first([SubsetVcf.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcf.subset_vcf_idx, vcf_idx])
        File contig_dbgap_vcf = select_first([SubsetDbGaPVcf.subset_vcf, RenameDbGaPContigs.renamed_vcf])
        File contig_dbgap_vcf_idx = select_first([SubsetDbGaPVcf.subset_vcf_idx, RenameDbGaPContigs.renamed_vcf_idx])

        call AnnotateDbGaPIds {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                dbgap_vcf = contig_dbgap_vcf,
                dbgap_vcf_idx = contig_dbgap_vcf_idx,
                prefix = "~{prefix}.~{contig}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = AnnotateDbGaPIds.annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.dbgap_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_dbgap = select_first([ConcatTsvs.concatenated_tsv, AnnotateDbGaPIds.annotations_tsv[0]])
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

        # Extract variants from input VCF
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            ~{vcf} \
            > input_variants.txt

        # Extract all variants from dbGaP VCF
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            ~{dbgap_vcf} \
            > dbgap_variants.txt

        # Match variants with dbGaP variants on CHROM, POS, REF and any matching ALT
        awk 'BEGIN{OFS="\t"} 
        NR==FNR {
            chrom = $1
            pos = $2
            ref = $3
            alt_string = $4
            id = $5
            n = split(alt_string, alts, ",")
            for (i = 1; i <= n; i++) {
                key = chrom "\t" pos "\t" ref "\t" alts[i]
                dbgap[key] = id
            }
            next
        } 
        {
            key = $1 "\t" $2 "\t" $3 "\t" $4
            if (key in dbgap) {
                print $1, $2, $3, $4, $5, dbgap[key]
            }
        }' dbgap_variants.txt input_variants.txt \
            > ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
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
