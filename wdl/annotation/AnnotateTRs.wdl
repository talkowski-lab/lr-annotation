version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx

        String prefix
        String tr_info
        String tr_filter
        String? tr_rename_ids_string
        Array[String] contigs

        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tr_vcf
        RuntimeAttr? runtime_attr_rename_variant_ids
        RuntimeAttr? runtime_attr_annotate_trs
        RuntimeAttr? runtime_attr_concat_vcf
    }

    call CheckSamplesMatch {
        input:
            vcf = vcf,
            tr_vcf = tr_vcf,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetVcf {
            input:
                vcf = vcf,
                vcf_index = vcf_idx,
                contig = contig,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.SubsetVcfToContig as SubsetTrVcf {
            input:
                vcf = tr_vcf,
                vcf_index = tr_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.tr",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_tr_vcf
        }

        if (defined(tr_rename_ids_string)) {
            call Helpers.RenameVariantIds {
                input:
                    vcf = SubsetTrVcf.subset_vcf,
                    vcf_index = SubsetTrVcf.subset_vcf_index,
                    id_format = select_first([tr_rename_ids_string]),
                    prefix = "~{prefix}.~{contig}.renamed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_variant_ids
            }
        }

        call AnnotateTRVariants {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_index,
                tr_vcf = select_first([RenameVariantIds.renamed_vcf, SubsetTrVcf.subset_vcf]),
                tr_vcf_idx = select_first([RenameVariantIds.renamed_vcf_index, SubsetTrVcf.subset_vcf_index]),
                tr_info = tr_info,
                tr_filter = tr_filter,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                check_passed = CheckSamplesMatch.check_passed,
                runtime_attr_override = runtime_attr_annotate_trs
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateTRVariants.annotated_vcf,
            vcfs_idx = AnnotateTRVariants.annotated_vcf_idx,
            merge_sort = false,
            prefix = prefix + ".annotated_trs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File tr_annotated_vcf = ConcatVcfs.concat_vcf
        File tr_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task CheckSamplesMatch {
    input {
        File vcf
        File tr_vcf
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} | sort > samples_vcf.txt
        bcftools query -l ~{tr_vcf} | sort > samples_tr_vcf.txt

        if ! diff samples_vcf.txt samples_tr_vcf.txt > /dev/null; then
            echo "ERROR: Sample lists do not match between VCF and TR VCF" >&2
            echo "VCF samples:" >&2
            cat samples_vcf.txt >&2
            echo "TR VCF samples:" >&2
            cat samples_tr_vcf.txt >&2
            exit 1
        fi

        echo "Sample lists match"
        echo "true" > check_passed.txt
    >>>

    output {
        Boolean check_passed = read_boolean("check_passed.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 5,
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
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateTRVariants {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        String tr_info
        String tr_filter
        String prefix
        String docker
        Boolean check_passed
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > sample_order.txt

        bcftools view -S sample_order.txt ~{tr_vcf} -Oz -o tr_reordered.vcf.gz
        tabix -p vcf tr_reordered.vcf.gz

        echo '##INFO=<ID=TR_CALLER,Number=1,Type=String,Description="Tandem repeat caller with variants that envelope this site">' > tr_caller_header.txt
        
        bcftools annotate \
            -h tr_caller_header.txt \
            -Oz -o tr_with_header.vcf.gz \
            tr_reordered.vcf.gz
        tabix -p vcf tr_with_header.vcf.gz

        bcftools view -h tr_with_header.vcf.gz > tr_tagged_header.txt
        bcftools view -H tr_with_header.vcf.gz | \
            awk -v tr_info="~{tr_info}" 'BEGIN{FS=OFS="\t"} {
                if ($8 == "." || $8 == "") {
                    $8 = "TR_CALLER=" tr_info
                } else {
                    $8 = $8 ";TR_CALLER=" tr_info
                }
                print
            }' | \
            cat tr_tagged_header.txt - | \
            bgzip -c > tr_tagged.vcf.gz
        tabix -p vcf tr_tagged.vcf.gz

        cat > header_additions.txt <<'HEADER_EOF'
##INFO=<ID=TR_CALLER,Number=1,Type=String,Description="Tandem repeat caller with variants that envelope this site">
##INFO=<ID=TR_CALLER_ID,Number=1,Type=String,Description="IDs of the tandem repeat variants that envelope this site">
##FILTER=<ID=~{tr_filter},Description="Variant is enveloped by a tandem repeat region">
HEADER_EOF

        bcftools view -h ~{tr_vcf} | \
            grep -E "^##(INFO|FORMAT|FILTER)=" | \
            grep -v "^##INFO=<ID=TR_CALLER," | \
            grep -v "^##FILTER=<ID=~{tr_filter}," \
            >> header_additions.txt || true
        
        bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\n' tr_tagged.vcf.gz > tr_regions.bed
        bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\n' ~{vcf} > vcf_regions.bed

        bedtools intersect \
            -a vcf_regions.bed \
            -b tr_regions.bed \
            -f 1.0 \
            -wa \
            -wb > enveloped_variants.bed
        
        if [ -s enveloped_variants.bed ]; then
            awk -v filter="~{tr_filter}" 'BEGIN{OFS="\t"} {print $1, $2, $5, $6, filter, $10}' \
                enveloped_variants.bed > filter_annotations.txt
        else
            touch filter_annotations.txt
        fi

        if [ -s filter_annotations.txt ]; then
            bgzip filter_annotations.txt
            
            tabix -s 1 -b 2 -e 2 filter_annotations.txt.gz
            
            bcftools annotate \
                -a filter_annotations.txt.gz \
                -c CHROM,POS,REF,ALT,FILTER,INFO/TR_CALLER_ID \
                -h header_additions.txt \
                -Oz -o vcf_with_filters.vcf.gz \
                ~{vcf}
            
            tabix -p vcf vcf_with_filters.vcf.gz
        else
            cp ~{vcf} vcf_with_filters.vcf.gz
            cp ~{vcf_idx} vcf_with_filters.vcf.gz.tbi
        fi

        bcftools annotate \
            -h header_additions.txt \
            -Oz -o tr_with_headers.vcf.gz \
            tr_tagged.vcf.gz
        tabix -p vcf tr_with_headers.vcf.gz

        bcftools concat \
            --allow-overlaps \
            vcf_with_filters.vcf.gz \
            tr_with_headers.vcf.gz \
            -Oz -o ~{prefix}.annotated.vcf.gz
                
        tabix -p vcf ~{prefix}.annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 5,
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
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
