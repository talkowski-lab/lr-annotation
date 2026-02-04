version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        Array[String] contigs
        String prefix

        String tr_caller
        String? tr_rename_ids_string

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
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.SubsetVcfToContig as SubsetTrVcf {
            input:
                vcf = tr_vcf,
                vcf_idx = tr_vcf_idx,
                contig = contig,
                extra_args = "--min-ac 1",
                prefix = "~{prefix}.~{contig}.tr",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_tr_vcf
        }

        if (defined(tr_rename_ids_string)) {
            call Helpers.RenameVariantIds {
                input:
                    vcf = SubsetTrVcf.subset_vcf,
                    vcf_idx = SubsetTrVcf.subset_vcf_idx,
                    id_format = select_first([tr_rename_ids_string]),
                    prefix = "~{prefix}.~{contig}.renamed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_variant_ids
            }
        }

        call AnnotateTRVariants {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                tr_vcf = select_first([RenameVariantIds.renamed_vcf, SubsetTrVcf.subset_vcf]),
                tr_vcf_idx = select_first([RenameVariantIds.renamed_vcf_idx, SubsetTrVcf.subset_vcf_idx]),
                tr_caller = tr_caller,
                prefix = "~{prefix}.~{contig}.annotated",
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
        File vcf = ConcatVcfs.concat_vcf
        File vcf_idx = ConcatVcfs.concat_vcf_idx
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

task AnnotateTRVariants {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        String tr_caller
        String prefix
        String docker
        Boolean check_passed
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > samples.txt
        bcftools view -S samples.txt ~{tr_vcf} -Oz -o tr_reordered.vcf.gz
        tabix -p vcf tr_reordered.vcf.gz

        echo '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of variant call">' \
            > tr_header.txt
        
        bcftools annotate \
            -h tr_header.txt \
            -Oz -o tr_header_added.vcf.gz \
            tr_reordered.vcf.gz

        bcftools view tr_header_added.vcf.gz \
            | awk -v info="~{tr_caller}" 'BEGIN{OFS="\t"} /^#/ {print; next} {
                $8 = $8 ";SOURCE=" info
                print
            }' | bgzip -c > tr_tagged.vcf.gz
        
        tabix -p vcf tr_tagged.vcf.gz

        bcftools query \
            -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\t%INFO/allele_type\t%INFO/allele_length\n' \
            ~{vcf} \
            | awk '$7 != "SNV" && $7 != "."' \
            > vcf.bed

        bcftools query \
            -f '%CHROM\t%POS\t%END\t%ID\t%REF\t%ALT\t%INFO/MOTIFS\n' \
            tr_tagged.vcf.gz \
            > tr.bed

        bedtools intersect \
            -f 1.0 \
            -wa \
            -wb \
            -a vcf.bed \
            -b tr.bed \
            > overlaps.bed

        awk -v filter="~{tr_caller}_OVERLAPPED" -v source="~{tr_caller}" 'BEGIN{OFS="\t"} {
            split($15, motifs, ",")
            min_motif = 1000000
            for (i in motifs) {
                len = length(motifs[i])
                if (len < min_motif) min_motif = len
            }

            allele_length = ($8 < 0) ? -$8 : $8
            if (allele_length < min_motif) next

            print $1, $2, $5, $6, filter, $12, source
        }' overlaps.bed | sort -k1,1 -k2,2n | bgzip -c > annotations.tsv.gz

        tabix -s 1 -b 2 -e 2 annotations.tsv.gz

        bcftools view \
            -h ~{vcf} \
        | grep -v 'ID=AL,' > vcf_no_al_header.txt

        bcftools view \
            -h ~{tr_vcf} \
        | grep 'ID=AL,' > al_header.txt || true

        bcftools reheader \
            -h vcf_no_al_header.txt \
            ~{vcf} \
            > vcf_stripped.vcf.gz

        cat <<EOF > new_header.txt
##INFO=<ID=TRID,Number=1,Type=String,Description="ID of enveloping tandem repeat">
##FILTER=<ID=~{tr_caller}_OVERLAPPED,Description="Variant enveloped by ~{tr_caller}">
EOF
        
        cat new_header.txt al_header.txt > merged_headers.txt

        bcftools annotate \
            -a annotations.tsv.gz \
            -c CHROM,POS,REF,ALT,FILTER,INFO/TRID,INFO/SOURCE \
            -h merged_headers.txt \
            -Oz -o vcf_annotated.vcf.gz \
            vcf_stripped.vcf.gz

        tabix -p vcf vcf_annotated.vcf.gz

        bcftools concat \
            --allow-overlaps \
            vcf_annotated.vcf.gz \
            tr_tagged.vcf.gz \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 20,
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
