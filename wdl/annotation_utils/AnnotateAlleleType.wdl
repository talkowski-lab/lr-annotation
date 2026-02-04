version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateAlleleType {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[String] contigs
        String prefix

        Array[String]? annotations_prefixes
        Array[String]? annotations_suffixes
        Array[Boolean]? annotations_lowercase

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        scatter (i in range(length(annotations_tsvs))) {
            call Helpers.SubsetTsvToContig as SubsetTsvs {
                input:
                    tsv = annotations_tsvs[i],
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.tsv~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_tsv
            }
        }

        call AnnotateSequentially {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                annotations_tsvs = SubsetTsvs.subset_tsv,
                annotations_prefixes = annotations_prefixes,
                annotations_suffixes = annotations_suffixes,
                annotations_lowercase = annotations_lowercase,
                prefix = "~{prefix}.~{contig}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateSequentially.annotated_vcf,
            vcf_idxs = AnnotateSequentially.annotated_vcf_idx,
            prefix = prefix + ".allele_type_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File allele_type_annotated_vcf = ConcatVcfs.concat_vcf
        File allele_type_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AnnotateSequentially {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[String]? annotations_prefixes
        Array[String]? annotations_suffixes
        Array[Boolean]? annotations_lowercase
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        current_vcf="~{vcf}"
        annotations_tsvs_array=(~{sep=' ' annotations_tsvs})
        prefixes_array=(~{sep=' ' select_first([annotations_prefixes, []])})
        suffixes_array=(~{sep=' ' select_first([annotations_suffixes, []])})
        lowercase_array=(~{sep=' ' select_first([annotations_lowercase, []])})

        for i in "${!annotations_tsvs_array[@]}"; do
            tsv_file="${annotations_tsvs_array[$i]}"
            prefix_val="${prefixes_array[$i]:-}"
            suffix_val="${suffixes_array[$i]:-}"
            lowercase_val="${lowercase_array[$i]:-false}"
            
            awk -v pre="$prefix_val" -v suf="$suffix_val" -v lower="$lowercase_val" '
                BEGIN {
                    FS = OFS = "\t"
                }
                {
                    gsub(/[\n\r]/, "", $6)
                    $6 = pre $6 suf
                    if (lower == "true") {
                        $6 = tolower($6)
                    }
                    print
                }
            ' "$tsv_file" > "modified_${i}.tsv"
            
            bgzip -c "modified_${i}.tsv" > "annotations_${i}.tsv.gz"
            tabix -s1 -b2 -e2 "annotations_${i}.tsv.gz"
            
            bcftools annotate \
                -a "annotations_${i}.tsv.gz" \
                -c CHROM,POS,REF,ALT,~ID,INFO/allele_type \
                -Oz -o "temp_${i}.vcf.gz" \
                "$current_vcf"
            current_vcf="temp_${i}.vcf.gz"
        done
        
        mv "$current_vcf" ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(annotations_tsvs, "GB")) + 20,
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
