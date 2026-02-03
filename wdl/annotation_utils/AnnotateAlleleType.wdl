version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateAlleleType {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[String] contigs
        String prefix

        Array[String]? annotations_prefixes
        Array[String]? annotations_suffixes

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
                prefix = "~{prefix}.~{contig}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateSequentially.annotated_vcf,
            vcfs_idx = AnnotateSequentially.annotated_vcf_idx,
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
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > prefixes.txt <<EOF
~{write_lines(select_first([annotations_prefixes, []]))}
EOF

        cat > suffixes.txt <<EOF
~{write_lines(select_first([annotations_suffixes, []]))}
EOF

        current_vcf="~{vcf}"
        i=0
        for tsv_file in ~{sep=' ' annotations_tsvs}; do
            prefix_val=$(sed -n "$((i + 1))p" prefixes.txt || echo "")
            suffix_val=$(sed -n "$((i + 1))p" suffixes.txt || echo "")
            
            awk -v pre="$prefix_val" -v suf="$suffix_val" 'BEGIN{FS=OFS="\t"} {$6=pre $6 suf; print}' "$tsv_file" > "modified_${i}.tsv"
            bgzip -c "modified_${i}.tsv" > "annotations_${i}.tsv.gz"
            tabix -s1 -b2 -e2 "annotations_${i}.tsv.gz"
            
            bcftools annotate \
                -a "annotations_${i}.tsv.gz" \
                -c CHROM,POS,REF,ALT,~ID,INFO/allele_type \
                -Oz -o "temp_${i}.vcf.gz" \
                "$current_vcf"
            current_vcf="temp_${i}.vcf.gz"

            i=$((i + 1))
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
