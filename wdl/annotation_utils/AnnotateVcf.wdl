version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateVcf {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[String] contigs
        String prefix

        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers

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
                info_names = info_names,
                info_descriptions = info_descriptions,
                info_types = info_types,
                info_numbers = info_numbers,
                prefix = "~{prefix}.~{contig}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateSequentially.annotated_vcf,
            vcfs_idx = AnnotateSequentially.annotated_vcf_idx,
            prefix = prefix + ".annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File annotated_vcf = ConcatVcfs.concat_vcf
        File annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AnnotateSequentially {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<EOF
import sys

info_names = [line.strip().split('\t') for line in open('~{write_tsv(info_names)}')]
info_descriptions = [line.strip().split('\t') for line in open('~{write_tsv(info_descriptions)}')]
info_types = [line.strip().split('\t') for line in open('~{write_tsv(info_types)}')]
info_numbers = [line.strip().split('\t') for line in open('~{write_tsv(info_numbers)}')]

if len(info_names) != len(info_descriptions) or len(info_names) != len(info_types) or len(info_names) != len(info_numbers):
    sys.stderr.write("Error: All info arrays must have the same length\n")
    sys.exit(1)

for i, (names, descs, types, numbers) in enumerate(zip(info_names, info_descriptions, info_types, info_numbers)):
    if len(names) != len(descs) or len(names) != len(types) or len(names) != len(numbers):
        sys.stderr.write(f"Error: info arrays at index {i} must have the same length\n")
        sys.exit(1)
    
    with open(f"header_{i}.txt", "w") as f:
        for name, desc, type_val, number in zip(names, descs, types, numbers):
            f.write(f'##INFO=<ID={name},Number={number},Type={type_val},Description="{desc}">\n')
    
    column_spec = ','.join(['CHROM', 'POS', 'REF', 'ALT', '~ID'] + [f'INFO/{name}' for name in names])
    with open(f"columns_{i}.txt", "w") as f:
        f.write(column_spec)

with open("num_tsvs.txt", "w") as f:
    f.write(str(len(info_names)))
EOF

        current_vcf="~{vcf}"
        
        i=0
        for tsv_file in ~{sep=' ' annotations_tsvs}; do
            bgzip -c "$tsv_file" > "annotations_${i}.tsv.gz"
            tabix -s1 -b2 -e2 "annotations_${i}.tsv.gz"
            
            COLUMN_SPEC=$(cat "columns_${i}.txt")
            
            bcftools annotate \
                -a "annotations_${i}.tsv.gz" \
                -h "header_${i}.txt" \
                -c "$COLUMN_SPEC" \
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
