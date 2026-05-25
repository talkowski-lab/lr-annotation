version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateVcfCleared {
    input {
        File vcf
        File vcf_idx
        File subset_untrimmed_vcf
        File subset_untrimmed_vcf_idx
        Array[File] annotations_tsvs
        Array[String] contigs
        String prefix

        Array[Boolean]? sort_tsvs
        Array[String]? subset_vcf_strings
        Array[String]? awk_tsv_conditions

        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_untrimmed_vcf
        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_clear_annotate_swap
        RuntimeAttr? runtime_attr_concat
    }

    Array[Boolean] sort_tsvs_resolved = select_first([sort_tsvs, []])
    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetMainVcf {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }

            call Helpers.SubsetVcfToContig as SubsetUntrimmedVcf {
                input:
                    vcf = subset_untrimmed_vcf,
                    vcf_idx = subset_untrimmed_vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.untrimmed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_untrimmed_vcf
            }

            scatter (i in range(length(annotations_tsvs))) {
                call Helpers.SubsetTsvToContig {
                    input:
                        tsv = annotations_tsvs[i],
                        contig = contig,
                        sort_output = if length(sort_tsvs_resolved) > 0 then sort_tsvs_resolved[i] else false,
                        prefix = "~{prefix}.~{contig}.tsv~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_tsv
                }
            }
        }

        File contig_vcf = select_first([SubsetMainVcf.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetMainVcf.subset_vcf_idx, vcf_idx])
        File contig_untrimmed_vcf = select_first([SubsetUntrimmedVcf.subset_vcf, subset_untrimmed_vcf])
        File contig_untrimmed_vcf_idx = select_first([SubsetUntrimmedVcf.subset_vcf_idx, subset_untrimmed_vcf_idx])
        Array[File] contig_tsvs = select_first([SubsetTsvToContig.subset_tsv, annotations_tsvs])

        call ClearAnnotateAndSwap {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                subset_untrimmed_vcf = contig_untrimmed_vcf,
                subset_untrimmed_vcf_idx = contig_untrimmed_vcf_idx,
                annotations_tsvs = contig_tsvs,
                subset_vcf_strings = subset_vcf_strings,
                awk_tsv_conditions = awk_tsv_conditions,
                info_names = info_names,
                info_descriptions = info_descriptions,
                info_types = info_types,
                info_numbers = info_numbers,
                prefix = "~{prefix}.~{contig}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_clear_annotate_swap
        }
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = ClearAnnotateAndSwap.annotated_vcf,
                vcf_idxs = ClearAnnotateAndSwap.annotated_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotated_vcf = select_first([ConcatVcfs.concat_vcf, ClearAnnotateAndSwap.annotated_vcf[0]])
        File annotated_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, ClearAnnotateAndSwap.annotated_vcf_idx[0]])
    }
}

task ClearAnnotateAndSwap {
    input {
        File vcf
        File vcf_idx
        File subset_untrimmed_vcf
        File subset_untrimmed_vcf_idx
        Array[File] annotations_tsvs
        Array[String]? subset_vcf_strings
        Array[String]? awk_tsv_conditions
        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Array[String] resolved_subset_strings = select_first([subset_vcf_strings, []])
    Array[String] resolved_awk_conditions = select_first([awk_tsv_conditions, []])

    command <<<
        set -euo pipefail

        # Create annotation files
        python3 <<EOF
import sys

info_names = [line.strip().split('\t') for line in open('~{write_tsv(info_names)}')]
info_descriptions = [line.strip().split('\t') for line in open('~{write_tsv(info_descriptions)}')]
info_types = [line.strip().split('\t') for line in open('~{write_tsv(info_types)}')]
info_numbers = [line.strip().split('\t') for line in open('~{write_tsv(info_numbers)}')]

if len(info_names) != len(info_descriptions) or len(info_names) != len(info_types) or len(info_names) != len(info_numbers):
    sys.stderr.write("Error: All info arrays must have the same length.\n")
    sys.exit(1)

all_info_to_clear = []
for i, (names, descs, types, numbers) in enumerate(zip(info_names, info_descriptions, info_types, info_numbers)):
    if len(names) != len(descs) or len(names) != len(types) or len(names) != len(numbers):
        sys.stderr.write(f"Error: info arrays at index {i} must have the same length.\n")
        sys.exit(1)

    with open(f"header_{i}.txt", "w") as f:
        for name, desc, type_val, number in zip(names, descs, types, numbers):
            f.write(f'##INFO=<ID={name},Number={number},Type={type_val},Description="{desc}">\n')

    column_spec = ','.join(['CHROM', 'POS', 'REF', 'ALT', '~ID'] + [f'INFO/{name}' for name in names])
    with open(f"columns_{i}.txt", "w") as f:
        f.write(column_spec)

    for name in names:
        if name not in all_info_to_clear:
            all_info_to_clear.append(name)

with open("info_to_clear.txt", "w") as f:
    for name in all_info_to_clear:
        f.write(name + "\n")
EOF

        # Clear INFO fields on the trimmed subset VCF
        python3 <<EOF
import pysam

with open("info_to_clear.txt") as f:
    info_to_clear = [line.strip() for line in f if line.strip()]

vcf_in = pysam.VariantFile("~{subset_untrimmed_vcf}")
vcf_out = pysam.VariantFile("cleared_subset.vcf.gz", "wz", header=vcf_in.header)

for rec in vcf_in:
    for name in info_to_clear:
        if name in rec.info:
            del rec.info[name]
    vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
EOF

        tabix -p vcf cleared_subset.vcf.gz

        # Sequentially annotate subset VCF
        current_vcf="cleared_subset.vcf.gz"
        SUBSET_FILE="~{write_lines(resolved_subset_strings)}"
        AWK_COND_FILE="~{write_lines(resolved_awk_conditions)}"
        i=0
        for tsv_file in ~{sep=' ' annotations_tsvs}; do
            AWK_ARG=""
            if [ -s "$AWK_COND_FILE" ]; then
                AWK_ARG=$(sed -n "$((i + 1))p" "$AWK_COND_FILE")
            fi

            if [ -n "$AWK_ARG" ]; then
                awk -F'\t' "$AWK_ARG" "$tsv_file" | bgzip -c > "annotations_${i}.tsv.gz"
            else
                bgzip -c "$tsv_file" > "annotations_${i}.tsv.gz"
            fi
            tabix -s1 -b2 -e2 "annotations_${i}.tsv.gz"

            COLUMN_SPEC=$(cat "columns_${i}.txt")

            SUBSET_ARG=""
            if [ -s "$SUBSET_FILE" ]; then
                SUBSET_ARG=$(sed -n "$((i + 1))p" "$SUBSET_FILE")
            fi

            bcftools annotate \
                -a "annotations_${i}.tsv.gz" \
                -h "header_${i}.txt" \
                -c "$COLUMN_SPEC" \
                $SUBSET_ARG \
                -Oz -o "temp_${i}.vcf.gz" \
                "$current_vcf"
            current_vcf="temp_${i}.vcf.gz"

            i=$((i + 1))
        done
        mv "$current_vcf" annotated_subset.vcf.gz
        tabix -p vcf annotated_subset.vcf.gz

        # Retrieve original variant IDs
        bcftools query -f '%INFO/original_ID\n' annotated_subset.vcf.gz \
            | awk 'NF && $1 != "."' \
            > original_ids_to_swap.txt

        # Remove original ID field
        bcftools annotate \
            -x INFO/original_ID \
            -Oz -o subset_for_concat.vcf.gz \
            annotated_subset.vcf.gz
        
        tabix -p vcf subset_for_concat.vcf.gz

        # Remove variants using original IDs
        if [ -s original_ids_to_swap.txt ]; then
            bcftools view \
                -e 'ID=@original_ids_to_swap.txt' \
                -Oz -o main_minus.vcf.gz \
                ~{vcf}
        else
            cp ~{vcf} main_minus.vcf.gz
        fi

        tabix -p vcf main_minus.vcf.gz

        # Concatenate retained original variants with subset variants 
        bcftools concat \
            --allow-overlaps \
            -Oz -o ~{prefix}.vcf.gz \
            main_minus.vcf.gz subset_for_concat.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * length(annotations_tsvs) * ceil(size(vcf, "GB") + size(annotations_tsvs, "GB")) + 10,
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
