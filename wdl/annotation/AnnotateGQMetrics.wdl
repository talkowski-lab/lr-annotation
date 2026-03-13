version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateGQMetrics {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Array[String] gq_fields
        Array[Array[Int]] gq_bins
        Array[String] gq_variant_filters
        Array[Boolean] gq_larger_field

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_generate_tsv
        RuntimeAttr? runtime_attr_apply_annotations
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = prefix + "." + contig,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        scatter (i in range(length(gq_fields))) {
            call GenerateGQAnnotationTsv {
                input:
                    vcf = SubsetVcfToContig.subset_vcf,
                    vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                    gq_field = gq_fields[i],
                    gq_bins = gq_bins[i],
                    gq_variant_filter = gq_variant_filters[i],
                    gq_larger_field = gq_larger_field[i],
                    prefix = prefix + "." + contig + "." + gq_fields[i],
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_generate_tsv
            }
        }

        call ApplyGQAnnotations {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                annotation_tsvs = GenerateGQAnnotationTsv.annotation_tsv,
                annotation_tsv_idxs = GenerateGQAnnotationTsv.annotation_tsv_idx,
                header_files = GenerateGQAnnotationTsv.header_file,
                col_spec_files = GenerateGQAnnotationTsv.col_spec_file,
                prefix = prefix + "." + contig + ".gq_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_apply_annotations
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = ApplyGQAnnotations.annotated_vcf,
            vcf_idxs = ApplyGQAnnotations.annotated_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = prefix + ".gq_annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File gq_annotated_vcf = ConcatVcfs.concat_vcf
        File gq_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task GenerateGQAnnotationTsv {
    input {
        File vcf
        File vcf_idx
        String gq_field
        Array[Int] gq_bins
        String gq_variant_filter
        Boolean gq_larger_field
        String prefix

        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        FILTER="~{gq_variant_filter}"
        FILTER_ARGS=""
        if [ -n "$FILTER" ] && [ "$FILTER" != "." ] && [ "$FILTER" != "None" ]; then
            echo "$FILTER" > filter.txt
            FILTER_ARGS="-i @filter.txt"
        fi

        bcftools view $FILTER_ARGS -Ov "~{vcf}" | python3 - <<'EOF'
import json
import pysam

field = "~{gq_field}"
with open("~{write_json(gq_bins)}") as _f:
    bins = json.load(_f)
larger_flag = ~{true="True" false="False" gq_larger_field}
prefix = "~{prefix}"

# Build and write header lines
bin_edges_str = "0|" + "|".join(map(str, bins))
col_names = [f"{field}_hist_all_bin_freq", f"{field}_hist_alt_bin_freq"]
header_lines = [
    f'##INFO=<ID={field}_hist_all_bin_freq,Number=1,Type=String,Description="Histogram for {field} in all individuals; bin edges are: {bin_edges_str}">',
    f'##INFO=<ID={field}_hist_alt_bin_freq,Number=1,Type=String,Description="Histogram for {field} in heterozygous/alt individuals; bin edges are: {bin_edges_str}">',
]
if larger_flag:
    col_names += [f"{field}_hist_all_n_larger", f"{field}_hist_alt_n_larger"]
    header_lines += [
        f'##INFO=<ID={field}_hist_all_n_larger,Number=1,Type=Integer,Description="Count of all genotypes with {field} > {bins[-1]}">',
        f'##INFO=<ID={field}_hist_alt_n_larger,Number=1,Type=Integer,Description="Count of alt genotypes with {field} > {bins[-1]}">',
    ]

with open(f"{prefix}.header.txt", "w") as f:
    f.write("\n".join(header_lines) + "\n")

# Column spec for bcftools annotate: match by CHROM/POS/REF/ALT, update ID, then INFO fields
col_spec = "CHROM,POS,REF,ALT,~ID," + ",".join(f"INFO/{c}" for c in col_names)
with open(f"{prefix}.col_spec.txt", "w") as f:
    f.write(col_spec)

def get_bin_index(val, bins):
    for j, b in enumerate(bins):
        if val <= b:
            return j
    return len(bins)

vcf_in = pysam.VariantFile("-", "r")
with open(f"{prefix}.tsv", "w") as out:
    for record in vcf_in:
        counts_all = [0] * (len(bins) + 1)
        counts_alt = [0] * (len(bins) + 1)

        for sample_data in record.samples.values():
            val = sample_data.get(field)
            if val is None:
                continue
            if isinstance(val, tuple):
                if not val or val[0] is None:
                    continue
                val = val[0]
            val = float(val)

            gt = sample_data.get("GT")
            is_alt = gt is not None and any(a is not None and a > 0 for a in gt)

            b_idx = get_bin_index(val, bins)
            counts_all[b_idx] += 1
            if is_alt:
                counts_alt[b_idx] += 1

        vid = record.id if record.id else "."
        row = [
            record.chrom, str(record.pos), record.ref,
            ",".join(record.alts) if record.alts else ".",
            vid,
            "|".join(map(str, counts_all)),
            "|".join(map(str, counts_alt)),
        ]
        if larger_flag:
            row += [str(counts_all[-1]), str(counts_alt[-1])]

        out.write("\t".join(row) + "\n")
EOF

        bgzip -c "~{prefix}.tsv" > "~{prefix}.tsv.gz"
        tabix -s1 -b2 -e2 "~{prefix}.tsv.gz"
    >>>

    output {
        File annotation_tsv = "~{prefix}.tsv.gz"
        File annotation_tsv_idx = "~{prefix}.tsv.gz.tbi"
        File header_file = "~{prefix}.header.txt"
        File col_spec_file = "~{prefix}.col_spec.txt"
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

task ApplyGQAnnotations {
    input {
        File vcf
        File vcf_idx
        Array[File] annotation_tsvs
        Array[File] annotation_tsv_idxs
        Array[File] header_files
        Array[File] col_spec_files
        String prefix

        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        annotation_tsvs=(~{sep=' ' annotation_tsvs})
        header_files=(~{sep=' ' header_files})
        col_spec_files=(~{sep=' ' col_spec_files})

        current_vcf="~{vcf}"
        for i in "${!annotation_tsvs[@]}"; do
            col_spec=$(cat "${col_spec_files[$i]}")

            bcftools annotate \
                -a "${annotation_tsvs[$i]}" \
                -h "${header_files[$i]}" \
                -c "$col_spec" \
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
