version 1.0

import "general/Structs.wdl"
import "general/Helpers.wdl"

workflow Annotate {
    input {
        File vcf
        File vcf_idx
        File annotations_tsv
        Array[String] contigs

        Array[String] info_names
        Array[String] info_descriptions
        Array[String] info_types
        Array[String] info_numbers

        String prefix
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
                vcf_index = vcf_idx,
                contig = contig,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.SubsetTsvToContig {
            input:
                tsv = annotations_tsv,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_tsv
        }

        call AnnotateVcfWithTsv {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_index,
                annotations_tsv = SubsetTsvToContig.subset_tsv,
                info_names = info_names,
                info_descriptions = info_descriptions,
                info_types = info_types,
                info_numbers = info_numbers,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateVcfWithTsv.annotated_vcf,
            vcfs_idx = AnnotateVcfWithTsv.annotated_vcf_idx,
            outfile_prefix = prefix + ".annotated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File annotated_vcf = ConcatVcfs.concat_vcf
        File annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AnnotateVcfWithTsv {
    input {
        File vcf
        File vcf_idx
        File annotations_tsv
        Array[String] info_names
        Array[String] info_descriptions
        Array[String] info_types
        Array[String] info_numbers
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import sys

info_names = "~{sep=',' info_names}".split(',')
info_descriptions = "~{sep=',' info_descriptions}".split(',')
info_types = "~{sep=',' info_types}".split(',')
info_numbers = "~{sep=',' info_numbers}".split(',')

if len(info_names) != len(info_descriptions) or len(info_names) != len(info_types) or len(info_names) != len(info_numbers):
    sys.stderr.write("Error: info_names, info_descriptions, info_types, and info_numbers must have the same length\n")
    sys.exit(1)

with open("header.txt", "w") as f:
    for name, desc, type_val, number in zip(info_names, info_descriptions, info_types, info_numbers):
        f.write(f'##INFO=<ID={name},Number={number},Type={type_val},Description="{desc}">\n')

columns = ['CHROM', 'POS', 'REF', 'ALT', 'ID'] + info_names
column_spec = ','.join(['CHROM', 'POS', 'REF', 'ALT', 'ID'] + [f'INFO/{name}' for name in info_names])

with open("columns.txt", "w") as f:
    f.write(column_spec)

EOF

        COLUMN_SPEC=$(cat columns.txt)

        if [ -s ~{annotations_tsv} ]; then
            bgzip -c ~{annotations_tsv} > annotations.tsv.gz
            tabix -s1 -b2 -e2 annotations.tsv.gz

            bcftools annotate \
                -a annotations.tsv.gz \
                -h header.txt \
                -c "$COLUMN_SPEC" \
                -Oz -o ~{prefix}.annotated.vcf.gz \
                ~{vcf}
        else
            cp ~{vcf} ~{prefix}.annotated.vcf.gz
        fi

        tabix -p vcf ~{prefix}.annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2*ceil(size(vcf, "GB") + size(annotations_tsv, "GB")) + 10,
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
