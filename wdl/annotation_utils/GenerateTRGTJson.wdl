version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow GenerateTRGTJson {
    input {
        File lps_tsv
        Array[String] contigs
        String prefix

        String stranalysis_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_convert
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call SubsetLpsTsvToContig {
            input:
                tsv = lps_tsv,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_tsv
        }

        call ConvertLPSTableToAFHistograms {
            input:
                lps_tsv = SubsetLpsTsvToContig.subset_tsv,
                prefix = "~{prefix}.~{contig}.af_histograms",
                docker = stranalysis_docker,
                runtime_attr_override = runtime_attr_convert
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = ConvertLPSTableToAFHistograms.af_histograms_tsv,
            sort_output = false,
            preserve_header = true,
            compressed_tsvs = true,
            prefix = "~{prefix}.af_histograms",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File af_histograms_tsv = ConcatTsvs.concatenated_tsv
    }
}

task SubsetLpsTsvToContig {
    input {
        File tsv
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        zcat -f ~{tsv} \
            | awk -F'\t' -v contig="~{contig}" '
                BEGIN {
                    if (contig ~ /^chr/) contig = substr(contig,4)
                }

                {
                    first = $1
                    split(first, arr, ",")
                    split(arr[1], parts, "-")
                    chr = parts[1]
                    if (chr ~ /^chr/) chr = substr(chr,4)
                    if (chr == contig) print
                }
            ' \
            | gzip -c \
            > ~{prefix}.tsv.gz
    >>>

    output {
        File subset_tsv = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 5,
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

task ConvertLPSTableToAFHistograms {
    input {
        File lps_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        python3 -m str_analysis.convert_multisample_LPS_table_to_allele_frequency_histograms \
            --input-table ~{lps_tsv} \
            --no-header

        mv "$(dirname ~{lps_tsv})"/*.per_locus_and_motif.*.tsv.gz ~{prefix}.tsv.gz
    >>>

    output {
        File af_histograms_tsv = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(lps_tsv, "GB")) + 20,
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
