version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateIndelTRs {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        Int min_tandem_repeat_length = 9
        Int min_repeats = 3
        Int min_repeat_unit_length = 1

        String stranalysis_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_filter
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

        call RunFilterVcfToTRs {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}.tr_annotations",
                min_tandem_repeat_length = min_tandem_repeat_length,
                min_repeats = min_repeats,
                min_repeat_unit_length = min_repeat_unit_length,
                docker = stranalysis_docker,
                runtime_attr_override = runtime_attr_filter
        }
    }

    call Helpers.ConcatTsvs as MergeAnnotations {
        input:
            tsvs = RunFilterVcfToTRs.tr_annotations_tsv,
            sort_output = false,
            prefix = "~{prefix}.tr_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File annotations_tsv_trs = MergeAnnotations.concatenated_tsv
    }
}

task RunFilterVcfToTRs {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        Int min_tandem_repeat_length
        Int min_repeats
        Int min_repeat_unit_length
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        bcftools view \
            -e 'INFO/allele_type="trv"' \
            -Oz -o filtered.vcf.gz \
            ~{vcf}
        
        tabix filtered.vcf.gz

        python3 -m str_analysis.filter_vcf_to_tandem_repeats catalog \
            -R ~{ref_fa} \
            --output-prefix ~{prefix} \
            --min-tandem-repeat-length ~{min_tandem_repeat_length} \
            --min-repeats ~{min_repeats} \
            --min-repeat-unit-length ~{min_repeat_unit_length} \
            --write-vcf \
            --trf-executable-path $(which trf) \
            --trf-threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            filtered.vcf.gz

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t1\n' \
            ~{prefix}.tandem_repeats.vcf.gz \
            > ~{prefix}.tsv
    >>>

    output {
        File tr_annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(ref_fa, "GB")) + 20,
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
