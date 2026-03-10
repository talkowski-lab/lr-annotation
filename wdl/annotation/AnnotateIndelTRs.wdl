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

        Boolean dont_run_trf = false
        Boolean dont_allow_interruptions = false
        Int min_tandem_repeat_length = 9
        Int min_repeats = 3
        Int min_repeat_unit_length = 1

        String stranalysis_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_filter
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

        call RunFilterVcfToTRs {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{contig}.filtered_trs",
                dont_run_trf = dont_run_trf,
                dont_allow_interruptions = dont_allow_interruptions,
                min_tandem_repeat_length = min_tandem_repeat_length,
                min_repeats = min_repeats,
                min_repeat_unit_length = min_repeat_unit_length,
                docker = stranalysis_docker,
                runtime_attr_override = runtime_attr_filter
        }

        call CreateTRAnnotationTsv {
            input:
                vcf = RunFilterVcfToTRs.filtered_vcf,
                vcf_idx = RunFilterVcfToTRs.filtered_vcf_idx,
                tr_tsv = RunFilterVcfToTRs.tr_tsv,
                prefix = "~{prefix}.~{contig}.tr_annotations",
                docker = stranalysis_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatTsvs as MergeAnnotations {
        input:
            tsvs = CreateTRAnnotationTsv.tr_annotations_tsv,
            sort_output = false,
            preserve_header = false,
            compressed_tsvs = false,
            prefix = "~{prefix}.tr_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File tr_annotations_tsv = MergeAnnotations.concatenated_tsv
    }
}

task RunFilterVcfToTRs {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        String prefix

        Boolean dont_run_trf
        Boolean dont_allow_interruptions
        Int min_tandem_repeat_length
        Int min_repeats
        Int min_repeat_unit_length

        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        bcftools view \
            -e 'INFO/allele_type="trv"' \
            -Oz -o ~{prefix}.filtered.vcf.gz \
            ~{vcf}
        
        tabix ~{prefix}.filtered.vcf.gz

        if [ "~{dont_run_trf}" == "true" ]; then
            TRF_ARGS="--dont-run-trf"
        else
            TRF_ARGS="--trf-executable-path $(which trf) --trf-threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}"
        fi

        python3 -m str_analysis.filter_vcf_to_tandem_repeats catalog \
            -R ~{ref_fa} \
            --output-prefix ~{prefix} \
            --min-tandem-repeat-length ~{min_tandem_repeat_length} \
            --min-repeats ~{min_repeats} \
            --min-repeat-unit-length ~{min_repeat_unit_length} \
            --write-tsv \
            ~{if dont_allow_interruptions then "--dont-allow-interruptions" else ""} \
            $TRF_ARGS \
            ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File tr_tsv = "~{prefix}.tandem_repeats.tsv.gz"
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

task CreateTRAnnotationTsv {
    input {
        File vcf
        File vcf_idx
        File tr_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        python3 <<CODE
import gzip
import pysam

vcf_file = pysam.VariantFile("~{vcf}")

annotation_cols = [
    "INS_or_DEL", "Motif", "CanonicalMotif", "MotifSize", "NumRepeatsInReference",
    "Locus", "LocusId", "SummaryString", "IsFoundInReference", "IsPureRepeat", "DetectionMode",
]

with gzip.open("~{tr_tsv}", "rt") as tsv_in, open("~{prefix}.tsv", "w") as out:
    col_idx = None
    for line in tsv_in:
        if col_idx is None:
            col_idx = {col: i for i, col in enumerate(line.strip().split("\t"))}
            continue
        fields = line.strip().split("\t")
        chrom = fields[col_idx["Chrom"]]
        vcf_pos = int(fields[col_idx["VcfPos"]])
        ins_or_del = fields[col_idx["INS_or_DEL"]]
        try:
            records = list(vcf_file.fetch(chrom, vcf_pos - 1, vcf_pos))
        except ValueError:
            continue
        for record in records:
            if record.pos != vcf_pos:
                continue
            allele_type = record.info.get("allele_type", "")
            if allele_type.upper() != ins_or_del:
                continue
            vid = record.id if record.id else "."
            anno_vals = [fields[col_idx[col]] for col in annotation_cols]
            out.write("\t".join([chrom, str(vcf_pos), record.ref, record.alts[0], vid] + anno_vals) + "\n")
CODE
    >>>

    output {
        File tr_annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(tr_tsv, "GB")) + 10,
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
