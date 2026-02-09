version 1.0

import "../utils/Helpers.wdl"
import "../tools/RepeatMasker.wdl"
import "../utils/Structs.wdl"

workflow AnnotateL1MEAID {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int min_length

        File? fa_override

        String utils_docker
        String repeatmasker_docker
        String l1meaid_docker
        String intact_mei_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_ins_to_fa
        RuntimeAttr? runtime_attr_repeat_masker
        RuntimeAttr? runtime_attr_limeaid
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_contigs
        RuntimeAttr? runtime_attr_concat_intactmei
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call RepeatMasker.RepeatMasker {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                min_length = min_length,
                prefix = "~{prefix}.~{contig}.rm",
                fa_override = fa_override,
                utils_docker = utils_docker,
                repeatmasker_docker = repeatmasker_docker,
                runtime_attr_ins_to_fa = runtime_attr_ins_to_fa,
                runtime_attr_repeat_masker = runtime_attr_repeat_masker
        }

        call L1MEAID {
            input:
                rm_fa = select_first([fa_override, RepeatMasker.rm_fa]),
                rm_out = RepeatMasker.rm_out,
                prefix = "~{prefix}.~{contig}.l1meaid",
                docker = l1meaid_docker,
                runtime_attr_override = runtime_attr_limeaid
        }

        call IntactMEI {
            input:
                l1meaid_output = L1MEAID.l1meaid_output,
                prefix = "~{prefix}.~{contig}.intactmei",
                docker = intact_mei_docker,
                runtime_attr_override = runtime_attr_filter
        }

        call GenerateAnnotationTable {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                filtered_tsv = IntactMEI.filtered_output,
                prefix = "~{prefix}.~{contig}.intactmei_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate
        }
    }

    call Helpers.ConcatTsvs as MergeIntactMeiOutput {
        input:
            tsvs = IntactMEI.filtered_output,
            skip_sort = true,
            prefix = prefix + ".intactmei_output",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_intactmei
    }

    call Helpers.ConcatTsvs as MergeAnnotations {
        input:
            tsvs = GenerateAnnotationTable.annotations_tsv,
            prefix = prefix + ".intactmei_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_contigs
    }

    output {
        File annotations_tsv_l1meaid = MergeAnnotations.concatenated_tsv
        File annotations_full_l1meaid = MergeIntactMeiOutput.concatenated_tsv
    }
}

task L1MEAID {
    input {
        File rm_fa
        File rm_out
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/src/L1ME-AID/limeaid.py \
            -i ~{rm_fa} \
            -r ~{rm_out} \
            -o ~{prefix}.txt
    >>>

    output {
        File l1meaid_output = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(rm_fa, "GB")) + 5,
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

task IntactMEI {
    input {
        File l1meaid_output
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        perl /opt/src/utility/limeaid.filter.pl \
            ~{l1meaid_output} \
            > ~{prefix}.tsv
    >>>

    output {
        File filtered_output = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(l1meaid_output, "GB")) + 5,
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

task GenerateAnnotationTable {
    input {
        File vcf
        File filtered_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ~{vcf} > vcf_lookup.tsv

        python3 <<EOF
import sys

vcf_lookup_file = "vcf_lookup.tsv"
input_tsv = "~{filtered_tsv}"
output_anno = "~{prefix}.tsv"

vcf_lookup = {}
with open(vcf_lookup_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) >= 5:
            chrom, pos, ref, alt, var_id = fields[0], fields[1], fields[2], fields[3], fields[4]
            vcf_lookup[(chrom, pos, alt)] = (ref, var_id)

with open(input_tsv, 'r') as f_in, open(output_anno, 'w') as f_out:
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) < 12: continue
        
        var_id_str = parts[0]
        sequence = parts[1]
        
        classification = parts[8]
        structure = parts[10]
        
        me_type = None
        
        if classification == "SINE/Alu" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "ALU"
        elif classification == "Retroposon/SVA" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "SVA"
        elif classification == "LINE/L1" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "LINE"
        
        if me_type:
            id_parts = var_id_str.split(';')
            if len(id_parts) >= 1:
                location_part = id_parts[0]
                chrom, pos = location_part.rsplit(':', 1)
                key = (chrom, pos, sequence)
                if key in vcf_lookup:
                    ref, var_id = vcf_lookup[key]
                    f_out.write(f"{chrom}\t{pos}\t{ref}\t{sequence}\t{var_id}\t{me_type}\n")
                else:
                    sys.stderr.write(f"Warning: No matching VCF record for {chrom}:{pos}...\n")
EOF
    >>>

    output {
        File annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(filtered_tsv, "GB")) + 5,
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
