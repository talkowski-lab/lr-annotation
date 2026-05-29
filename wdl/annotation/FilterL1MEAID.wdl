version 1.0

workflow FilterL1MEAID {
    input {
        File vcf
        File raw_intact_mei_tsv
        String prefix

        String utils_docker
        RuntimeAttr? runtime_attr_annotate
    }

    # This workflow is for datasets where RepeatMasker, L1ME-AID, and the
    # IntactMEI filtering step have already been run. Re-run only the current
    # `GenerateAnnotationTable` logic against the existing IntactMEI TSV.
    call GenerateAnnotationTable {
        input:
            vcf = vcf,
            filtered_tsv = raw_intact_mei_tsv,
            prefix = "~{prefix}.intactmei_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_annotate
    }

    output {
        File annotations_tsv_l1meaid = GenerateAnnotationTable.annotations_tsv
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
import ast


def get_mei_len(element_hits, subfam, classification):
    try:
        hits = ast.literal_eval(element_hits)
    except (SyntaxError, ValueError):
        return 0

    intervals = []
    id_to_intervals = {}
    for hit in hits:
        fields = str(hit).split()
        if len(fields) < 14:
            continue

        anno = fields[9]
        anno_class = fields[10]
        if classification == "Retroposon/SVA":
            keep = (anno_class == classification)
        else:
            keep = (anno == subfam)
        if not keep:
            continue

        try:
            start = int(fields[5])
            end = int(fields[6])
        except ValueError:
            continue
        if start > end:
            start, end = end, start
        if classification == "Retroposon/SVA":
            intervals.append((start, end))
        else:
            hit_id = fields[-1]
            id_to_intervals.setdefault(hit_id, []).append((start, end))

    if classification != "Retroposon/SVA":
        best_len = 0
        for hit_id in id_to_intervals:
            curr_intervals = sorted(id_to_intervals[hit_id])
            merged = [curr_intervals[0]]
            for start, end in curr_intervals[1:]:
                prev_start, prev_end = merged[-1]
                if start <= prev_end + 1:
                    merged[-1] = (prev_start, max(prev_end, end))
                else:
                    merged.append((start, end))
            curr_len = sum(end - start + 1 for start, end in merged)
            if curr_len > best_len:
                best_len = curr_len
        return best_len

    if not intervals:
        return 0

    intervals.sort()
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))

    return sum(end - start + 1 for start, end in merged)

vcf_lookup_file = "vcf_lookup.tsv"
input_tsv = "~{filtered_tsv}"
output_anno = "~{prefix}.tsv"

vcf_lookup = {}
with open(vcf_lookup_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) >= 5:
            chrom, pos, ref, alt, var_id = fields[0], fields[1], fields[2], fields[3], fields[4]
            vcf_lookup[(chrom, pos, ref, alt)] = var_id

with open(input_tsv, 'r') as f_in, open(output_anno, 'w') as f_out:
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) < 12: 
            continue
        
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
            full_id = parts[0]
            full_id_parts = full_id.split(';')
            chrom = full_id_parts[0].split(':')[0]
            pos = full_id_parts[0].split(':')[1]
            ref = full_id_parts[1].split('_')[0]
            alt = parts[1]
            subfam = parts[4]
            ori = parts[9]
            mei_len = get_mei_len(parts[2], subfam, classification)
            query_len = int(parts[3])
            mei_fraction = float(mei_len) / query_len if query_len > 0 else 0.0
            if me_type == "ALU":
                if mei_fraction < 0.7:
                    continue
            elif mei_fraction < 0.2:
                continue
            key = (chrom, pos, ref, alt)
            if key in vcf_lookup:
                f_out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{vcf_lookup[key]}\t{me_type}\t{subfam}\t{mei_len}\t{ori}\t{structure}\n")
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}