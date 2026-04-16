version 1.0

import "../utils/Structs.wdl"

workflow FillFormatFields {
    input {
        File unfilled_vcf
        File unfilled_vcf_idx
        File filled_vcf
        File filled_vcf_idx
        String prefix

        Array[String] format_fields
        String? include_field
        String? include_value
        Boolean fill_alt_gts = false
        Boolean fill_ref_gts = false

        String utils_docker

        RuntimeAttr? runtime_attr_fill
    }

    call FillVcfFormatFields {
        input:
            unfilled_vcf = unfilled_vcf,
            unfilled_vcf_idx = unfilled_vcf_idx,
            filled_vcf = filled_vcf,
            filled_vcf_idx = filled_vcf_idx,
            format_fields = format_fields,
            include_field = include_field,
            include_value = include_value,
            fill_alt_gts = fill_alt_gts,
            fill_ref_gts = fill_ref_gts,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_fill
    }

    output {
        File refilled_vcf = FillVcfFormatFields.output_vcf
        File refilled_vcf_idx = FillVcfFormatFields.output_vcf_idx
    }
}

task FillVcfFormatFields {
    input {
        File unfilled_vcf
        File unfilled_vcf_idx
        File filled_vcf
        File filled_vcf_idx
        Array[String] format_fields
        String? include_field
        String? include_value
        Boolean fill_alt_gts = false
        Boolean fill_ref_gts = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -sf "~{filled_vcf_idx}" "~{filled_vcf}.tbi"

        python3 <<CODE
import pysam

with open("~{write_lines(format_fields)}") as fh:
    format_fields = [l.strip() for l in fh if l.strip()]

assert "GT" not in format_fields, "GT must not be passed in format_fields; use fill_alt_gts/fill_ref_gts instead"

include_field = "~{default="" include_field}" or None
include_value = "~{default="" include_value}" or None
fill_alt_gts = ~{true="True" false="False" fill_alt_gts}
fill_ref_gts = ~{true="True" false="False" fill_ref_gts}

unfilled_in = pysam.VariantFile("~{unfilled_vcf}")
filled_in = pysam.VariantFile("~{filled_vcf}")

out_header = unfilled_in.header.copy()
for field in format_fields:
    if field in filled_in.header.formats and field not in out_header.formats:
        out_header.add_record(filled_in.header.formats[field].record)

filled_samples = set(filled_in.header.samples)
common_samples = [s for s in out_header.samples if s in filled_samples]

out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=out_header)

def variant_passes_include(rec):
    if include_field is None:
        return True
    val = rec.info.get(include_field)
    if isinstance(val, tuple):
        val = val[0] if val else None
    return str(val) == include_value

def gt_is_alt(gt):
    return any(a is not None and a > 0 for a in gt)

for unfilled_rec in unfilled_in:
    unfilled_rec.translate(out_header)
    match = None
    for cand in filled_in.fetch(unfilled_rec.chrom, unfilled_rec.start, unfilled_rec.stop):
        if cand.id == unfilled_rec.id:
            match = cand
            break

    if match and variant_passes_include(unfilled_rec):
        for sample in common_samples:
            for field in format_fields:
                if field in match.format:
                    try:
                        unfilled_rec.samples[sample][field] = match.samples[sample][field]
                    except Exception as e:
                        unfilled_val = unfilled_rec.samples[sample][field] if field in unfilled_rec.samples[sample] else "N/A"
                        filled_val = match.samples[sample][field] if field in match.samples[sample] else "N/A"
                        print(f"[{unfilled_rec.id}] Couldn't set {field} for {sample}: {unfilled_val} --> {filled_val}")
                        pass

            if fill_alt_gts or fill_ref_gts:
                src_gt = match.samples[sample].get("GT")
                cur_gt = unfilled_rec.samples[sample].get("GT")
                if src_gt is not None and cur_gt is not None:
                    cur_is_alt = gt_is_alt(cur_gt)
                    if (fill_alt_gts and cur_is_alt) or (fill_ref_gts and not cur_is_alt):
                        unfilled_rec.samples[sample]["GT"] = src_gt
                        unfilled_rec.samples[sample].phased = match.samples[sample].phased

    out.write(unfilled_rec)

out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(unfilled_vcf, "GB") + size(filled_vcf, "GB")) + 10,
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