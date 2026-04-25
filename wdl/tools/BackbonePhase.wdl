version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow BackbonePhase {
    input {
        File vcf
        File vcf_idx
        Array[File] base_vcfs
        Array[File] base_vcf_idxs
        String contig
        String prefix
        Boolean reassign_overlapped = true

        File? swap_samples_base

        String docker

        RuntimeAttr? runtime_attr_swap_sample_ids
        RuntimeAttr? runtime_attr_subset_base_vcf
        RuntimeAttr? runtime_attr_prepare_base_vcf
        RuntimeAttr? runtime_attr_assign_samples
        RuntimeAttr? runtime_attr_compute_flips
        RuntimeAttr? runtime_attr_apply_flips
    }

    scatter (i in range(length(base_vcfs))) {
        if (defined(swap_samples_base)) {
            call Helpers.SwapSampleIds {
                input:
                    vcf = base_vcfs[i],
                    vcf_idx = base_vcf_idxs[i],
                    sample_swap_list = select_first([swap_samples_base]),
                    prefix = "~{prefix}.base_~{i}.swapped",
                    docker = docker,
                    runtime_attr_override = runtime_attr_swap_sample_ids
            }
        }

        call Helpers.SubsetVcfToContig {
            input:
                vcf = select_first([SwapSampleIds.swapped_vcf, base_vcfs[i]]),
                vcf_idx = select_first([SwapSampleIds.swapped_vcf_idx, base_vcf_idxs[i]]),
                contig = contig,
                prefix = "~{prefix}.base_~{i}.contig",
                docker = docker,
                runtime_attr_override = runtime_attr_subset_base_vcf
        }

        call PrepareBaseVcf {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.base_~{i}.prepared",
                docker = docker,
                runtime_attr_override = runtime_attr_prepare_base_vcf
        }
    }

    call AssignSamplesToBaseVcfs {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            base_vcfs = PrepareBaseVcf.prepared_vcf,
            prefix = "~{prefix}.assignments",
            docker = docker,
            runtime_attr_override = runtime_attr_assign_samples
    }

    scatter (i in range(length(base_vcfs))) {
        call ComputePhaseFlips {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                base_vcf = PrepareBaseVcf.prepared_vcf[i],
                base_vcf_idx = PrepareBaseVcf.prepared_vcf_idx[i],
                assignment_tsv = AssignSamplesToBaseVcfs.assignment_tsv,
                base_vcf_index = i,
                reassign_overlapped = reassign_overlapped,
                prefix = "~{prefix}.base_~{i}.flips",
                docker = docker,
                runtime_attr_override = runtime_attr_compute_flips
        }
    }

    call ApplyPhaseFlips {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            flip_tsvs = ComputePhaseFlips.flips_tsv,
            missing_samples_file = AssignSamplesToBaseVcfs.missing_samples,
            prefix = "~{prefix}.transferred",
            docker = docker,
            runtime_attr_override = runtime_attr_apply_flips
    }

    output {
        File transferred_vcf = ApplyPhaseFlips.transferred_vcf
        File transferred_vcf_idx = ApplyPhaseFlips.transferred_vcf_idx
        File missing_samples = ApplyPhaseFlips.missing_samples
    }
}

task PrepareBaseVcf {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools norm -m-any ~{vcf} \
            | bcftools view -v snps \
            -Oz -o ~{prefix}.vcf.gz

        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File prepared_vcf = "~{prefix}.vcf.gz"
        File prepared_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
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

task AssignSamplesToBaseVcfs {
    input {
        File vcf
        File vcf_idx
        Array[File] base_vcfs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

with pysam.VariantFile("~{vcf}", "r") as vcf_in:
    vcf_sample_set = set(vcf_in.header.samples)

sample_to_base_idx = {}
base_vcf_paths = []
with open("~{write_lines(base_vcfs)}") as f:
    for line in f:
        path = line.strip()
        if path:
            base_vcf_paths.append(path)

for base_idx, path in enumerate(base_vcf_paths):
    with pysam.VariantFile(path, "r") as base_in:
        for sample in base_in.header.samples:
            if sample in vcf_sample_set and sample not in sample_to_base_idx:
                sample_to_base_idx[sample] = base_idx

with open("~{prefix}.tsv", "w") as out:
    for sample in sorted(sample_to_base_idx):
        out.write(f"{sample}\t{sample_to_base_idx[sample]}\n")

missing = sorted(vcf_sample_set - set(sample_to_base_idx))
with open("~{prefix}.missing_samples.txt", "w") as out:
    for sample in missing:
        out.write(sample + "\n")
CODE
    >>>

    output {
        File assignment_tsv = "~{prefix}.tsv"
        File missing_samples = "~{prefix}.missing_samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 + ceil(size(vcf, "GiB")) + ceil(size(base_vcfs, "GiB")),
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

task ComputePhaseFlips {
    input {
        File vcf
        File vcf_idx
        File base_vcf
        File base_vcf_idx
        File assignment_tsv
        Int base_vcf_index
        Boolean reassign_overlapped
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import os
import subprocess
from collections import defaultdict

import pysam


def normalize_allele(ref, alt):
    r, a = ref.upper(), alt.upper()
    while len(r) > 1 and len(a) > 1 and r[-1] == a[-1]:
        r = r[:-1]
        a = a[:-1]
    offset = 0
    while len(r) > 1 and len(a) > 1 and r[0] == a[0]:
        r = r[1:]
        a = a[1:]
        offset += 1
    return r, a, offset


def extract_sample(input_vcf, sample, output_vcf):
    subprocess.run(
        [
            "bcftools",
            "view",
            "-s",
            sample,
            "--min-ac",
            "1",
            "-Oz",
            "-o",
            output_vcf,
            input_vcf,
        ],
        check=True,
    )
    subprocess.run(["bcftools", "index", "-t", output_vcf], check=True)


def phase_gt_by_support(gt, support_hap):
    a, b = gt
    if support_hap not in (0, 1):
        return gt
    if a == 0 and b != 0:
        alt = b
        return (alt, 0) if support_hap == 0 else (0, alt)
    if a != 0 and b == 0:
        alt = a
        return (alt, 0) if support_hap == 0 else (0, alt)
    if a != 0 and b != 0:
        return (a, b) if support_hap == 0 else (b, a)
    return gt


def build_overlap_reassignments(sample_vcf, sample):
    unphased_rows = []
    phased_rows = []
    record_gts = {}

    with pysam.VariantFile(sample_vcf, "r") as sample_in:
        for rec in sample_in:
            gt_data = rec.samples[sample]
            gt = gt_data.get("GT")
            if gt is None or None in gt or len(gt) != 2:
                continue
            if gt[0] == gt[1]:
                continue

            record_gts[rec.id] = gt
            ps = gt_data.get("PS")

            for allele_order_idx in range(2):
                allele_idx = gt[allele_order_idx]
                if allele_idx == 0 or allele_idx > len(rec.alts):
                    continue

                norm_ref, norm_alt, offset = normalize_allele(rec.ref, rec.alts[allele_idx - 1])
                start0 = rec.pos + offset - 1
                end = start0 + len(norm_ref)
                span = end - start0

                if gt_data.phased:
                    phased_rows.append((
                        rec.contig,
                        start0,
                        end,
                        rec.id,
                        allele_order_idx,
                        ps if ps is not None else -1,
                        span,
                    ))
                else:
                    unphased_rows.append((
                        rec.contig,
                        start0,
                        end,
                        rec.id,
                        allele_order_idx,
                        span,
                    ))

    if not unphased_rows or not phased_rows:
        return {}

    unphased_rows.sort(key=lambda row: (row[0], row[1], row[2], row[3], row[4]))
    phased_rows.sort(key=lambda row: (row[0], row[1], row[2], row[3], row[4]))

    with open("unphased_overlap.bed", "w") as out:
        for row in unphased_rows:
            out.write("\t".join(map(str, row)) + "\n")

    with open("phased_overlap.bed", "w") as out:
        for row in phased_rows:
            out.write("\t".join(map(str, row)) + "\n")

    with open("overlap_matches.tsv", "w") as overlap_out:
        subprocess.run(
            [
                "bedtools",
                "intersect",
                "-wa",
                "-wb",
                "-f",
                "1.0",
                "-a",
                "unphased_overlap.bed",
                "-b",
                "phased_overlap.bed",
            ],
            check=True,
            stdout=overlap_out,
        )

    best_support = {}
    with open("overlap_matches.tsv") as matches_in:
        for line in matches_in:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 13:
                continue

            variant_id = parts[3]
            support_hap = int(parts[10])
            support_ps = int(parts[11])
            support_span = int(parts[12])

            current = best_support.get(variant_id)
            if current is None or support_span > current["span"]:
                best_support[variant_id] = {
                    "span": support_span,
                    "support_hap": support_hap,
                    "support_ps": support_ps,
                }

    reassigned = {}
    for variant_id, support in best_support.items():
        new_gt = phase_gt_by_support(record_gts[variant_id], support["support_hap"])
        if support["support_ps"] >= 0:
            reassigned[variant_id] = {
                "gt": new_gt,
                "ps": support["support_ps"],
            }

    for path in [
        "unphased_overlap.bed",
        "phased_overlap.bed",
        "overlap_matches.tsv",
    ]:
        if os.path.exists(path):
            os.remove(path)

    return reassigned


base_vcf_index = ~{base_vcf_index}
reassign_overlapped = "~{reassign_overlapped}" == "true"
assigned_samples = []
with open("~{assignment_tsv}") as f:
    for line in f:
        sample, idx = line.rstrip("\n").split("\t")
        if int(idx) == base_vcf_index:
            assigned_samples.append(sample)

with open("~{prefix}.tsv", "w") as out:
    out.write("VARIANT_ID\tSAMPLE\tNEW_GT\tNEW_PS\n")

    for sample_idx, sample in enumerate(assigned_samples):
        sample_vcf = f"sample_{sample_idx}.input.vcf.gz"
        sample_base_vcf = f"sample_{sample_idx}.base.vcf.gz"

        extract_sample("~{vcf}", sample, sample_vcf)
        extract_sample("~{base_vcf}", sample, sample_base_vcf)

        overlap_reassignments = build_overlap_reassignments(sample_vcf, sample) if reassign_overlapped else {}

        base_gts = {}
        with pysam.VariantFile(sample_base_vcf, "r") as base_in:
            for rec in base_in:
                if not rec.alts or len(rec.alts) != 1:
                    continue
                gt_data = rec.samples[sample]
                gt = gt_data.get("GT")
                if gt is None or None in gt or len(gt) != 2:
                    continue
                if gt[0] == gt[1] or not gt_data.phased:
                    continue
                key = (rec.contig, rec.pos, rec.ref.upper(), rec.alts[0].upper())
                base_gts[key] = (gt[0], gt[1])

        tally = defaultdict(lambda: [0, 0])
        matches = 0
        with pysam.VariantFile(sample_vcf, "r") as sample_in:
            for rec in sample_in:
                if not rec.alts:
                    continue
                gt_data = rec.samples[sample]
                gt = gt_data.get("GT")
                phased = gt_data.phased
                ps = gt_data.get("PS")
                if rec.id in overlap_reassignments:
                    gt = overlap_reassignments[rec.id]["gt"]
                    phased = True
                    ps = overlap_reassignments[rec.id]["ps"]
                if gt is None or None in gt or len(gt) != 2:
                    continue
                if gt[0] == gt[1] or not phased:
                    continue
                if ps is None:
                    continue
                for hap_idx in range(2):
                    allele_idx = gt[hap_idx]
                    if allele_idx == 0 or allele_idx > len(rec.alts):
                        continue
                    norm_ref, norm_alt, offset = normalize_allele(rec.ref, rec.alts[allele_idx - 1])
                    if len(norm_ref) != 1 or len(norm_alt) != 1:
                        continue
                    key = (rec.contig, rec.pos + offset, norm_ref, norm_alt)
                    if key not in base_gts:
                        continue
                    group = (1, 0) if hap_idx == 0 else (0, 1)
                    if base_gts[key] == group:
                        tally[ps][0] += 1
                    else:
                        tally[ps][1] += 1
                    matches += 1

        flip_set = set()
        for ps_group, (conc, disc) in tally.items():
            if disc > conc:
                flip_set.add(ps_group)

        with pysam.VariantFile(sample_vcf, "r") as sample_in:
            for rec in sample_in:
                gt_data = rec.samples[sample]
                orig_gt = gt_data.get("GT")
                gt = orig_gt
                phased = gt_data.phased
                ps = gt_data.get("PS")
                overlap_reassigned = rec.id in overlap_reassignments
                if overlap_reassigned:
                    gt = overlap_reassignments[rec.id]["gt"]
                    phased = True
                    ps = overlap_reassignments[rec.id]["ps"]
                if gt is None or None in gt or len(gt) != 2:
                    continue
                if gt[0] == gt[1] or not phased:
                    continue
                if ps is None:
                    continue
                a, b = gt[0], gt[1]
                should_flip = ps in flip_set

                final_gt = (b, a) if should_flip else (a, b)
                gt_changed = orig_gt is not None and tuple(orig_gt) != final_gt
                ps_changed = overlap_reassigned and gt_data.get("PS") != ps
                if overlap_reassigned or should_flip or gt_changed or ps_changed:
                    out.write(f"{rec.id}\t{sample}\t{final_gt[0]}|{final_gt[1]}\t{ps}\n")

        for path in [sample_vcf, sample_vcf + ".tbi", sample_base_vcf, sample_base_vcf + ".tbi"]:
            os.remove(path)
CODE
    >>>

    output {
        File flips_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 + ceil(size(vcf, "GiB")) + ceil(size(base_vcf, "GiB")) * 2,
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

task ApplyPhaseFlips {
    input {
        File vcf
        File vcf_idx
        Array[File] flip_tsvs
        File missing_samples_file
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from collections import defaultdict

import pysam

with open("~{missing_samples_file}") as src, open("~{prefix}.missing_samples.txt", "w") as dst:
    missing = [line.rstrip("\n") for line in src if line.strip()]
    for sample in missing:
        dst.write(sample + "\n")

flips = defaultdict(dict)
flip_tsv_paths = []
with open("~{write_lines(flip_tsvs)}") as f:
    for line in f:
        path = line.strip()
        if path:
            flip_tsv_paths.append(path)

for path in flip_tsv_paths:
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) not in (3, 4):
                continue
            variant_id, sample, new_gt = parts[:3]
            new_ps = parts[3] if len(parts) == 4 else None
            if sample in flips[variant_id]:
                continue
            a, b = new_gt.split("|")
            flips[variant_id][sample] = ((int(a), int(b)), None if new_ps in (None, "", ".") else int(new_ps))

total_flips = sum(len(sample_map) for sample_map in flips.values())
with pysam.VariantFile("~{vcf}", "r") as vcf_in:
    with pysam.VariantFile("transferred.vcf", "w", header=vcf_in.header) as vcf_out:
        for rec in vcf_in:
            if rec.id in flips:
                new_rec = rec.copy()
                for sample, (new_gt, new_ps) in flips[rec.id].items():
                    new_rec.samples[sample]["GT"] = new_gt
                    new_rec.samples[sample].phased = True
                    if new_ps is not None:
                        new_rec.samples[sample]["PS"] = new_ps
                vcf_out.write(new_rec)
            else:
                vcf_out.write(rec)
CODE

        bgzip transferred.vcf

        bcftools index -t transferred.vcf.gz

        mv transferred.vcf.gz ~{prefix}.vcf.gz

        mv transferred.vcf.gz.tbi ~{prefix}.vcf.gz.tbi
    >>>

    output {
        File transferred_vcf = "~{prefix}.vcf.gz"
        File transferred_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File missing_samples = "~{prefix}.missing_samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4 + ceil(size(vcf, "GiB") * 3),
        disk_gb: 10 + ceil(size(vcf, "GiB") * 3),
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
