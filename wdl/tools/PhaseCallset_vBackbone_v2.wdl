version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PhaseCallset_vBackbone {
    input {
        File vcf
        File vcf_idx
        Array[File] base_vcfs
        Array[File] base_vcf_idxs
        String prefix

        File? swap_samples_base

        String docker

        RuntimeAttr? runtime_attr_swap_sample_ids
        RuntimeAttr? runtime_attr_prepare_base_vcf
        RuntimeAttr? runtime_attr_compute_flips
        RuntimeAttr? runtime_attr_apply_flips
    }

    # Prepare each base_vcf: optional sample swap, then split multiallelics + filter to SNVs
    # Then compute phase flips for all samples assigned to this base_vcf
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

        call PrepareBaseVcf {
            input:
                vcf = select_first([SwapSampleIds.swapped_vcf, base_vcfs[i]]),
                vcf_idx = select_first([SwapSampleIds.swapped_vcf_idx, base_vcf_idxs[i]]),
                prefix = "~{prefix}.base_~{i}.prepared",
                docker = docker,
                runtime_attr_override = runtime_attr_prepare_base_vcf
        }

        call ComputePhaseFlips {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                base_vcf = PrepareBaseVcf.prepared_vcf,
                base_vcf_idx = PrepareBaseVcf.prepared_vcf_idx,
                base_vcf_index = i,
                prefix = "~{prefix}.base_~{i}.flips",
                docker = docker,
                runtime_attr_override = runtime_attr_compute_flips
        }
    }

    # Apply all flips in a single pass through the original VCF
    call ApplyPhaseFlips {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            flip_tsvs = ComputePhaseFlips.flips_tsv,
            missing_samples_files = ComputePhaseFlips.assigned_samples,
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

        # split multiallelics, then filter to biallelic SNVs only
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

task ComputePhaseFlips {
    input {
        File vcf
        File vcf_idx
        File base_vcf
        File base_vcf_idx
        Int base_vcf_index
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
from collections import defaultdict

def normalize_allele(ref, alt):
    """Trim shared suffix then prefix from a REF/ALT pair."""
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

base_vcf_index = ~{base_vcf_index}

# Find overlapping samples between vcf and this base_vcf
vcf_in = pysam.VariantFile("~{vcf}", "r")
vcf_samples = set(vcf_in.header.samples)
vcf_in.close()

base_in = pysam.VariantFile("~{base_vcf}", "r")
base_samples = set(base_in.header.samples)
base_in.close()

overlapping = sorted(vcf_samples & base_samples)
print(f"base_vcf_{base_vcf_index}: {len(overlapping)} overlapping samples")

with open("~{prefix}.samples.txt", "w") as f:
    for s in overlapping:
        f.write(s + "\n")

# Load base_vcf het phased SNP genotypes for all overlapping samples at once
base_gts = defaultdict(dict)
base_in = pysam.VariantFile("~{base_vcf}", "r")
for rec in base_in:
    if not rec.alts or len(rec.alts) != 1:
        continue
    key = (rec.contig, rec.pos, rec.ref.upper(), rec.alts[0].upper())
    for sample in overlapping:
        gt_data = rec.samples[sample]
        gt = gt_data.get("GT")
        if gt is None or None in gt or len(gt) != 2:
            continue
        if gt[0] == gt[1] or not gt_data.phased:
            continue
        base_gts[sample][key] = (gt[0], gt[1])
base_in.close()
for s in overlapping:
    print(f"  {s}: {len(base_gts.get(s, {}))} het phased SNVs in base_vcf")

# First pass through vcf: tally concordant/discordant per (sample, ps, haplotype_group)
tally = defaultdict(lambda: defaultdict(lambda: [0, 0]))
vcf_in = pysam.VariantFile("~{vcf}", "r")
for rec in vcf_in:
    if not rec.alts:
        continue
    for sample in overlapping:
        gt_data = rec.samples[sample]
        gt = gt_data.get("GT")
        if gt is None or None in gt or len(gt) != 2:
            continue
        if gt[0] == gt[1] or not gt_data.phased:
            continue
        ps = gt_data.get("PS")
        if ps is None:
            continue
        sample_base = base_gts.get(sample)
        if not sample_base:
            continue
        for hap_idx in range(2):
            allele_idx = gt[hap_idx]
            if allele_idx == 0 or allele_idx > len(rec.alts):
                continue
            norm_ref, norm_alt, offset = normalize_allele(rec.ref, rec.alts[allele_idx - 1])
            if len(norm_ref) != 1 or len(norm_alt) != 1:
                continue
            key = (rec.contig, rec.pos + offset, norm_ref, norm_alt)
            if key not in sample_base:
                continue
            group = (1, 0) if hap_idx == 0 else (0, 1)
            base_gt = sample_base[key]
            if base_gt == group:
                tally[sample][(ps, group)][0] += 1
            else:
                tally[sample][(ps, group)][1] += 1
vcf_in.close()

# Determine which (sample, ps, group) combinations should be flipped
flip_sets = {}
for sample in overlapping:
    sample_flips = set()
    for (ps, group), (conc, disc) in tally[sample].items():
        if disc > conc:
            sample_flips.add((ps, group))
    flip_sets[sample] = sample_flips
    matched = sum(c + d for c, d in tally[sample].values())
    flipped = len(sample_flips)
    total = len(tally[sample])
    print(f"  {sample}: {matched} matches, flipping {flipped}/{total} haplotype groups")

# Second pass through vcf: write out flips
with open("~{prefix}.tsv", "w") as out:
    out.write("VARIANT_ID\tSAMPLE\tNEW_GT\n")
    vcf_in = pysam.VariantFile("~{vcf}", "r")
    for rec in vcf_in:
        for sample in overlapping:
            flip_set = flip_sets.get(sample)
            if not flip_set:
                continue
            gt_data = rec.samples[sample]
            gt = gt_data.get("GT")
            if gt is None or None in gt or len(gt) != 2:
                continue
            if not gt_data.phased or gt[0] == gt[1]:
                continue
            ps = gt_data.get("PS")
            if ps is None:
                continue
            a, b = gt[0], gt[1]
            should_flip = False
            if a != 0 and b == 0:
                should_flip = (ps, (1, 0)) in flip_set
            elif a == 0 and b != 0:
                should_flip = (ps, (0, 1)) in flip_set
            else:
                should_flip = (ps, (1, 0)) in flip_set and (ps, (0, 1)) in flip_set
            if should_flip:
                out.write(f"{rec.id}\t{sample}\t{b}|{a}\n")
    vcf_in.close()
CODE
    >>>

    output {
        File flips_tsv = "~{prefix}.tsv"
        File assigned_samples = "~{prefix}.samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4 + ceil(size(base_vcf, "GiB") * 2),
        disk_gb: ceil(size(vcf, "GiB")) + ceil(size(base_vcf, "GiB")) + 10,
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
        Array[File] missing_samples_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
from collections import defaultdict

# Determine missing samples: vcf samples not found in any base_vcf
vcf_in = pysam.VariantFile("~{vcf}", "r")
vcf_samples = set(vcf_in.header.samples)
vcf_in.close()

assigned_samples = set()
samples_files = []
with open("~{write_lines(missing_samples_files)}") as f:
    for line in f:
        line = line.strip()
        if line:
            samples_files.append(line)

for path in samples_files:
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s:
                assigned_samples.add(s)

missing = sorted(vcf_samples - assigned_samples)
with open("~{prefix}.missing_samples.txt", "w") as f:
    for s in missing:
        f.write(s + "\n")
print(f"{len(assigned_samples)} assigned samples, {len(missing)} missing")

# Load all flip TSVs: {variant_id: {sample: (a, b)}}
flips = defaultdict(dict)
flip_tsv_paths = []
with open("~{write_lines(flip_tsvs)}") as f:
    for line in f:
        line = line.strip()
        if line:
            flip_tsv_paths.append(line)

for path in flip_tsv_paths:
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                continue
            vid, sample, new_gt = parts
            # First-match: skip if this (variant, sample) was already assigned by an earlier base_vcf
            if sample in flips.get(vid, {}):
                continue
            alleles = new_gt.split("|")
            flips[vid][sample] = (int(alleles[0]), int(alleles[1]))

total_flips = sum(len(v) for v in flips.values())
print(f"Loaded {total_flips} flips across {len(flips)} variants")

# Single pass: apply flips
vcf_in = pysam.VariantFile("~{vcf}", "r")
vcf_out = pysam.VariantFile("transferred.vcf", "w", header=vcf_in.header)
applied = 0
for rec in vcf_in:
    if rec.id in flips:
        new_rec = rec.copy()
        sample_flips = flips[rec.id]
        for sample, new_gt in sample_flips.items():
            new_rec.samples[sample]["GT"] = new_gt
            new_rec.samples[sample].phased = True
            applied += 1
        vcf_out.write(new_rec)
    else:
        vcf_out.write(rec)
vcf_in.close()
vcf_out.close()
print(f"Applied {applied} genotype flips")
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
        preemptible_tries: 1,
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
