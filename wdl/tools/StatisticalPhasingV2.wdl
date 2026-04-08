version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "StatisticalPhasing.wdl" as StatPhasing

workflow StatisticalPhasingV2 {
    input {
        File vcf
        File vcf_idx
        File base_vcf
        File base_vcf_idx
        String prefix

        Int operation
        String weight_tag
        Int is_weight_format_field
        Float default_weight
        Boolean remove_duplicates_by_phased_fraction

        String variant_filter_args = "-i 'MAC>=2'"

        File fix_variant_collisions_java

        String utils_docker
        String pysam_docker

        RuntimeAttr? runtime_attr_get_overlapping_samples
        RuntimeAttr? runtime_attr_subset_main_vcf
        RuntimeAttr? runtime_attr_filter_vcf
        RuntimeAttr? runtime_attr_uqids
        RuntimeAttr? runtime_attr_remove_duplicates
        RuntimeAttr? runtime_attr_fill_vcf_tags
        RuntimeAttr? runtime_attr_fix_variant_collisions
        RuntimeAttr? runtime_attr_prepare_base_vcf
        RuntimeAttr? runtime_attr_transfer_haplotypes
    }

    call GetOverlappingSamples {
        input:
            vcf1 = vcf,
            vcf2 = base_vcf,
            prefix = "~{prefix}.overlapping_samples",
            docker = pysam_docker,
            runtime_attr_override = runtime_attr_get_overlapping_samples
    }

    Array[String] overlapping_samples = read_lines(GetOverlappingSamples.samples_file)

    call Helpers.SubsetVcfToSamples as SubsetMainVcf {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            samples = overlapping_samples,
            filter_to_sample = false,
            prefix = "~{prefix}.subset",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_main_vcf
    }

    call StatPhasing.SplitAndFilterVcf {
        input:
            vcf = SubsetMainVcf.subset_vcf,
            vcf_idx = SubsetMainVcf.subset_vcf_idx,
            prefix = "~{prefix}.filtered",
            filter_args = variant_filter_args,
            runtime_attr_override = runtime_attr_filter_vcf
    }

    call StatPhasing.EnsureUniqueIDsAndNormalize {
        input:
            vcf = SplitAndFilterVcf.filtered_vcf,
            vcf_idx = SplitAndFilterVcf.filtered_vcf_idx,
            prefix = "~{prefix}.filtered.uqids",
            docker = pysam_docker,
            runtime_attr_override = runtime_attr_uqids
    }

    if (remove_duplicates_by_phased_fraction) {
        call StatPhasing.RemoveDuplicatesByPhasedFraction {
            input:
                vcf = EnsureUniqueIDsAndNormalize.uqids_vcf,
                vcf_idx = EnsureUniqueIDsAndNormalize.uqids_vcf_idx,
                prefix = "~{prefix}.filtered.uqids",
                docker = pysam_docker,
                runtime_attr_override = runtime_attr_remove_duplicates
        }
    }

    call StatPhasing.FillVcfTags {
        input:
            vcf = select_first([RemoveDuplicatesByPhasedFraction.dups_removed_vcf, EnsureUniqueIDsAndNormalize.uqids_vcf]),
            vcf_idx = select_first([RemoveDuplicatesByPhasedFraction.dups_removed_vcf_idx, EnsureUniqueIDsAndNormalize.uqids_vcf_idx]),
            prefix = "~{prefix}.preprocessed",
            runtime_attr_override = runtime_attr_fill_vcf_tags
    }

    call StatPhasing.FixVariantCollisions {
        input:
            phased_vcf = FillVcfTags.filled_tag_vcf,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            default_weight = default_weight,
            prefix = "~{prefix}.collisionless",
            runtime_attr_override = runtime_attr_fix_variant_collisions
    }

    call PrepareBaseVcf {
        input:
            vcf = base_vcf,
            vcf_idx = base_vcf_idx,
            samples = overlapping_samples,
            prefix = "~{prefix}.base.prepared",
            runtime_attr_override = runtime_attr_prepare_base_vcf
    }

    call TransferHaplotypesFromBaseVcf {
        input:
            original_vcf = FixVariantCollisions.phased_collisionless_vcf,
            original_vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx,
            base_vcf = PrepareBaseVcf.prepared_vcf,
            base_vcf_idx = PrepareBaseVcf.prepared_vcf_idx,
            prefix = "~{prefix}.base_transferred",
            docker = pysam_docker,
            runtime_attr_override = runtime_attr_transfer_haplotypes
    }

    output {
        File uqids_split_vcf = EnsureUniqueIDsAndNormalize.uqids_vcf
        File uqids_split_vcf_idx = EnsureUniqueIDsAndNormalize.uqids_vcf_idx
        File? removed_duplicates_split_vcf = select_first([FillVcfTags.filled_tag_vcf, RemoveDuplicatesByPhasedFraction.dups_removed_vcf])
        File? removed_duplicates_split_vcf_idx = select_first([FillVcfTags.filled_tag_vcf_idx, RemoveDuplicatesByPhasedFraction.dups_removed_vcf_idx])
        File collisionless_split_vcf = FixVariantCollisions.phased_collisionless_vcf
        File collisionless_split_vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx
        File base_prepared_vcf = PrepareBaseVcf.prepared_vcf
        File base_prepared_vcf_idx = PrepareBaseVcf.prepared_vcf_idx
        File base_transferred_vcf = TransferHaplotypesFromBaseVcf.transferred_vcf
        File base_transferred_vcf_idx = TransferHaplotypesFromBaseVcf.transferred_vcf_idx
    }
}

task PrepareBaseVcf {
    input {
        File vcf
        File vcf_idx
        Array[String] samples
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        printf '%s\n' ~{sep=' ' samples} > samples.txt

        # split multiallelics, filter to SNVs only, then subset to overlapping samples
        bcftools norm -m-any ~{vcf} \
            | bcftools view -v snps \
            | bcftools view --samples-file samples.txt \
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GetOverlappingSamples {
    input {
        File vcf1
        File vcf2
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import pysam

        vcf1 = pysam.VariantFile("~{vcf1}", "r")
        vcf2 = pysam.VariantFile("~{vcf2}", "r")
        samples1 = set(vcf1.header.samples)
        samples2 = set(vcf2.header.samples)
        vcf1.close()
        vcf2.close()

        overlapping = sorted(samples1 & samples2)
        with open("~{prefix}.txt", "w") as f:
            for s in overlapping:
                f.write(s + "\n")
        print(f"Found {len(overlapping)} overlapping samples")
        CODE
    >>>

    output {
        File samples_file = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task TransferHaplotypesFromBaseVcf {
    input {
        File original_vcf
        File original_vcf_idx
        File base_vcf
        File base_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import pysam
        from collections import defaultdict

        def transfer_haplotypes_from_base_vcf(original_vcf, base_vcf, output_vcf):
            # load base_vcf: for each sample, index phased het variants by (chrom, pos, ref, alt)
            base_in = pysam.VariantFile(base_vcf, "r")
            samples = list(base_in.header.samples)
            # sample -> variant_key -> (allele0, allele1, is_phased)
            base_gts = defaultdict(dict)
            for rec in base_in:
                if not rec.alts:
                    continue
                key = (rec.contig, rec.pos, rec.ref.upper(), rec.alts[0].upper())
                for sample in samples:
                    gt_data = rec.samples[sample]
                    gt = gt_data.get("GT")
                    if gt is None or None in gt or len(gt) != 2:
                        continue
                    # only store het genotypes
                    if gt[0] == gt[1]:
                        continue
                    base_gts[sample][key] = (gt[0], gt[1], gt_data.phased)
            base_in.close()

            # first pass over original_vcf: collect all records and build per-sample PS block votes
            orig_in = pysam.VariantFile(original_vcf, "r")
            records = []
            # sample -> ps_id -> list of votes (+1 keep, -1 flip)
            sample_ps_votes = defaultdict(lambda: defaultdict(list))

            for rec in orig_in:
                records.append(rec.copy())
                if not rec.alts:
                    continue

                # only use phased het SNVs (allele_type=snv) for voting
                allele_type = rec.info.get("allele_type")
                if allele_type is None or allele_type.lower() != "snv":
                    continue

                key = (rec.contig, rec.pos, rec.ref.upper(), rec.alts[0].upper())

                for sample in samples:
                    gt_data = rec.samples[sample]
                    gt = gt_data.get("GT")
                    if gt is None or None in gt or len(gt) != 2:
                        continue
                    # must be het and phased in original_vcf
                    if gt[0] == gt[1] or not gt_data.phased:
                        continue
                    ps = gt_data.get("PS")
                    if ps is None:
                        continue
                    # check matching het in base_vcf
                    if key not in base_gts[sample]:
                        continue
                    base_a0, base_a1, base_phased = base_gts[sample][key]
                    # only use orientation votes when base_vcf is phased
                    if not base_phased:
                        continue
                    orig_a0, orig_a1 = gt[0], gt[1]
                    if (orig_a0, orig_a1) == (base_a0, base_a1):
                        sample_ps_votes[sample][ps].append(1)
                    elif (orig_a0, orig_a1) == (base_a1, base_a0):
                        sample_ps_votes[sample][ps].append(-1)

            orig_in.close()

            # determine flip decision per (sample, ps_block): flip if majority vote is negative
            flip_map = defaultdict(dict)
            for sample, ps_blocks in sample_ps_votes.items():
                for ps, votes in ps_blocks.items():
                    flip_map[sample][ps] = sum(votes) < 0

            # second pass: write output with flipped haplotypes where needed
            vcf_out = pysam.VariantFile(output_vcf, "w", header=records[0].header if records else pysam.VariantFile(original_vcf, "r").header)
            for rec in records:
                new_rec = rec.copy()
                for sample in samples:
                    gt_data = new_rec.samples[sample]
                    gt = gt_data.get("GT")
                    if gt is None or None in gt or len(gt) != 2:
                        continue
                    ps = gt_data.get("PS")
                    if ps is None or not gt_data.phased:
                        continue
                    if flip_map[sample].get(ps, False):
                        gt_data["GT"] = (gt[1], gt[0])
                        gt_data.phased = True
                vcf_out.write(new_rec)
            vcf_out.close()

            flipped_blocks = sum(v for sample_map in flip_map.values() for v in sample_map.values())
            total_blocks = sum(len(sample_map) for sample_map in flip_map.values())
            print(f"Flipped {flipped_blocks} / {total_blocks} phase blocks (across all samples)")

        transfer_haplotypes_from_base_vcf("~{original_vcf}", "~{base_vcf}", "transferred.vcf")
        CODE

        bgzip transferred.vcf
        bcftools index -t transferred.vcf.gz
        mv transferred.vcf.gz ~{prefix}.vcf.gz
        mv transferred.vcf.gz.tbi ~{prefix}.vcf.gz.tbi
    >>>

    output {
        File transferred_vcf = "~{prefix}.vcf.gz"
        File transferred_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(original_vcf, "GiB") * 3) + ceil(size(base_vcf, "GiB") * 3),
        disk_gb: 15 + ceil(size(original_vcf, "GiB") * 3) + ceil(size(base_vcf, "GiB") * 3),
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
