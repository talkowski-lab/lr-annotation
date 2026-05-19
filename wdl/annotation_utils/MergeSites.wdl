version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeSites {
    input {
        File vcf
        File vcf_idx
        String prefix

        Float del_reciprocal_overlap = 0.5
        Float del_sequence_similarity = 0.5
        Float del_sample_similarity = 0.5
        Int del_size_max = 1000000
        Int del_size_min = 0
        Int del_size_filt = 0

        Float ins_size_similarity = 2.0
        Int ins_breakpoint_distance = 10
        Float ins_sequence_similarity = 0.5
        Float ins_sample_similarity = 0.5
        Int ins_size_max = 1000000
        Int ins_size_min = 0
        Int ins_size_filt = 0

        String utils_docker

        RuntimeAttr? runtime_attr_subset_by_type
        RuntimeAttr? runtime_attr_collapse_dels
        RuntimeAttr? runtime_attr_collapse_ins
        RuntimeAttr? runtime_attr_concat
    }

    # Split VCF into DEL, INS and OTHER
    call Helpers.SubsetVcfByArgs as SubsetDels {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = 'INFO/allele_type~"del"',
            prefix = "~{prefix}.dels",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_by_type
    }

    call Helpers.SubsetVcfByArgs as SubsetIns {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = 'INFO/allele_type~"ins" || INFO/allele_type~"dup"',
            prefix = "~{prefix}.ins",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_by_type
    }

    call Helpers.SubsetVcfByArgs as SubsetOther {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = 'INFO/allele_type!~"del" && INFO/allele_type!~"ins" && INFO/allele_type!~"dup"',
            prefix = "~{prefix}.other",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_by_type
    }

    # Collapse DELs
    call CollapseSites as CollapseDels {
        input:
            vcf = SubsetDels.subset_vcf,
            vcf_idx = SubsetDels.subset_vcf_idx,
            pctovl = del_reciprocal_overlap,
            pctseq = del_sequence_similarity,
            sample_similarity = del_sample_similarity,
            pctsize = 0.0,
            refdist = 0,
            sizemax = del_size_max,
            sizemin = del_size_min,
            sizefilt = del_size_filt,
            variant_type = "DEL",
            prefix = "~{prefix}.dels.collapsed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_collapse_dels
    }

    # Collapse INS
    call CollapseSites as CollapseIns {
        input:
            vcf = SubsetIns.subset_vcf,
            vcf_idx = SubsetIns.subset_vcf_idx,
            pctovl = 0.0,
            pctseq = ins_sequence_similarity,
            sample_similarity = ins_sample_similarity,
            pctsize = ins_size_similarity,
            refdist = ins_breakpoint_distance,
            sizemax = ins_size_max,
            sizemin = ins_size_min,
            sizefilt = ins_size_filt,
            variant_type = "INS",
            prefix = "~{prefix}.ins.collapsed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_collapse_ins
    }

    # Concat collapsed DELs + INS + OTHER
    call Helpers.ConcatVcfs {
        input:
            vcfs = [CollapseDels.collapsed_vcf, CollapseIns.collapsed_vcf, SubsetOther.subset_vcf],
            vcf_idxs = [CollapseDels.collapsed_vcf_idx, CollapseIns.collapsed_vcf_idx, SubsetOther.subset_vcf_idx],
            allow_overlaps = true,
            naive = false,
            prefix = "~{prefix}.merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File merged_vcf = ConcatVcfs.concat_vcf
        File merged_vcf_idx = ConcatVcfs.concat_vcf_idx
        File del_collapse_summary = CollapseDels.summary_tsv
        File ins_collapse_summary = CollapseIns.summary_tsv
    }
}

task CollapseSites {
    input {
        File vcf
        File vcf_idx
        Float pctovl
        Float pctseq
        Float pctsize
        Float sample_similarity
        Int refdist
        Int sizemax
        Int sizemin
        Int sizefilt
        String variant_type
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        truvari collapse \
            -i ~{vcf} \
            -o collapsed.vcf \
            -c removed.vcf \
            --keep common \
            --sizemin 0 \
            --sizemax -1 \
            --pctovl ~{pctovl} \
            --pctseq ~{pctseq} \
            --pctsize ~{pctsize} \
            --refdist ~{refdist} \
            --sizemax ~{sizemax} \
            --sizemin ~{sizemin} \
            --sizefilt ~{sizefilt}

        bgzip -f collapsed.vcf
        bgzip -f removed.vcf

        tabix -f -p vcf collapsed.vcf.gz

        # Post-process: check sample overlap and consolidate genotypes
        python3 <<CODE
import pysam
from collections import defaultdict

def get_nonref_samples(record):
    nonref = set()
    for sample in record.samples:
        gt = record.samples[sample]['GT']
        if gt is not None and any(a is not None and a != 0 for a in gt):
            nonref.add(sample)
    return nonref

def sample_overlap(set_a, set_b):
    if not set_a and not set_b:
        return 1.0
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)

sample_sim_threshold = ~{sample_similarity}

# Index original input VCF by variant ID for clean genotypes
original_records = {}
with pysam.VariantFile("~{vcf}") as orig_vcf:
    orig_header = orig_vcf.header.copy()
    for record in orig_vcf:
        original_records[record.id] = record

# Build mapping of CollapseId -> list of removed record IDs
# MatchId is Number=. so pysam returns it as a tuple
collapse_groups = defaultdict(list)
with pysam.VariantFile("removed.vcf.gz") as removed_vcf:
    for record in removed_vcf:
        match_id = record.info.get("MatchId", None)
        if match_id:
            key = match_id[0] if isinstance(match_id, (tuple, list)) else match_id
            collapse_groups[key].append(record.id)

# Process kept variants: check sample overlap, consolidate or reject
kept_vcf = pysam.VariantFile("collapsed.vcf.gz")
out_vcf = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=orig_header)

for kept_record in kept_vcf:
    collapse_id = kept_record.info.get("CollapseId", None)

    if collapse_id and collapse_id in collapse_groups:
        removed_ids = collapse_groups[collapse_id]
        orig_kept = original_records.get(kept_record.id)
        kept_nonref = get_nonref_samples(orig_kept) if orig_kept else get_nonref_samples(kept_record)

        for removed_id in removed_ids:
            orig_removed = original_records.get(removed_id)
            if not orig_removed:
                continue
            removed_nonref = get_nonref_samples(orig_removed)
            overlap = sample_overlap(kept_nonref, removed_nonref)

            if overlap >= sample_sim_threshold:
                # Valid collapse: consolidate non-ref genotypes into kept
                for sample in orig_kept.samples:
                    kept_gt = orig_kept.samples[sample]['GT']
                    collapsed_gt = orig_removed.samples[sample]['GT']
                    if (kept_gt is None or all(a is None or a == 0 for a in kept_gt)):
                        if collapsed_gt is not None and any(a is not None and a != 0 for a in collapsed_gt):
                            orig_kept.samples[sample]['GT'] = collapsed_gt
            else:
                # Failed sample overlap: write removed variant back as separate site
                out_vcf.write(orig_removed)

    # Write the kept variant (with consolidated GTs) using original record
    if kept_record.id in original_records:
        out_vcf.write(original_records[kept_record.id])

kept_vcf.close()
out_vcf.close()
CODE

        # Sort and index the final output
        bcftools sort -T . -Oz -o sorted.vcf.gz ~{prefix}.vcf.gz
        mv sorted.vcf.gz ~{prefix}.vcf.gz
        tabix -f -p vcf ~{prefix}.vcf.gz

        n_output=$(bcftools view -H ~{prefix}.vcf.gz | wc -l | awk '{print $1}')
        n_collapsed=$((n_input - n_output))

        echo -e "variant_type\tn_input\tn_output\tn_collapsed" > ~{prefix}.summary.tsv
        echo -e "~{variant_type}\t${n_input}\t${n_output}\t${n_collapsed}" >> ~{prefix}.summary.tsv
    >>>

    output {
        File collapsed_vcf = "~{prefix}.vcf.gz"
        File collapsed_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File summary_tsv = "~{prefix}.summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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
