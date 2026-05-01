version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeVcfs {
    input {
        Array[File] contig_vcfs
        Array[File] contig_vcf_idxs
        String contig

        Int min_truvari_match = 50
        Int? shard_bin_size
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_baseline_merge
        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_shard
        RuntimeAttr? runtime_attr_merge_per_shard
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_shard_summaries
    }

    if (defined(shard_bin_size)) {
        call Helpers.CreateContigShards {
            input:
                vcfs = contig_vcfs,
                vcf_idxs = contig_vcf_idxs,
                contig = contig,
                shard_bin_size = select_first([shard_bin_size]),
                prefix = "~{prefix}.shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_shards
        }

        scatter (j in range(length(CreateContigShards.shard_regions))) {
            scatter (i in range(length(contig_vcfs))) {
                call Helpers.SubsetVcfToRegion as SubsetCallsetToShard {
                    input:
                        vcf = contig_vcfs[i],
                        vcf_idx = contig_vcf_idxs[i],
                        region = CreateContigShards.shard_regions[j],
                        prefix = "~{prefix}.shard_~{j}.callset_~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_shard
                }
            }

            call BaselineMergeVcfs as BaselineMergeShard {
                input:
                    vcfs = SubsetCallsetToShard.subset_vcf,
                    vcf_idxs = SubsetCallsetToShard.subset_vcf_idx,
                    prefix = "~{prefix}.shard_~{j}.baseline",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_baseline_merge
            }

            call MergeVcfsPerShard as MergeShard {
                input:
                    vcf = BaselineMergeShard.merged_vcf,
                    vcf_idx = BaselineMergeShard.merged_vcf_idx,
                    contig = contig,
                    min_truvari_match = min_truvari_match,
                    prefix = "~{prefix}.shard_~{j}.merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_per_shard
            }
        }

        call Helpers.ConcatVcfs as ConcatShards {
            input:
                vcfs = MergeShard.merged_vcf,
                vcf_idxs = MergeShard.merged_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shards
        }

        call Helpers.ConcatTsvs as ConcatShardSummaries {
            input:
                tsvs = MergeShard.summary_tsv,
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.merge_summary",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shard_summaries
        }
    }

    if (!defined(shard_bin_size)) {
        call BaselineMergeVcfs as BaselineMergeContig {
            input:
                vcfs = contig_vcfs,
                vcf_idxs = contig_vcf_idxs,
                prefix = "~{prefix}.baseline",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_baseline_merge
        }

        call MergeVcfsPerShard as MergeWholeContig {
            input:
                vcf = BaselineMergeContig.merged_vcf,
                vcf_idx = BaselineMergeContig.merged_vcf_idx,
                contig = contig,
                min_truvari_match = min_truvari_match,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_per_shard
        }
    }

    output {
        File merged_vcf = select_first([ConcatShards.concat_vcf, MergeWholeContig.merged_vcf])
        File merged_vcf_idx = select_first([ConcatShards.concat_vcf_idx, MergeWholeContig.merged_vcf_idx])
        File merge_summary_tsv = select_first([ConcatShardSummaries.concatenated_tsv, MergeWholeContig.summary_tsv])
    }
}

task BaselineMergeVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        paste ~{write_lines(vcfs)} ~{write_lines(vcf_idxs)} > vcf_pairs.tsv

        i=0
        while IFS=$'\t' read -r vcf vcf_idx; do
            if [[ "$vcf_idx" != "${vcf}.tbi" ]]; then
                ln -sf "$vcf_idx" "${vcf}.tbi"
            fi

            bcftools annotate \
                -x FORMAT/AL \
                -Oz -o "cleaned_${i}.vcf.gz" \
                "$vcf"

            tabix -f -p vcf "cleaned_${i}.vcf.gz"
            echo "cleaned_${i}.vcf.gz" >> cleaned_vcfs.list
            i=$((i + 1))
        done < vcf_pairs.tsv

        bcftools merge \
            -m none \
            -Oz -o merged.unsorted.vcf.gz \
            -l cleaned_vcfs.list

        bcftools sort \
            -T . \
            -Oz -o ~{prefix}.vcf.gz \
            merged.unsorted.vcf.gz

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 4 * ceil(size(vcfs, "GB")) + 20,
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

task MergeVcfsPerShard {
    input {
        File vcf
        File vcf_idx
        String contig
        Int min_truvari_match
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        VCF="~{vcf}"
        VCF_IDX="~{vcf_idx}"
        if [[ "$VCF_IDX" != "${VCF}.tbi" ]]; then
            ln -sf "$VCF_IDX" "${VCF}.tbi"
        fi

        bcftools view -i 'INFO/allele_type=="trv"' -Oz -o trv.vcf.gz "$VCF"
        bcftools view -i 'INFO/allele_type!="trv"' -Oz -o non_trv.vcf.gz "$VCF"
        tabix -f -p vcf trv.vcf.gz
        tabix -f -p vcf non_trv.vcf.gz

        python3 <<'PYEOF'
import os
import pysam
import subprocess
from collections import defaultdict

CONTIG = "~{contig}"
MIN_TRUVARI = ~{min_truvari_match}
PREFIX = "~{prefix}"
TRV_VCF = "trv.vcf.gz"
NON_TRV_VCF = "non_trv.vcf.gz"

def read_records(path):
    with pysam.VariantFile(path) as vcf_in:
        header = vcf_in.header.copy()
        samples = list(vcf_in.header.samples)
        records = [record.copy() for record in vcf_in]
    return header, samples, records

def get_allele_length(record):
    val = record.info.get("allele_length", 0)
    if isinstance(val, (list, tuple)):
        val = val[0]
    if val in (None, ".", ""):
        return 0
    return int(val)

def get_allele_type(record):
    val = record.info.get("allele_type", "")
    if isinstance(val, (list, tuple)):
        val = val[0]
    return str(val)

def merge_info(records):
    merged = {}
    for record in records:
        for key in record.info:
            if key not in merged:
                val = record.info[key]
                merged[key] = list(val) if isinstance(val, tuple) else val

    return merged

def recalc_ac_an_af(record):
    an = 0
    ac_counts = defaultdict(int)
    for sample in record.samples:
        gt = record.samples[sample]["GT"]
        if gt is None:
            continue
        for allele in gt:
            if allele is None:
                continue
            an += 1
            if allele > 0:
                ac_counts[allele] += 1

    n_alts = len(record.alts) if record.alts else 1
    ac = [ac_counts.get(i + 1, 0) for i in range(n_alts)]
    af = [c / an if an > 0 else 0.0 for c in ac]

    if "AN" in record.header.info:
        record.info["AN"] = an
    if "AC" in record.header.info:
        record.info["AC"] = tuple(ac)
    if "AF" in record.header.info:
        record.info["AF"] = tuple(af)

def group_by_key(records, key_fn):
    groups = defaultdict(list)
    order = []
    seen = set()
    for rec in records:
        k = key_fn(rec)
        if k not in seen:
            order.append(k)
            seen.add(k)
        groups[k].append(rec)
    return order, groups

def build_merged_record(group, header, samples, combined_alts=None):
    rep = group[0]

    if combined_alts is not None:
        all_alts = combined_alts
    else:
        all_alts = list(rep.alts) if rep.alts else ["."]

    out = header.new_record()
    out.chrom = rep.chrom
    out.pos = rep.pos
    out.id = rep.id
    out.ref = rep.ref
    out.alts = tuple(all_alts)
    out.qual = rep.qual
    out.filter.add("PASS")

    merged_info = merge_info(group)
    for key, val in merged_info.items():
        if key in header.info:
            try:
                if isinstance(val, list):
                    out.info[key] = tuple(val)
                else:
                    out.info[key] = val
            except Exception:
                pass

    alt_index_maps = []
    for rec in group:
        if rec.alts:
            idx_map = {}
            for i, a in enumerate(rec.alts):
                if a in all_alts:
                    idx_map[i + 1] = all_alts.index(a) + 1
                else:
                    idx_map[i + 1] = i + 1
        else:
            idx_map = {}
        alt_index_maps.append(idx_map)

    sample_gts = {}
    for rec, idx_map in zip(group, alt_index_maps):
        for sample in rec.samples:
            gt = rec.samples[sample].get("GT")
            if gt is None or all(a is None for a in gt):
                continue
            remapped = tuple(
                (idx_map.get(a, a) if a is not None and a > 0 else a)
                for a in gt
            )
            sample_gts[sample] = remapped

    for sample in samples:
        gt = sample_gts.get(sample, (None, None))
        try:
            out.samples[sample]["GT"] = gt
        except Exception:
            pass

    recalc_ac_an_af(out)
    return out

header, samples, trv_records = read_records(TRV_VCF)
_, _, non_trv_records = read_records(NON_TRV_VCF)

n_input = len(trv_records) + len(non_trv_records)

def trv_key(rec):
    return (rec.chrom, rec.pos, rec.ref)

trv_order, trv_groups = group_by_key(trv_records, trv_key)

merged_trv = []
for k in trv_order:
    group = trv_groups[k]
    seen_alts = []
    seen_set = set()
    for rec in group:
        for a in (rec.alts or []):
            if a not in seen_set:
                seen_alts.append(a)
                seen_set.add(a)
    merged_trv.append(build_merged_record(group, header, samples, combined_alts=seen_alts))

def exact_key(rec):
    alts = ",".join(rec.alts) if rec.alts else "."
    return (rec.chrom, rec.pos, rec.ref, alts)

exact_order, exact_groups = group_by_key(non_trv_records, exact_key)

exact_merged_non_trv = []
for k in exact_order:
    exact_merged_non_trv.append(build_merged_record(exact_groups[k], header, samples))

small_non_trv = []
truvari_inputs = []
for record in exact_merged_non_trv:
    if abs(get_allele_length(record)) >= MIN_TRUVARI:
        truvari_inputs.append(record)
    else:
        small_non_trv.append(record)

merged_truvari = []
n_truvari_input = len(truvari_inputs)

if truvari_inputs:
    tmp_vcf_path = "tmp_truvari_input.vcf.gz"
    tmp_out = pysam.VariantFile(tmp_vcf_path, "w", header=header)
    for rec in truvari_inputs:
        tmp_out.write(rec)
    tmp_out.close()
    subprocess.run(["tabix", "-f", "-p", "vcf", tmp_vcf_path], check=True)

    kept_out = "tmp_truvari_kept.vcf"
    collapsed_out = "tmp_truvari_removed.vcf"
    subprocess.run([
        "truvari", "collapse",
        "-i", tmp_vcf_path,
        "-o", kept_out,
        "-c", collapsed_out,
        "--sizemin", "0",
    ], check=True)

    subprocess.run(["bgzip", "-f", kept_out], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", kept_out + ".gz"], check=True)

    orig_by_id = {}
    for rec in truvari_inputs:
        if rec.id:
            orig_by_id[rec.id] = rec

    collapse_to_kept = defaultdict(list)
    if os.path.exists(collapsed_out):
        c_vcf = pysam.VariantFile(collapsed_out)
        for rec in c_vcf:
            cid = rec.info.get("CollapseId", None)
            if cid is not None and rec.id:
                collapse_to_kept[str(cid)].append(rec.id)
        c_vcf.close()

    kept_vcf = pysam.VariantFile(kept_out + ".gz")
    for rec in kept_vcf:
        cid = rec.info.get("CollapseId", None)
        member_recs = []
        if rec.id and rec.id in orig_by_id:
            member_recs.append(orig_by_id[rec.id])
        if cid is not None:
            for mid in collapse_to_kept.get(str(cid), []):
                if mid in orig_by_id:
                    member_recs.append(orig_by_id[mid])
        if not member_recs:
            member_recs = [rec]

        merged_info = merge_info(member_recs)
        for key, val in merged_info.items():
            if key in header.info:
                try:
                    if isinstance(val, list):
                        rec.info[key] = tuple(val)
                    else:
                        rec.info[key] = val
                except Exception:
                    pass
        recalc_ac_an_af(rec)
        merged_truvari.append(rec.copy())
    kept_vcf.close()

all_merged = merged_trv + small_non_trv + merged_truvari
all_merged.sort(key=lambda r: (r.chrom, r.pos, r.ref, ",".join(r.alts or ["."])))

out_vcf_path = PREFIX + ".vcf.gz"
vcf_out = pysam.VariantFile(out_vcf_path, "w", header=header)
for rec in all_merged:
    vcf_out.write(rec)
vcf_out.close()

subprocess.run(["tabix", "-f", "-p", "vcf", out_vcf_path], check=True)

def count_by_type(records):
    counts = defaultdict(int)
    for rec in records:
        counts[get_allele_type(rec)] += 1
    return counts

n_trv_input = len(trv_records)
n_trv_merged_out = len(merged_trv)
n_trv_collapsed = n_trv_input - n_trv_merged_out

n_exact_input = len(non_trv_records)
n_exact_merged_out = len(exact_merged_non_trv)
n_exact_collapsed = n_exact_input - n_exact_merged_out

n_truvari_merged = len(merged_truvari)
n_truvari_collapsed = n_truvari_input - n_truvari_merged

with open(PREFIX + ".merge_summary.tsv", "w") as f:
    f.write("\t".join([
        "contig",
        "n_input_variants",
        "n_output_variants",
        "n_collapsed_total",
        "n_trv_input",
        "n_trv_output",
        "n_trv_collapsed",
        "n_exact_input",
        "n_exact_output",
        "n_exact_collapsed",
        "n_truvari_input",
        "n_truvari_output",
        "n_truvari_collapsed",
        "input_by_allele_type",
    ]) + "\n")

    combined_type_counts = count_by_type(trv_records + non_trv_records)
    type_str = ";".join(f"{t}:{c}" for t, c in sorted(combined_type_counts.items()))

    f.write("\t".join([
        CONTIG,
        str(n_input),
        str(len(all_merged)),
        str(n_input - len(all_merged)),
        str(n_trv_input),
        str(n_trv_merged_out),
        str(n_trv_collapsed),
        str(n_exact_input),
        str(n_exact_merged_out),
        str(n_exact_collapsed),
        str(n_truvari_input),
        str(n_truvari_merged),
        str(n_truvari_collapsed),
        type_str,
    ]) + "\n")

print(f"[MergeVcfsPerShard] {CONTIG}: {n_input} input -> {len(all_merged)} output variants")
PYEOF
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File summary_tsv = "~{prefix}.merge_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 4 * ceil(size(vcf, "GB")) + 20,
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
