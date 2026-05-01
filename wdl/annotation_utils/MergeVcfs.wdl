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

        RuntimeAttr? runtime_attr_bcftools_merge
        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_shard
        RuntimeAttr? runtime_attr_merge_per_contig
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_shard_summaries
    }

    call Helpers.MergeVcfs as MergeCallsets {
        input:
            vcfs = contig_vcfs,
            vcf_idxs = contig_vcf_idxs,
            contig = contig,
            extra_args = "-m none",
            sort_merged = true,
            prefix = "~{prefix}.callsets_merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_bcftools_merge
    }

    if (defined(shard_bin_size)) {
        call Helpers.CreateContigShards {
            input:
                vcfs = [MergeCallsets.merged_vcf],
                vcf_idxs = [MergeCallsets.merged_vcf_idx],
                contig = contig,
                shard_bin_size = select_first([shard_bin_size]),
                prefix = "~{prefix}.shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_shards
        }

        scatter (j in range(length(CreateContigShards.shard_regions))) {
            call Helpers.SubsetVcfToRegion as SubsetShard {
                input:
                    vcf = MergeCallsets.merged_vcf,
                    vcf_idx = MergeCallsets.merged_vcf_idx,
                    region = CreateContigShards.shard_regions[j],
                    prefix = "~{prefix}.shard_~{j}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_shard
            }

            call MergeVcfsPerContig as MergeVcfsPerShard {
                input:
                    vcf = SubsetShard.subset_vcf,
                    vcf_idx = SubsetShard.subset_vcf_idx,
                    contig = contig,
                    min_truvari_match = min_truvari_match,
                    prefix = "~{prefix}.shard_~{j}.merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_per_contig
            }
        }

        call Helpers.ConcatVcfs as ConcatShards {
            input:
                vcfs = MergeVcfsPerShard.merged_vcf,
                vcf_idxs = MergeVcfsPerShard.merged_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shards
        }

        call Helpers.ConcatTsvs as ConcatShardSummaries {
            input:
                tsvs = MergeVcfsPerShard.summary_tsv,
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.merge_summary",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shard_summaries
        }
    }

    if (!defined(shard_bin_size)) {
        call MergeVcfsPerContig {
            input:
                vcf = MergeCallsets.merged_vcf,
                vcf_idx = MergeCallsets.merged_vcf_idx,
                contig = contig,
                min_truvari_match = min_truvari_match,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_per_contig
        }
    }

    output {
        File merged_vcf = select_first([ConcatShards.concat_vcf, MergeVcfsPerContig.merged_vcf])
        File merged_vcf_idx = select_first([ConcatShards.concat_vcf_idx, MergeVcfsPerContig.merged_vcf_idx])
        File merge_summary_tsv = select_first([ConcatShardSummaries.concatenated_tsv, MergeVcfsPerContig.summary_tsv])
    }
}

task MergeVcfsPerContig {
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

        # Ensure index is co-located
        VCF="~{vcf}"
        VCF_IDX="~{vcf_idx}"
        if [[ "$VCF_IDX" != "${VCF}.tbi" ]]; then
            ln -sf "$VCF_IDX" "${VCF}.tbi"
        fi

        python3 <<'PYEOF'
import pysam
import subprocess
import os
import sys
from collections import defaultdict

VCF_PATH = "~{vcf}"
CONTIG = "~{contig}"
MIN_TRUVARI = ~{min_truvari_match}
PREFIX = "~{prefix}"

# ── helpers ─────────────────────────────────────────────────────────────────

BASIC_TYPES = {"ins", "del", "snv"}

def get_allele_length(record):
    val = record.info.get("allele_length", 0)
    if isinstance(val, (list, tuple)):
        val = val[0]
    return int(val)

def get_allele_type(record):
    val = record.info.get("allele_type", "")
    if isinstance(val, (list, tuple)):
        val = val[0]
    return str(val)

def pick_allele_type(types):
    """Return a non-basic type if present, else any."""
    for t in types:
        if t.lower() not in BASIC_TYPES:
            return t
    return types[0]

def merge_info(records, header):
    """Merge INFO fields across records, with allele_length and allele_type logic."""
    merged = {}
    for record in records:
        for key in record.info:
            if key not in merged:
                val = record.info[key]
                merged[key] = list(val) if isinstance(val, tuple) else val

    # allele_length: largest value
    lengths = []
    for record in records:
        if "allele_length" in record.info:
            lengths.append(get_allele_length(record))
    if lengths:
        merged["allele_length"] = max(lengths, key=abs)

    # allele_type: non-basic wins
    types = [get_allele_type(r) for r in records if "allele_type" in r.info]
    if types:
        merged["allele_type"] = pick_allele_type(types)

    return merged

def recalc_ac_an_af(record):
    """Recalculate AC, AN, AF from genotypes in-place."""
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

# ── read all records, tag with source index ─────────────────────────────────

vcf_in = pysam.VariantFile(VCF_PATH)
header = vcf_in.header
samples = list(header.samples)

all_records = []
for record in vcf_in.fetch(CONTIG):
    all_records.append(record.copy())
vcf_in.close()

n_input = len(all_records)

# Partition
trv_records = []       # allele_type == "trv"
exact_records = []     # non-trv, abs(allele_length) < MIN_TRUVARI
truvari_records = []   # non-trv, abs(allele_length) >= MIN_TRUVARI

for rec in all_records:
    at = get_allele_type(rec)
    al = abs(get_allele_length(rec))
    if at == "trv":
        trv_records.append(rec)
    elif al < MIN_TRUVARI:
        exact_records.append(rec)
    else:
        truvari_records.append(rec)

# ── merge helpers ────────────────────────────────────────────────────────────

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

def build_merged_record(group, header, combined_alts=None):
    """
    Merge a group of records into one. combined_alts overrides ALT list (for trv).
    All sample GTs are placed into the single output record.
    Returns a new VariantRecord.
    """
    rep = group[0]

    # Build unified allele list
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

    # Merge INFO
    merged_info = merge_info(group, header)
    for key, val in merged_info.items():
        if key in header.info:
            try:
                if isinstance(val, list):
                    out.info[key] = tuple(val)
                else:
                    out.info[key] = val
            except Exception:
                pass

    # Build per-record alt-index mapping (only needed for trv where alts differ)
    alt_index_maps = []
    for rec in group:
        if rec.alts:
            idx_map = {}
            for i, a in enumerate(rec.alts):
                if a in all_alts:
                    idx_map[i + 1] = all_alts.index(a) + 1
                else:
                    idx_map[i + 1] = i + 1  # fallback
        else:
            idx_map = {}
        alt_index_maps.append(idx_map)

    # Collect sample GTs: last non-missing value wins
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
        gt = sample_gts.get(sample, (None,))
        try:
            out.samples[sample]["GT"] = gt
        except Exception:
            pass

    recalc_ac_an_af(out)
    return out

# ── process trv (match on CHROM+POS+REF) ─────────────────────────────────────

def trv_key(rec):
    return (rec.chrom, rec.pos, rec.ref)

trv_order, trv_groups = group_by_key(trv_records, trv_key)

merged_trv = []
for k in trv_order:
    group = trv_groups[k]
    # Combined ALTs: union preserving order
    seen_alts = []
    seen_set = set()
    for rec in group:
        for a in (rec.alts or []):
            if a not in seen_set:
                seen_alts.append(a)
                seen_set.add(a)
    merged_trv.append(build_merged_record(group, header, combined_alts=seen_alts))

# ── process exact non-trv (match on CHROM+POS+REF+ALT) ───────────────────────

def exact_key(rec):
    alts = ",".join(rec.alts) if rec.alts else "."
    return (rec.chrom, rec.pos, rec.ref, alts)

exact_order, exact_groups = group_by_key(exact_records, exact_key)

merged_exact = []
for k in exact_order:
    merged_exact.append(build_merged_record(exact_groups[k], header))

# ── process truvari non-trv (truvari collapse) ────────────────────────────────

merged_truvari = []
n_truvari_input = len(truvari_records)
n_truvari_merged = 0

if truvari_records:
    # Write subset VCF for truvari collapse
    tmp_vcf_path = "tmp_truvari_input.vcf.gz"
    tmp_out = pysam.VariantFile(tmp_vcf_path, "w", header=header)
    for rec in truvari_records:
        tmp_out.write(rec)
    tmp_out.close()
    subprocess.run(["tabix", "-p", "vcf", tmp_vcf_path], check=True)

    # Run truvari collapse; -o receives kept variants (with merged GTs), -c receives collapsed-away variants
    # --sizemin 0 / --sizefilt 0: disable truvari's internal size floor (default 50bp) so all
    # variants we pass in are analyzed, consistent with TruvariMatch.wdl pattern in this repo
    kept_out = "tmp_truvari_kept.vcf"
    collapsed_out = "tmp_truvari_removed.vcf"
    subprocess.run([
        "truvari", "collapse",
        "-i", tmp_vcf_path,
        "-o", kept_out,
        "-c", collapsed_out,
        "--sizemin", "0",
        "--sizefilt", "0",
    ], check=True)

    # Compress and index kept output
    subprocess.run(["bgzip", kept_out], check=True)
    subprocess.run(["tabix", "-p", "vcf", kept_out + ".gz"], check=True)

    # Build index from original records by ID for INFO merging
    orig_by_id = {}
    for rec in truvari_records:
        if rec.id:
            orig_by_id[rec.id] = rec

    # Read collapsed-away file to get CollapseId -> kept variant mapping
    # CollapseId in the collapsed file matches CollapseId in the kept file
    collapse_to_kept = defaultdict(list)
    if os.path.exists(collapsed_out):
        c_vcf = pysam.VariantFile(collapsed_out)
        for rec in c_vcf.fetch():
            cid = rec.info.get("CollapseId", None)
            if cid is not None:
                collapse_to_kept[str(cid)].append(rec.id)
        c_vcf.close()

    # Read kept variants; merge INFO from all originally matching records
    kept_vcf = pysam.VariantFile(kept_out + ".gz")
    for rec in kept_vcf.fetch():
        cid = rec.info.get("CollapseId", None)
        # Gather original members: the kept variant + any collapsed-away variants
        member_recs = []
        if rec.id and rec.id in orig_by_id:
            member_recs.append(orig_by_id[rec.id])
        if cid is not None:
            for mid in collapse_to_kept.get(str(cid), []):
                if mid in orig_by_id:
                    member_recs.append(orig_by_id[mid])
        if not member_recs:
            member_recs = [rec]

        # Apply INFO merge logic on top of the truvari-collapsed record (which already has merged GTs)
        merged_info = merge_info(member_recs, header)
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

    n_truvari_merged = len(merged_truvari)

# ── write merged VCF ─────────────────────────────────────────────────────────

all_merged = merged_trv + merged_exact + merged_truvari

# Sort by pos
all_merged.sort(key=lambda r: (r.chrom, r.pos))

out_vcf_path = PREFIX + ".vcf.gz"
vcf_out = pysam.VariantFile(out_vcf_path, "w", header=header)
for rec in all_merged:
    vcf_out.write(rec)
vcf_out.close()

subprocess.run(["tabix", "-p", "vcf", out_vcf_path], check=True)

# ── summary statistics ────────────────────────────────────────────────────────

def count_by_type(records):
    counts = defaultdict(int)
    for rec in records:
        counts[get_allele_type(rec)] += 1
    return counts

# "exclusive to vcf_0" = variants that are not merged with any other record
# (groups of size 1 after merging). We track the source vcf index via INFO source tag
# Since the input is a single merged VCF (distinct samples already combined by caller),
# we compute per-allele-type merge counts.

n_trv_input = len(trv_records)
n_trv_merged_out = len(merged_trv)
n_trv_collapsed = n_trv_input - n_trv_merged_out  # merged away records

n_exact_input = len(exact_records)
n_exact_merged_out = len(merged_exact)
n_exact_collapsed = n_exact_input - n_exact_merged_out

n_truvari_collapsed = n_truvari_input - n_truvari_merged if truvari_records else 0

# Per-allele-type input breakdown
trv_type_counts = count_by_type(trv_records)
exact_type_counts = count_by_type(exact_records)
truvari_type_counts = count_by_type(truvari_records)

# Write summary
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

    # Aggregate allele_type input counts
    combined_type_counts = defaultdict(int)
    for d in [trv_type_counts, exact_type_counts, truvari_type_counts]:
        for t, c in d.items():
            combined_type_counts[t] += c
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

print(f"[MergeVcfsPerContig] {CONTIG}: {n_input} input -> {len(all_merged)} output variants")
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
