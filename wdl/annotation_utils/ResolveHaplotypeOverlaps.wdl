version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ResolveHaplotypeOverlaps {
    input {
        File vcf
        File vcf_idx
        String contig
        String prefix

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_get_samples
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_detect_overlaps
        RuntimeAttr? runtime_attr_concat_tsvs
        RuntimeAttr? runtime_attr_shard_vcf
        RuntimeAttr? runtime_attr_clear_overlaps
        RuntimeAttr? runtime_attr_concat_shards
    }

    call Helpers.GetSamplesFromVcf {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_get_samples
    }

    scatter (sample in GetSamplesFromVcf.samples) {
        call Helpers.ExtractSample {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                sample = sample,
                extra_args = "-e 'INFO/allele_type=\"trv\" || INFO/TR_ENVELOPED=1'",
                prefix = "~{prefix}.~{sample}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_sample
        }

        call DetectSampleOverlaps {
            input:
                vcf = ExtractSample.subset_vcf,
                vcf_idx = ExtractSample.subset_vcf_idx,
                sample = sample,
                prefix = "~{prefix}.~{sample}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_detect_overlaps
        }
    }

    call Helpers.ConcatTsvs as ConcatOverlapTsvs {
        input:
            tsvs = DetectSampleOverlaps.overlaps_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.overlaps",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    call Helpers.ConcatTsvs as ConcatClearedTsvs {
        input:
            tsvs = DetectSampleOverlaps.cleared_calls_tsv,
            sort_output = false,
            preserve_header = false,
            prefix = "~{prefix}.cleared_calls",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    if (defined(records_per_shard)) {
        call Helpers.ShardVcfByRecords {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                records_per_shard = select_first([records_per_shard]),
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard_vcf
        }

        scatter (shard_pair in zip(ShardVcfByRecords.shards, ShardVcfByRecords.shard_idxs)) {
            call ClearOverlapCalls as ClearShard {
                input:
                    vcf = shard_pair.left,
                    vcf_idx = shard_pair.right,
                    cleared_calls_tsv = ConcatClearedTsvs.concatenated_tsv,
                    prefix = "~{prefix}.~{contig}.shard.cleared",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_clear_overlaps
            }
        }

        call Helpers.ConcatVcfs as ConcatShards {
            input:
                vcfs = ClearShard.output_vcf,
                vcf_idxs = ClearShard.output_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.~{contig}.cleared",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shards
        }
    }

    if (!defined(records_per_shard)) {
        call ClearOverlapCalls as ClearAll {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                cleared_calls_tsv = ConcatClearedTsvs.concatenated_tsv,
                prefix = "~{prefix}.~{contig}.cleared",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_clear_overlaps
        }
    }

    output {
        File overlap_resolved_vcf = select_first([ConcatShards.concat_vcf, ClearAll.output_vcf])
        File overlap_resolved_vcf_idx = select_first([ConcatShards.concat_vcf_idx, ClearAll.output_vcf_idx])
        File overlap_tsv = ConcatOverlapTsvs.concatenated_tsv
    }
}


task DetectSampleOverlaps {
    input {
        File vcf
        File vcf_idx
        String sample
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import subprocess

TYPE_RANK = {'DEL': 2, 'INS': 1, 'SNV': 0, 'OTHER': -1}


def classify_type(at, al_str):
    at = str(at).lower()
    try:
        al = abs(int(al_str))
    except (ValueError, TypeError):
        al = None
    if al == 0:
        return 'SNV'
    if 'del' in at:
        return 'DEL'
    if 'ins' in at or 'dup' in at:
        return 'INS'
    return 'OTHER'


def get_length_bin(al_str):
    try:
        al = abs(int(al_str))
    except (ValueError, TypeError):
        return 'unknown'
    if al == 0:
        return 'SNV'
    elif al <= 29:
        return '1-29bp'
    elif al <= 49:
        return '30-49bp'
    elif al <= 500:
        return '50-500bp'
    elif al <= 5000:
        return '500bp-5kbp'
    elif al <= 50000:
        return '5kbp-50kbp'
    else:
        return '50kbp+'


def wins(iv_i, iv_j):
    """Full tie-breaking cascade: ref span → GQ → allele_length → type rank → QUAL → file order."""
    _, _, t_i, _, _, ref_i, gq_i, al_i, qual_i, idx_i = iv_i
    _, _, t_j, _, _, ref_j, gq_j, al_j, qual_j, idx_j = iv_j
    if len(ref_i) != len(ref_j):
        return len(ref_i) > len(ref_j)
    if gq_i != gq_j:
        return gq_i > gq_j
    if al_i != al_j:
        return al_i > al_j
    ri, rj = TYPE_RANK.get(t_i, -1), TYPE_RANK.get(t_j, -1)
    if ri != rj:
        return ri > rj
    if qual_i != qual_j:
        return qual_i > qual_j
    return idx_i < idx_j


sample_name = '~{sample}'

proc = subprocess.Popen(
    ['bcftools', 'query',
     '-f', '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/allele_type\t%INFO/allele_length\t%QUAL[\t%GT\t%GQ]\n',
     '~{vcf}'],
    stdout=subprocess.PIPE,
    text=True,
    bufsize=1
)

hap0 = []
hap1 = []
n_records = 0

for line in proc.stdout:
    row = line.rstrip('\n').split('\t')
    if len(row) < 10:
        continue
    pos1 = int(row[1])
    ref = row[2].upper()
    vid = row[4]
    at = row[5]
    al_str = row[6]
    qual_str = row[7]
    gt_str = row[8]
    gq_str = row[9]
    if '|' not in gt_str:
        continue
    parts = gt_str.split('|')
    try:
        a0 = int(parts[0]) if parts[0] not in ('.', '') else -1
        a1 = int(parts[1]) if parts[1] not in ('.', '') else -1
    except (ValueError, IndexError):
        continue
    try:
        gq = int(gq_str)
    except (ValueError, TypeError):
        gq = 0
    try:
        al_int = abs(int(al_str))
    except (ValueError, TypeError):
        al_int = 0
    try:
        qual = float(qual_str)
    except (ValueError, TypeError):
        qual = 0.0
    iv = (pos1, pos1 + len(ref) - 1, classify_type(at, al_str), get_length_bin(al_str),
          vid, ref, gq, al_int, qual, n_records)
    if a0 > 0:
        hap0.append(iv)
    if a1 > 0:
        hap1.append(iv)
    n_records += 1

proc.wait()

overlaps_rows = []
cleared_set = set()


def find_overlaps(intervals, haplotype):
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    n = len(sorted_ivs)
    for i in range(n):
        s_i, e_i = sorted_ivs[i][0], sorted_ivs[i][1]
        for j in range(i + 1, n):
            if sorted_ivs[j][0] > e_i:
                break
            vid_i, t_i, lb_i = sorted_ivs[i][4], sorted_ivs[i][2], sorted_ivs[i][3]
            vid_j, t_j, lb_j = sorted_ivs[j][4], sorted_ivs[j][2], sorted_ivs[j][3]
            if wins(sorted_ivs[i], sorted_ivs[j]):
                winner_id, winner_type, winner_lb = vid_i, t_i, lb_i
                loser_id, loser_type, loser_lb = vid_j, t_j, lb_j
            else:
                winner_id, winner_type, winner_lb = vid_j, t_j, lb_j
                loser_id, loser_type, loser_lb = vid_i, t_i, lb_i
            cleared_set.add(loser_id)
            overlaps_rows.append([sample_name, haplotype,
                                   winner_id, winner_type, winner_lb,
                                   loser_id, loser_type, loser_lb])


if len(hap0) >= 2:
    find_overlaps(hap0, 1)
if len(hap1) >= 2:
    find_overlaps(hap1, 2)

with open('~{prefix}.overlaps.tsv', 'w') as f:
    f.write('sample\thaplotype\tvariant_id_retained\tvar_type_retained\tsize_bin_retained\tvariant_id_cleared\tvar_type_cleared\tsize_bin_cleared\n')
    for row in overlaps_rows:
        f.write('\t'.join(str(x) for x in row) + '\n')

with open('~{prefix}.cleared_calls.tsv', 'w') as f:
    for vid in sorted(cleared_set):
        f.write(f'{vid}\t{sample_name}\n')

print(f'{sample_name}: {n_records} records, {len(overlaps_rows)} pairs, {len(cleared_set)} calls to clear')
CODE
    >>>

    output {
        File overlaps_tsv = "~{prefix}.overlaps.tsv"
        File cleared_calls_tsv = "~{prefix}.cleared_calls.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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


task ClearOverlapCalls {
    input {
        File vcf
        File vcf_idx
        File cleared_calls_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

to_clear = {}
with open('~{cleared_calls_tsv}') as f:
    for line in f:
        line = line.rstrip('\n')
        if not line:
            continue
        vid, sample = line.split('\t')
        if vid not in to_clear:
            to_clear[vid] = set()
        to_clear[vid].add(sample)

vcf_in = pysam.VariantFile('~{vcf}')
vcf_out = pysam.VariantFile('~{prefix}.vcf.gz', 'w', header=vcf_in.header)
n_cleared = 0

for rec in vcf_in:
    if rec.id in to_clear:
        for sample_name in to_clear[rec.id]:
            if sample_name not in rec.samples:
                continue
            s = rec.samples[sample_name]
            s['GT'] = (None, None)
            s.phased = False
            if 'GQ' in rec.format:
                s['GQ'] = None
            if 'EV' in rec.format:
                s['EV'] = '.'
            if 'BEV' in rec.format:
                s['BEV'] = '.'
            for field in ['AD', 'PL']:
                if field in rec.format:
                    try:
                        current = s.get(field)
                        if current is not None:
                            s[field] = tuple(None for _ in current)
                        else:
                            s[field] = None
                    except Exception:
                        pass
            n_cleared += 1
    vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
print(f'Cleared {n_cleared:,} sample calls')
CODE

        tabix -p vcf '~{prefix}.vcf.gz'
    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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
