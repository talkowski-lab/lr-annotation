version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateDbVaR {
    input {
        File vcf
        File vcf_idx
        File dbvar_vcf
        File dbvar_vcf_idx
        Array[String] contigs
        String prefix

        Int min_length

        Int? records_per_shard

        Float del_size_similarity = 0.8
        Float del_reciprocal_overlap = 0.8
        Int del_breakpoint_window = 500
        Float dup_size_similarity = 0.8
        Float dup_reciprocal_overlap = 0.8
        Int dup_breakpoint_window = 500
        Float ins_size_similarity = 0.8
        Int ins_breakpoint_window = 100

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_dbvar_vcf
        RuntimeAttr? runtime_attr_convert_symbolic
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByArgs as FilterVcf {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                include_args = "abs(INFO/allele_length) >= ~{min_length} && (INFO/allele_type = \"del\" || INFO/allele_type = \"ins\" || INFO/allele_type = \"dup\")",
                extra_args = if single_contig then "" else "--regions " + contig,
                prefix = "~{prefix}.~{contig}.filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.ConvertToSymbolic {
            input:
                vcf = FilterVcf.subset_vcf,
                vcf_idx = FilterVcf.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.symbolic",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert_symbolic
        }

        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetDbVaRVcf {
                input:
                    vcf = dbvar_vcf,
                    vcf_idx = dbvar_vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.dbvar",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_dbvar_vcf
            }
        }

        File contig_dbvar_vcf = select_first([SubsetDbVaRVcf.subset_vcf, dbvar_vcf])
        File contig_dbvar_vcf_idx = select_first([SubsetDbVaRVcf.subset_vcf_idx, dbvar_vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = ConvertToSymbolic.processed_vcf,
                    vcf_idx = ConvertToSymbolic.processed_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [ConvertToSymbolic.processed_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [ConvertToSymbolic.processed_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call AnnotateDbVaRIds {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    dbvar_vcf = contig_dbvar_vcf,
                    dbvar_vcf_idx = contig_dbvar_vcf_idx,
                    del_size_similarity = del_size_similarity,
                    del_reciprocal_overlap = del_reciprocal_overlap,
                    del_breakpoint_window = del_breakpoint_window,
                    dup_size_similarity = dup_size_similarity,
                    dup_reciprocal_overlap = dup_reciprocal_overlap,
                    dup_breakpoint_window = dup_breakpoint_window,
                    ins_size_similarity = ins_size_similarity,
                    ins_breakpoint_window = ins_breakpoint_window,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.dbvar_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = AnnotateDbVaRIds.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.dbvar_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, AnnotateDbVaRIds.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotations {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.dbvar_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_dbvar = select_first([MergeAnnotations.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task AnnotateDbVaRIds {
    input {
        File vcf
        File vcf_idx
        File dbvar_vcf
        File dbvar_vcf_idx
        Float del_size_similarity
        Float del_reciprocal_overlap
        Int del_breakpoint_window
        Float dup_size_similarity
        Float dup_reciprocal_overlap
        Int dup_breakpoint_window
        Float ins_size_similarity
        Int ins_breakpoint_window
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import pysam

DEL_SIZE_SIM = float(~{del_size_similarity})
DEL_REC_OVL = float(~{del_reciprocal_overlap})
DEL_BP_WIN = int(~{del_breakpoint_window})
DUP_SIZE_SIM = float(~{dup_size_similarity})
DUP_REC_OVL = float(~{dup_reciprocal_overlap})
DUP_BP_WIN = int(~{dup_breakpoint_window})
INS_SIZE_SIM = float(~{ins_size_similarity})
INS_BP_WIN = int(~{ins_breakpoint_window})

PARAMS = {
    'DEL': {'allowed': {'DEL', 'CNV'}, 'size_sim': DEL_SIZE_SIM, 'rec_ovl': DEL_REC_OVL, 'bp_win': DEL_BP_WIN},
    'DUP': {'allowed': {'DUP', 'CNV'}, 'size_sim': DUP_SIZE_SIM, 'rec_ovl': DUP_REC_OVL, 'bp_win': DUP_BP_WIN},
    'INS': {'allowed': {'INS'}, 'size_sim': INS_SIZE_SIM, 'rec_ovl': None, 'bp_win': INS_BP_WIN},
}

def parse_svlen(rec):
    raw = rec.info.get('SVLEN', None)
    if raw is not None:
        if isinstance(raw, (list, tuple)):
            raw = raw[0]
        try:
            return abs(int(raw))
        except (ValueError, TypeError):
            pass
    return max(0, rec.stop - rec.start)

def passes_del_dup(qs, qe, ql, cs, ce, cl, size_sim, rec_ovl, bp_win):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < size_sim:
        return False
    overlap = max(0, min(qe, ce) - max(qs, cs))
    if overlap / min(ql, cl) < rec_ovl:
        return False
    if abs(qs - cs) > bp_win and abs(qe - ce) > bp_win:
        return False
    return True

def passes_ins(qs, ql, cs, cl, size_sim, bp_win):
    if ql == 0 or cl == 0:
        return False
    if min(ql, cl) / max(ql, cl) < size_sim:
        return False
    if abs(qs - cs) > bp_win:
        return False
    return True

dbvar = pysam.VariantFile("~{dbvar_vcf}")
query = pysam.VariantFile("~{vcf}")

with open("~{prefix}.annotations.tsv", 'w') as out:
    for rec in query:
        svtype = rec.info.get('SVTYPE')
        if svtype not in PARAMS:
            continue
        params = PARAMS[svtype]
        chrom = rec.chrom
        qs = rec.start
        qe = rec.stop
        ql = parse_svlen(rec)
        ref = rec.ref if rec.ref else '.'
        alt = rec.alts[0] if rec.alts else '.'
        var_id = rec.id if rec.id else '.'

        bp_win = params['bp_win']
        size_sim = params['size_sim']
        if svtype == 'INS':
            region_start = max(0, qs - bp_win)
            region_end = qs + bp_win + 1
        else:
            margin = bp_win + int(ql / size_sim) + 1
            region_start = max(0, qs - margin)
            region_end = qe + margin

        try:
            cands = dbvar.fetch(chrom, region_start, region_end)
        except (ValueError, KeyError):
            continue

        matched = []
        for cand in cands:
            ctype = cand.info.get('SVTYPE')
            if ctype not in params['allowed']:
                continue
            cs = cand.start
            ce = cand.stop
            cl = parse_svlen(cand)
            cand_id = cand.id if cand.id else cand.info.get('DBVARID', '.')
            if isinstance(cand_id, (list, tuple)):
                cand_id = cand_id[0]

            if svtype == 'INS':
                ok = passes_ins(qs, ql, cs, cl, size_sim, bp_win)
            else:
                ok = passes_del_dup(qs, qe, ql, cs, ce, cl, size_sim, params['rec_ovl'], bp_win)
            if ok:
                matched.append(str(cand_id))

        if matched:
            seen = set()
            uniq = [m for m in matched if not (m in seen or seen.add(m))]
            out.write(f"{chrom}\t{qs + 1}\t{ref}\t{alt}\t{var_id}\t{','.join(uniq)}\n")

query.close()
dbvar.close()
EOF
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, dbvar_vcf], "GB")) + 10,
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
