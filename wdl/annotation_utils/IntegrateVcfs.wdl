version 1.0

import "../utils/Helpers.wdl" as Helpers
import "../utils/Structs.wdl"

workflow IntegrateVcfs {
    input {
        File snv_indel_vcf
        File snv_indel_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids
        Int min_size_sv = 50
        
        String utils_docker

        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_filter_snv_indel
        RuntimeAttr? runtime_attr_filter_sv
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_annotate_svlen_svtype
    }

    call CheckSampleConsistency {
        input:
            snv_indel_vcf = snv_indel_vcf,
            sv_vcf = sv_vcf,
            sample_ids = sample_ids,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_check_samples
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToSampleList as SubsetSnvIndel {
            input:
                vcf = snv_indel_vcf,
                vcf_idx = snv_indel_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.snv_indel.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call SplitMultiallelics {
            input:
                vcf = SubsetSnvIndel.subset_vcf,
                vcf_idx = SubsetSnvIndel.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.split",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call AnnotateSvlenSvtype as AnnotateSnvIndel {
            input:
                vcf = SplitMultiallelics.split_vcf,
                vcf_idx = SplitMultiallelics.split_vcf_idx,
                prefix = "~{prefix}.~{contig}.snv_indel.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_svlen_svtype
        }

        call Helpers.SubsetVcfBySize as FilterSnvIndel {
            input:
                vcf = AnnotateSnvIndel.annotated_vcf,
                vcf_idx = AnnotateSnvIndel.annotated_vcf_idx,
                locus = contig,
                max_size = min_size_sv - 1,
                prefix = "~{prefix}.~{contig}.snv_indel",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_snv_indel
        }

        call Helpers.SubsetVcfToSampleList as SubsetSv {
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                samples = sample_ids,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call AnnotateSvlenSvtype {
            input:
                vcf = SubsetSv.subset_vcf,
                vcf_idx = SubsetSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_svlen_svtype
        }

        call Helpers.SubsetVcfBySize as FilterSv {
            input:
                vcf = AnnotateSvlenSvtype.annotated_vcf,
                vcf_idx = AnnotateSvlenSvtype.annotated_vcf_idx,
                locus = contig,
                min_size = min_size_sv,
                prefix = "~{prefix}.~{contig}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_sv
        }

        call Helpers.ConcatVcfs as MergeContigVcfs {
            input:
                vcfs = [FilterSnvIndel.subset_vcf, FilterSv.subset_vcf],
                vcfs_idx = [FilterSnvIndel.subset_vcf_idx, FilterSv.subset_vcf_idx],
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeContigVcfs.concat_vcf,
            vcfs_idx = MergeContigVcfs.concat_vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File integrated_vcf = ConcatVcfs.concat_vcf
        File integrated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task CheckSampleConsistency {
    input {
        File snv_indel_vcf
        File sv_vcf
        Array[String] sample_ids
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        printf '%s\n' ~{sep=' ' sample_ids} > requested_samples.txt

        bcftools query -l ~{snv_indel_vcf} | sort > snv_indel_samples.txt
        bcftools query -l ~{sv_vcf} | sort > sv_samples.txt

        missing_in_snv=""
        missing_in_sv=""
        
        while read sample; do
            if ! grep -q "^${sample}$" snv_indel_samples.txt; then
                missing_in_snv="${missing_in_snv} ${sample}"
            fi
            if ! grep -q "^${sample}$" sv_samples.txt; then
                missing_in_sv="${missing_in_sv} ${sample}"
            fi
        done < requested_samples.txt

        if [ -n "${missing_in_snv}" ]; then
            echo "ERROR: The following samples are missing from SNV/indel VCF:${missing_in_snv}"
            exit 1
        fi

        if [ -n "${missing_in_sv}" ]; then
            echo "ERROR: The following samples are missing from SV VCF:${missing_in_sv}"
            exit 1
        fi

        echo "All requested samples are present in both VCFs"
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([snv_indel_vcf, sv_vcf], "GB")) + 10,
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

task SplitMultiallelics {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

def parse_id(id_str):
    parts = id_str.split('_')
    alt = parts[-1]
    ref = parts[-2]
    pos = int(parts[-3])
    chrom = "_".join(parts[:-3])
    
    if len(ref) == 1 and len(alt) == 1:
        vtype = 'SNV'
    elif len(alt) > len(ref):
        vtype = 'INS'
    elif len(ref) > len(alt):
        vtype = 'DEL'
    else:
        vtype = 'COMPLEX'
        
    return chrom, pos, ref, alt, vtype

vcf_in = pysam.VariantFile("~{vcf}")
vcf_out = pysam.VariantFile("temp_unsorted.vcf", "w", header=vcf_in.header)

for rec in vcf_in:
    alts = rec.alts
    if len(alts) < 2:
        vcf_out.write(rec)
        continue

    ref_len = len(rec.ref)
    if (ref_len == 1 and all(len(a) == 1 for a in alts)) or \
       all(len(a) > ref_len for a in alts) or \
       all(len(a) < ref_len for a in alts):
        if ';' in rec.id: 
            rec.id = rec.id.replace(';', '-')
        vcf_out.write(rec)
        continue

    ids = rec.id.split(';')
    groups = {}

    for i, id_str in enumerate(ids):
        chrom, pos, ref, alt, vtype = parse_id(id_str)
        key = (chrom, pos, ref, vtype)
        if key not in groups: 
            groups[key] = []
        groups[key].append((i, alt, id_str))

    for key in sorted(groups.keys(), key=lambda k: k[1]):
        chrom, pos, ref, _ = key
        items = groups[key]
        
        new_rec = rec.copy()
        new_rec.chrom = chrom
        new_rec.pos = pos
        new_rec.ref = ref
        new_rec.alts = tuple(x[1] for x in items)
        new_rec.id = "-".join(x[2] for x in items)
        idx_map = {x[0]+1: new_i+1 for new_i, x in enumerate(items)}

        for sample in new_rec.samples:
            old_gt = new_rec.samples[sample]['GT']
            new_gt = []
            for allele in old_gt:
                if allele is None: 
                    new_gt.append(None)
                elif allele == 0: 
                    new_gt.append(0)
                elif allele in idx_map: 
                    new_gt.append(idx_map[allele])
                else: 
                    new_gt.append(0)
            new_rec.samples[sample]['GT'] = tuple(new_gt)
            
        vcf_out.write(new_rec)

vcf_out.close()
vcf_in.close()
CODE

        bcftools sort temp_unsorted.vcf -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File split_vcf = "~{prefix}.vcf.gz"
        File split_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task AnnotateSvlenSvtype {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch new_headers.txt
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=SVLEN'; then
            echo '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=SVTYPE'; then
            echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">' >> new_headers.txt
        fi

        if [ -s new_headers.txt ]; then
            bcftools annotate -h new_headers.txt ~{vcf} -Oz -o temp.vcf.gz
            tabix -p vcf temp.vcf.gz
        else
            ln -s ~{vcf} temp.vcf.gz
            ln -s ~{vcf_idx} temp.vcf.gz.tbi
        fi

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/SVTYPE\n' temp.vcf.gz | \
        awk -F'\t' '{
            split($4, alleles, ",")
            ref_len = length($3)
            max_diff = -1
            is_all_snv = 1
            best_type = "SNV"

            for (i in alleles) {
                alt_len = length(alleles[i])
                diff = (alt_len > ref_len) ? (alt_len - ref_len) : (ref_len - alt_len)
                if (diff > max_diff) {
                    max_diff = diff
                    best_type = (alt_len > ref_len) ? "INS" : (alt_len < ref_len ? "DEL" : "SNV")
                }
                if (alt_len != 1 || ref_len != 1) { is_all_snv = 0 }
            }

            if (is_all_snv == 1) {
                calc_svtype = "SNV"
                calc_svlen = 1
            } else {
                calc_svtype = best_type
                calc_svlen = max_diff
            }

            final_svtype = ($6 != ".") ? $6 : calc_svtype
            final_svlen = (is_all_snv == 1) ? 1 : (($5 != ".") ? $5 : calc_svlen)

            print $1"\t"$2"\t"$3"\t"$4"\t"final_svtype"\t"final_svlen
        }' | bgzip -c > annot.txt.gz
        
        tabix -s1 -b2 -e2 annot.txt.gz

        bcftools annotate -a annot.txt.gz \
            -c CHROM,POS,REF,ALT,INFO/SVTYPE,INFO/SVLEN \
            temp.vcf.gz -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 5,
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

