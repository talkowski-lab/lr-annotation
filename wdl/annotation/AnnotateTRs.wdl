version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        Array[String] contigs
        String prefix

        Array[String] sample_ids

        Array[File] tr_catalogs
        Array[String] tr_catalog_ids

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig_base
        RuntimeAttr? runtime_attr_subset_contig_tr
        RuntimeAttr? runtime_attr_subset_samples_base
        RuntimeAttr? runtime_attr_subset_samples_tr
        RuntimeAttr? runtime_attr_check_sample_consistency
        RuntimeAttr? runtime_attr_set_missing_filters
        RuntimeAttr? runtime_attr_subset_catalog
        RuntimeAttr? runtime_attr_tag_tr_vcf
        RuntimeAttr? runtime_attr_set_tr_ids
        RuntimeAttr? runtime_attr_annotate_vcf
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetContigBase {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                extra_args = "--min-ac 1",
                prefix = "~{prefix}.~{contig}.integrated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_base
        }

        call Helpers.SubsetVcfToSamples as SubsetSamplesBase {
            input:
                vcf = SubsetContigBase.subset_vcf,
                vcf_idx = SubsetContigBase.subset_vcf_idx,
                samples = sample_ids,
                prefix = "~{prefix}.~{contig}.integrated_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_samples_base
        }

        call Helpers.SubsetVcfToContig as SubsetContigTr {
            input:
                vcf = tr_vcf,
                vcf_idx = tr_vcf_idx,
                contig = contig,
                extra_args = "--min-ac 1",
                prefix = "~{prefix}.~{contig}.tr",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_tr
        }

        call Helpers.SubsetVcfToSamples as SubsetSamplesTr {
            input:
                vcf = SubsetContigTr.subset_vcf,
                vcf_idx = SubsetContigTr.subset_vcf_idx,
                samples = sample_ids,
                prefix = "~{prefix}.~{contig}.tr_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_samples_tr
        }

        call Helpers.CheckSampleConsistency {
            input:
                vcfs = [SubsetSamplesBase.subset_vcf, SubsetSamplesTr.subset_vcf],
                vcf_idxs = [SubsetSamplesBase.subset_vcf_idx, SubsetSamplesTr.subset_vcf_idx],
                sample_ids = sample_ids,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_check_sample_consistency
        }

        call Helpers.SetMissingFiltersToPass {
            input:
                vcf = SubsetSamplesTr.subset_vcf,
                vcf_idx = SubsetSamplesTr.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_subset.refiltered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_missing_filters
        }

        scatter (i in range(length(tr_catalogs))) {
            call Helpers.SubsetTsvToContig {
                input:
                    tsv = tr_catalogs[i],
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.catalog~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_catalog
            }
        }

        call TagTrVcfWithCatalogs {
            input:
                vcf = SetMissingFiltersToPass.filtered_vcf,
                vcf_idx = SetMissingFiltersToPass.filtered_vcf_idx,
                catalogs = SubsetTsvToContig.subset_tsv,
                catalog_ids = tr_catalog_ids,
                prefix = "~{prefix}.~{contig}.tr_subset.tagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_tag_tr_vcf
        }

        call SetTrVariantIds {
            input:
                vcf = TagTrVcfWithCatalogs.tagged_vcf,
                vcf_idx = TagTrVcfWithCatalogs.tagged_vcf_idx,
                prefix = "~{prefix}.~{contig}.tr_subset.ids",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_set_tr_ids
        }

        call AnnotateVcfWithTRs {
            input:
                vcf = SubsetSamplesBase.subset_vcf,
                vcf_idx = SubsetSamplesBase.subset_vcf_idx,
                tr_vcf = SetTrVariantIds.renamed_vcf,
                tr_vcf_idx = SetTrVariantIds.renamed_vcf_idx,
                prefix = "~{prefix}.~{contig}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_annotate_vcf
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = AnnotateVcfWithTRs.annotated_vcf,
            vcf_idxs = AnnotateVcfWithTRs.annotated_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.annotated_trs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File tr_annotated_vcf = ConcatVcfs.concat_vcf
        File tr_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task TagTrVcfWithCatalogs {
    input {
        File vcf
        File vcf_idx
        Array[File] catalogs
        Array[String] catalog_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pysam import VariantFile

catalog_files = "~{sep=',' catalogs}".split(',')
catalog_ids_list = "~{sep=',' catalog_ids}".split(',')

lookup = {}
for cat_file, cat_id in zip(catalog_files, catalog_ids_list):
    with open(cat_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            trid = None
            for field in parts[3].split(';'):
                if field.startswith('ID='):
                    trid = field[3:]
                    break
            if trid:
                lookup[trid] = cat_id

vcf_in = VariantFile("~{vcf}")
header = vcf_in.header.copy()
if 'allele_type' not in header.info:
    header.add_meta('INFO', items=[('ID', 'allele_type'), ('Number', 1), ('Type', 'String'), ('Description', 'Allele type')])
if 'SOURCE' not in header.info:
    header.add_meta('INFO', items=[('ID', 'SOURCE'), ('Number', 1), ('Type', 'String'), ('Description', 'Source of variant call')])

vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=header)
for record in vcf_in:
    record.info['allele_type'] = 'trv'
    record.info['SOURCE'] = lookup.get(record.info.get('TRID'))
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File tagged_vcf = "~{prefix}.vcf.gz"
        File tagged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 5,
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

task SetTrVariantIds {
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
from pysam import VariantFile
from collections import defaultdict

vcf_in = VariantFile("~{vcf}")
id_counts = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)}"
    id_counts[new_id] += 1
vcf_in.close()

vcf_in = VariantFile("~{vcf}")
vcf_out = VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
id_seen = defaultdict(int)
for record in vcf_in:
    new_id = f"{record.chrom}-{record.pos}-TRV-{len(record.ref)}"
    if id_counts[new_id] > 1:
        id_seen[new_id] += 1
        record.id = f"{new_id}_{id_seen[new_id]}"
    else:
        record.id = new_id
    vcf_out.write(record)
vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
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

task AnnotateVcfWithTRs {
    input {
        File vcf
        File vcf_idx
        File tr_vcf
        File tr_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            ~{tr_vcf} \
            | awk 'BEGIN{OFS="\t"} {
                end = $2 + length($4)
                print $1, $2, end, $3, $4, $5
            }' > tr.bed
        
        bcftools query \
            -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
            ~{vcf} \
            | awk 'BEGIN{OFS="\t"} {
                end = $2 + length($4)
                print $1, $2, end, $3, $4, $5
            }' > vcf.bed

        bedtools intersect \
            -f 1.0 \
            -wa \
            -wb \
            -a vcf.bed \
            -b tr.bed \
            > overlaps.bed

        rm -f tr.bed vcf.bed

        awk 'BEGIN{OFS="\t"} {
            print $1, $2, $5, $6, $4, "1", $10
        }' overlaps.bed \
            | sort -k1,1 -k2,2n \
            | bgzip -c > annotations.tsv.gz

        tabix -s 1 -b 2 -e 2 annotations.tsv.gz
        
        rm -f overlaps.bed

        cat <<EOF > new_header.txt
##INFO=<ID=TR_OVERLAPPED,Number=0,Type=Flag,Description="Variant enveloped by tandem repeat">
##INFO=<ID=TRID,Number=1,Type=String,Description="ID of enveloping tandem repeat">
EOF

        bcftools annotate \
            -a annotations.tsv.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/TR_OVERLAPPED,INFO/TRID \
            -h new_header.txt \
            -Oz -o vcf_annotated.vcf.gz \
            ~{vcf}

        tabix -p vcf vcf_annotated.vcf.gz
        
        rm -f annotations.tsv.gz annotations.tsv.gz.tbi new_header.txt

        bcftools concat \
            --allow-overlaps \
            -Oz -o ~{prefix}.vcf.gz \
            vcf_annotated.vcf.gz \
            ~{tr_vcf}

        tabix -p vcf ~{prefix}.vcf.gz
        
        rm -f vcf_annotated.vcf.gz vcf_annotated.vcf.gz.tbi
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10 * ceil(size(vcf, "GB") + size(tr_vcf, "GB")) + 20,
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
