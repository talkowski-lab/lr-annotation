version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow Kanpig {
    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[String] sample_ids
        Array[File] bams
        Array[File] bais
        Array[String] sexes
        Array[String] contigs
        String prefix

        File ref_fa
        File ref_fai
        File ploidy_bed_male
        File ploidy_bed_female

        String merge_args = "--merge id"
        String kanpig_params = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"

        String kanpig_docker
        String utils_docker

        File? swap_samples

        RuntimeAttr? runtime_attr_swap_samples
        RuntimeAttr? runtime_attr_subset_to_contig
        RuntimeAttr? runtime_attr_subset_to_sample
        RuntimeAttr? runtime_attr_run_kanpig
        RuntimeAttr? runtime_attr_merge_genotypes
        RuntimeAttr? runtime_attr_concat_raw_per_sample
        RuntimeAttr? runtime_attr_concat_merged_per_sample
        RuntimeAttr? runtime_attr_merge_raw_vcfs
        RuntimeAttr? runtime_attr_merge_processed_vcfs
    }

    if (defined(swap_samples)) {
        call Helpers.SwapSampleIds {
            input:
                vcf = cohort_vcf,
                vcf_idx = cohort_vcf_idx,
                sample_swap_list = select_first([swap_samples]),
                prefix = "~{prefix}.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples
        }
    }

    File final_cohort_vcf = select_first([SwapSampleIds.swapped_vcf, cohort_vcf])
    File final_cohort_vcf_idx = select_first([SwapSampleIds.swapped_vcf_idx, cohort_vcf_idx])

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = final_cohort_vcf,
                vcf_idx = final_cohort_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_to_contig
        }
    }

    scatter (i in range(length(sample_ids))) {
        File ploidy_bed = if sexes[i] == "M" then ploidy_bed_male else ploidy_bed_female

        scatter (j in range(length(contigs))) {
            call Helpers.SubsetVcfToSamples {
                input:
                    vcf = SubsetVcfToContig.subset_vcf[j],
                    vcf_idx = SubsetVcfToContig.subset_vcf_idx[j],
                    samples = [sample_ids[i]],
                    filter_to_sample = false,
                    prefix = "~{prefix}.~{sample_ids[i]}.~{contigs[j]}.subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_to_sample
            }

            call RunKanpig {
                input:
                    input_vcf = SubsetVcfToSamples.subset_vcf,
                    input_vcf_idx = SubsetVcfToSamples.subset_vcf_idx,
                    bam = bams[i],
                    bai = bais[i],
                    sample_id = sample_ids[i],
                    ploidy_bed = ploidy_bed,
                    kanpig_params = kanpig_params,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}.~{sample_ids[i]}.~{contigs[j]}.kanpig",
                    docker = kanpig_docker,
                    runtime_attr_override = runtime_attr_run_kanpig
            }

            call MergeGenotypes {
                input:
                    base_vcf = SubsetVcfToSamples.subset_vcf,
                    base_vcf_idx = SubsetVcfToSamples.subset_vcf_idx,
                    kanpig_vcf = RunKanpig.regenotyped_vcf,
                    kanpig_vcf_idx = RunKanpig.regenotyped_vcf_idx,
                    prefix = "~{prefix}.~{sample_ids[i]}.~{contigs[j]}.kanpig_merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_genotypes
            }
        }

        call Helpers.ConcatVcfs as ConcatRawPerSample {
            input:
                vcfs = RunKanpig.regenotyped_vcf,
                vcf_idxs = RunKanpig.regenotyped_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.~{sample_ids[i]}.kanpig_raw",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_raw_per_sample
        }

        call Helpers.ConcatVcfs as ConcatMergedPerSample {
            input:
                vcfs = MergeGenotypes.merged_vcf,
                vcf_idxs = MergeGenotypes.merged_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.~{sample_ids[i]}.kanpig_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_merged_per_sample
        }
    }

    call Helpers.MergeVcfs as MergeRaw {
        input:
            vcfs = ConcatRawPerSample.concat_vcf,
            vcf_idxs = ConcatRawPerSample.concat_vcf_idx,
            prefix = "~{prefix}.kanpig_raw",
            extra_args = merge_args,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_raw_vcfs
    }

    call Helpers.MergeVcfs as MergeProcessed {
        input:
            vcfs = ConcatMergedPerSample.concat_vcf,
            vcf_idxs = ConcatMergedPerSample.concat_vcf_idx,
            prefix = "~{prefix}.kanpig_merged",
            extra_args = merge_args,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_processed_vcfs
    }

    output {
        File kanpig_raw_vcf = MergeRaw.merged_vcf
        File kanpig_raw_vcf_idx = MergeRaw.merged_vcf_idx
        File kanpig_regenotyped_vcf = MergeProcessed.merged_vcf
        File kanpig_regenotyped_vcf_idx = MergeProcessed.merged_vcf_idx
    }
}

task RunKanpig {
    input {
        File input_vcf
        File input_vcf_idx
        File bam
        File bai
        String sample_id
        File ploidy_bed
        String kanpig_params
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: { localization_optional: true }
        bai: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt

        kanpig gt \
            --input ~{input_vcf} \
            --out ~{prefix}.kanpig.vcf \
            --reads ~{bam} \
            --reference ~{ref_fa} \
            --ploidy-bed ~{ploidy_bed} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --sample ~{sample_id} \
            ~{kanpig_params}

        bcftools sort \
            -Oz -o ~{prefix}.vcf.gz \
            ~{prefix}.kanpig.vcf
        
        rm ~{prefix}.kanpig.vcf

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File regenotyped_vcf = "~{prefix}.vcf.gz"
        File regenotyped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(input_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: 2
        memory: 4 + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeGenotypes {
    input {
        File base_vcf
        File base_vcf_idx
        File kanpig_vcf
        File kanpig_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

kp_vcf = pysam.VariantFile("~{kanpig_vcf}")
base_vcf = pysam.VariantFile("~{base_vcf}")

header = base_vcf.header.copy()
if 'AD' not in header.formats:
    header.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">')
if 'GQ' not in header.formats:
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
if 'DP' not in header.formats:
    header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">')

out = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)
sample = list(base_vcf.header.samples)[0]

def is_non_ref(gt):
    return any(a is not None and a > 0 for a in gt)

def is_missing(gt):
    return all(a is None for a in gt)

for rec in base_vcf:
    rec.translate(out.header)
    base_gt = rec.samples[sample]['GT']
    kp_rec = next(
        r for r in kp_vcf.fetch(rec.chrom, rec.pos - 1, rec.pos)
        if r.ref == rec.ref and r.alts == rec.alts and r.id == rec.id
    )
    kp_gt = kp_rec.samples[sample]['GT']

    # Ref/missing in base and ref in kanpig → set FORMAT fields to kanpig
    if not is_non_ref(base_gt) and not is_non_ref(kp_gt) and not is_missing(kp_gt):
        for field in ['GT', 'AD', 'GQ', 'DP']:
            val = kp_rec.samples[sample][field]
            if field == 'AD' and (val is None or len(val) != len(rec.alts) + 1):
                continue
            rec.samples[sample][field] = val
    
    out.write(rec)

base_vcf.close()
kp_vcf.close()
out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(base_vcf, "GB") + size(kanpig_vcf, "GB")) + 10,
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
