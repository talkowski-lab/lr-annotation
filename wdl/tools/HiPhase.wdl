version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow HiPhase {
    input {
        File bam
        File bai
        File small_vcf
        File small_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        File? trgt_vcf
        File? trgt_vcf_idx
        Array[String] contigs
        String prefix

        String? hiphase_extra_args

        File ref_fa
        File ref_fai

        String sv_base_mini_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_preprocess_vcf
        RuntimeAttr? runtime_attr_sync_contigs
        RuntimeAttr? runtime_attr_hiphase
        RuntimeAttr? runtime_attr_hiphase_trgt
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call PreprocessVCF { 
            input:
                vcf = sv_vcf,
                vcf_idx = sv_vcf_idx,
                prefix = "~{prefix}.~{contig}.preprocessed",
                runtime_attr_override = runtime_attr_preprocess_vcf
        }

        call Helpers.SubsetVcfToContig as SubsetVcfShort { 
            input:
                vcf = small_vcf,
                vcf_idx = small_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.short",
                docker = sv_pipeline_base_docker
        }

        call Helpers.SubsetVcfToContig as SubsetVcfSv {
            input:
                vcf = PreprocessVCF.preprocessed_vcf,
                vcf_idx = PreprocessVCF.preprocessed_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv",
                docker = sv_pipeline_base_docker
        }

        call SyncContigs {
            input:
                small_vcf = SubsetVcfShort.subset_vcf,
                small_vcf_idx = SubsetVcfShort.subset_vcf_idx,
                sv_vcf = SubsetVcfSv.subset_vcf,
                sv_vcf_idx = SubsetVcfSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.synced",
                docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_sync_contigs
        }

        if (defined(trgt_vcf) && defined(trgt_vcf_idx)) {
            call Helpers.SubsetVcfToContig as SubsetVcfTRGT { 
                input:
                    vcf = select_first([trgt_vcf]),
                    vcf_idx = select_first([trgt_vcf_idx]),
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.trgt",
                    docker = sv_pipeline_base_docker
            }

            call HiPhaseTRGT { 
                input:
                    bam = bam,
                    bai = bai,
                    unphased_sv_vcf = SyncContigs.synced_vcf,
                    unphased_sv_idx = SyncContigs.synced_idx,
                    unphased_snp_vcf = SubsetVcfShort.subset_vcf,
                    unphased_snp_idx = SubsetVcfShort.subset_vcf_idx,
                    unphased_trgt_vcf = SubsetVcfTRGT.subset_vcf,
                    unphased_trgt_idx = SubsetVcfTRGT.subset_vcf_idx,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = prefix,
                    extra_args = hiphase_extra_args,
                    runtime_attr_override = runtime_attr_hiphase_trgt
            }

            call Helpers.ConcatVcfsLR as ConcatWithTRGT {
                input:
                    vcfs = [HiPhaseTRGT.phased_snp_vcf, HiPhaseTRGT.phased_sv_vcf, HiPhaseTRGT.phased_trgt_vcf],
                    vcf_idxs = [HiPhaseTRGT.phased_snp_vcf_idx, HiPhaseTRGT.phased_sv_vcf_idx, HiPhaseTRGT.phased_trgt_vcf_idx],
                    prefix = "~{prefix}.~{contig}.hiphased",
                    sv_base_mini_docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_concat
            }
        }

        if (!defined(trgt_vcf)) {
            call HiPhase { 
                input:
                    bam = bam,
                    bai = bai,
                    unphased_sv_vcf = SyncContigs.synced_vcf,
                    unphased_sv_idx = SyncContigs.synced_idx,
                    unphased_snp_vcf = SubsetVcfShort.subset_vcf,
                    unphased_snp_idx = SubsetVcfShort.subset_vcf_idx,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = prefix,
                    extra_args = hiphase_extra_args,
                    runtime_attr_override = runtime_attr_hiphase
            }

            call Helpers.ConcatVcfsLR as ConcatWithoutTRGT {
                input:
                    vcfs = [HiPhase.phased_snp_vcf, HiPhase.phased_sv_vcf],
                    vcf_idxs = [HiPhase.phased_snp_vcf_idx, HiPhase.phased_sv_vcf_idx],
                    prefix = "~{prefix}.~{contig}.hiphased",
                    sv_base_mini_docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_concat
            }
        }
    }

    call Helpers.ConcatVcfsLR {
        input:
            vcfs = select_all(flatten([ConcatWithTRGT.concat_vcf, ConcatWithoutTRGT.concat_vcf])),
            vcf_idxs = select_all(flatten([ConcatWithTRGT.concat_vcf_idx, ConcatWithoutTRGT.concat_vcf_idx])),
            prefix = "~{prefix}.phased",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File hiphase_vcf = ConcatVcfsLR.concat_vcf
        File hiphase_vcf_idx = ConcatVcfsLR.concat_vcf_idx
        Array[File] hiphase_haplotag_files = select_all(flatten([HiPhaseTRGT.haplotag_file, HiPhase.haplotag_file]))
        Array[File] hiphase_stats = select_all(flatten([HiPhaseTRGT.hiphase_stats, HiPhase.hiphase_stats]))
        Array[File] hiphase_blocks = select_all(flatten([HiPhaseTRGT.hiphase_blocks, HiPhase.hiphase_blocks]))
        Array[File] hiphase_summary = select_all(flatten([HiPhaseTRGT.hiphase_summary, HiPhase.hiphase_summary]))
    }
}

task PreprocessVCF {
    input {
        File vcf
        File vcf_idx
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    String docker_dir = "/truvari_intrasample"

    command <<<
        set -euo pipefail

        cp ~{docker_dir}/convert_lower_case.py .

        python convert_lower_case.py -i ~{vcf} -o ~{prefix}.uppercase.vcf
        bgzip ~{prefix}.uppercase.vcf ~{prefix}.uppercase.vcf.gz
        tabix -p vcf ~{prefix}.uppercase.vcf.gz

        bcftools +setGT ~{prefix}.uppercase.vcf.gz --no-version -Oz -o ~{prefix}.vcf.gz -- --target-gt a --new-gt u
        bcftools index -t ~{prefix}.vcf.gz

        ls -l
    >>>

    output {
        File preprocessed_vcf = "~{prefix}.vcf.gz"
        File preprocessed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
        docker: "hangsuunc/cleanvcf:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task SyncContigs {
    input {
        File small_vcf
        File small_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view --header-only ~{sv_vcf} | grep -v "^##contig" > sv_header_no_contigs.txt
        bcftools view --header-only ~{small_vcf} | grep "^##contig" > small_vcf_contigs.txt

        grep "^##fileformat" sv_header_no_contigs.txt > new_header.txt
        cat small_vcf_contigs.txt >> new_header.txt
        grep -v "^##fileformat" sv_header_no_contigs.txt >> new_header.txt

        bcftools reheader --header new_header.txt --output ~{prefix}.vcf.gz ~{sv_vcf}
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File synced_vcf = "~{prefix}.vcf.gz"
        File synced_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(small_vcf, "GB") + size(sv_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0,
        docker: docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task HiPhase {
    input {
        File bam
        File bai
        File unphased_snp_vcf
        File unphased_snp_idx
        File unphased_sv_vcf
        File unphased_sv_idx
        String? extra_args
        File ref_fa
        File ref_fai
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 6,
        disk_gb: 2 * ceil(size(bam, "GB")) + 50,
        boot_disk_gb: 100,
        preemptible_tries: 2,
        max_retries: 0,
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/hiphase:v1.5.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        hiphase \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --bam ~{bam} \
            --reference ~{ref_fa} \
            --global-realignment-cputime 300 \
            --vcf ~{unphased_snp_vcf} \
            --output-vcf ~{prefix}_phased_snp.vcf.gz \
            --vcf ~{unphased_sv_vcf} \
            --output-vcf ~{prefix}_phased_sv.vcf.gz \
            --haplotag-file ~{prefix}_phased_sv_haplotag.tsv \
            --stats-file ~{prefix}.stats.csv \
            --blocks-file ~{prefix}.blocks.tsv \
            --summary-file ~{prefix}.summary.tsv \
            --verbose \
            ~{if defined(extra_args) then extra_args else ""}
    >>>

    output {
        File phased_snp_vcf = "~{prefix}_phased_snp.vcf.gz"
        File phased_snp_vcf_idx = "~{prefix}_phased_snp.vcf.gz.tbi"
        File phased_sv_vcf   = "~{prefix}_phased_sv.vcf.gz"
        File phased_sv_vcf_idx = "~{prefix}_phased_sv.vcf.gz.tbi"
        File haplotag_file = "~{prefix}_phased_sv_haplotag.tsv"
        File hiphase_stats = "~{prefix}.stats.csv"
        File hiphase_blocks = "~{prefix}.blocks.tsv"
        File hiphase_summary = "~{prefix}.summary.tsv"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: 12 + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}

task HiPhaseTRGT {
    input {
        File bam
        File bai
        File unphased_snp_vcf
        File unphased_snp_idx
        File unphased_sv_vcf
        File unphased_sv_idx
        File unphased_trgt_vcf
        File unphased_trgt_idx
        String? extra_args
        File ref_fa
        File ref_fai
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 6,
        disk_gb: 2 * ceil(size(bam, "GB")) + 50,
        boot_disk_gb: 100,
        preemptible_tries: 2,
        max_retries: 0,
        docker: "us.gcr.io/broad-dsp-lrma/hangsuunc/hiphase:v1.5.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        hiphase \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --bam ~{bam} \
            --reference ~{ref_fa} \
            --global-realignment-cputime 300 \
            --vcf ~{unphased_snp_vcf} \
            --output-vcf ~{prefix}_phased_snp.vcf.gz \
            --vcf ~{unphased_sv_vcf} \
            --output-vcf ~{prefix}_phased_sv.vcf.gz \
            --vcf ~{unphased_trgt_vcf} \
            --output-vcf ~{prefix}_phased_trgt.vcf.gz \
            --haplotag-file ~{prefix}_phased_sv_haplotag.tsv \
            --stats-file ~{prefix}.stats.csv \
            --blocks-file ~{prefix}.blocks.tsv \
            --summary-file ~{prefix}.summary.tsv \
            --verbose \
            ~{if defined(extra_args) then extra_args else ""}
    >>>

    output {
        File phased_snp_vcf = "~{prefix}_phased_snp.vcf.gz"
        File phased_snp_vcf_idx = "~{prefix}_phased_snp.vcf.gz.tbi"
        File phased_sv_vcf   = "~{prefix}_phased_sv.vcf.gz"
        File phased_sv_vcf_idx = "~{prefix}_phased_sv.vcf.gz.tbi"
        File phased_trgt_vcf   = "~{prefix}_phased_trgt.vcf.gz"
        File phased_trgt_vcf_idx = "~{prefix}_phased_trgt.vcf.gz.tbi"
        File haplotag_file = "~{prefix}_phased_sv_haplotag.tsv"
        File hiphase_stats = "~{prefix}.stats.csv"
        File hiphase_blocks = "~{prefix}.blocks.tsv"
        File hiphase_summary = "~{prefix}.summary.tsv"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: 12 + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: select_first([runtime_attr.docker, default_attr.docker])
    }
}
