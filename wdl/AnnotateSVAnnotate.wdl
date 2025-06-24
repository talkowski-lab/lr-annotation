version 1.0

workflow AnnotateSVAnnotate {
    input {
        File vcf
        File vcf_idx
        String prefix
        File coding_gtf
        File noncoding_bed
        Int? min_svlen

        Array[String] contigs

        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_preprocess
        RuntimeAttr? runtime_attr_annotate_func
    }

    scatter (contig in contigs) {
        call SubsetVcf {
            input:
                vcf = vcf,
                vcf_index = vcf_idx,
                prefix = "~{prefix}.~{contig}",
                locus = contig,
                min_svlen = min_svlen,
                runtime_attr_override = runtime_attr_subset_vcf
          }

        call PreprocessMergedVcf {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_tbi,
                prefix = "~{prefix}.~{contig}.doubled",
                runtime_attr_override = runtime_attr_preprocess
        }

        call AnnotateFunctionalConsequences {
            input:
                vcf = PreprocessMergedVcf.processed_vcf,
                vcf_idx = PreprocessMergedVcf.processed_tbi,
                noncoding_bed = noncoding_bed,
                coding_gtf = coding_gtf,
                prefix = "~{prefix}.~{contig}.doubled.anno_func",
                runtime_attr_override = runtime_attr_annotate_func
        }
    }

    call ConcatVcfs as ConcatUnannotated {
        input:
            vcfs=SubsetVcf.unannotated_vcf,
            vcfs_idx=SubsetVcf.unannotated_tbi,
            allow_overlaps=true,
            outfile_prefix="~{prefix}.unannotated.concat",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_concat
    }

    call ConcatVcfs {
        input:
            vcfs=AnnotateFunctionalConsequences.anno_vcf,
            vcfs_idx=AnnotateFunctionalConsequences.anno_tbi,
            allow_overlaps=true,
            outfile_prefix="~{prefix}.concat",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_concat
    }

    call MergeVcf {
        input:
            annotated_vcf = ConcatVcfs.concat_vcf,
            annotated_tbi = ConcatVcfs.concat_vcf_idx,
            unannotated_vcf = ConcatUnannotated.concat_vcf,
            unannotated_tbi = ConcatUnannotated.concat_vcf_idx,
            prefix = prefix,
            runtime_attr_override = runtime_attr_concat
    }

    call RevertSymbolicAlts {
        input:
            annotated_vcf = MergeVcf.merged_vcf,
            annotated_tbi = MergeVcf.merged_vcf_idx,
            original_vcf = vcf,
            original_tbi = vcf_idx,
            prefix = prefix,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File svannotated_vcf = MergeVcf.merged_vcf
        File svannotated_vcf_index = MergeVcf.merged_vcf_idx
        File svannotated_reverted_vcf = RevertSymbolicAlts.reverted_vcf
        File svannotated_reverted_vcf_index = RevertSymbolicAlts.reverted_tbi
    }
}

task SubsetVcf {
    input {
        File vcf
        File vcf_index
        String locus
        String prefix = "subset"
        Int? min_svlen

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([vcf, vcf_index], "GB")) + 2

    command <<<
        set -euxo pipefail

        VCF_FILE="~{vcf}"
        MIN_SVLEN="~{min_svlen}"

        if [ -n "$MIN_SVLEN" ]; then
            bcftools view "$VCF_FILE" --regions "~{locus}" --include "abs(INFO/SVLEN)>=$MIN_SVLEN" | bgzip > "~{prefix}.vcf.gz"
            tabix -p vcf "~{prefix}.vcf.gz"

            bcftools view "$VCF_FILE" --regions "~{locus}" --include "abs(INFO/SVLEN)<$MIN_SVLEN" | bgzip > "~{prefix}.unannotated.vcf.gz"
            tabix -p vcf "~{prefix}.unannotated.vcf.gz"
        else
            bcftools view "$VCF_FILE" --regions "~{locus}" | bgzip > "~{prefix}.vcf.gz"
            tabix -p vcf "~{prefix}.vcf.gz"
            
            touch "~{prefix}.unannotated.vcf.gz"
            touch "~{prefix}.unannotated.vcf.gz.tbi"
        fi
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_tbi = "~{prefix}.vcf.gz.tbi"
        File unannotated_vcf = "~{prefix}.unannotated.vcf.gz"
        File unannotated_tbi = "~{prefix}.unannotated.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr  runtime_default = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/ymostovoy/lr-utils-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }
}

task PreprocessMergedVcf {
    input {
        File vcf
        File vcf_idx
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(10 + size(vcf, "GB") * 2),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: "quay.io/ymostovoy/lr-process-mendelian:latest"
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        # make all DEL svlengths positive
        # make all ALTs symbolic (required for SVAnnotate)

        python /symbalts.py ~{vcf} | \
            python /abs_svlen.py /dev/stdin | \
            bcftools view -Oz > ~{prefix}.vcf.gz

        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_tbi = "~{prefix}.vcf.gz.tbi"
    }
}

task AnnotateFunctionalConsequences {
    input {
        File vcf
        File vcf_idx
        File noncoding_bed
        File coding_gtf
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 8,
        disk_gb: ceil(10 + size(vcf, "GB") * 5),
        cpu_cores: 2,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: "quay.io/ymostovoy/lr-svannotate:latest"
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)

    command <<<
        set -euxo pipefail

        gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate \
            -V ~{vcf} \
            --non-coding-bed ~{noncoding_bed} \
            --protein-coding-gtf ~{coding_gtf} \
            -O ~{prefix}.vcf
        bcftools view -Oz ~{prefix}.vcf > ~{prefix}.vcf.gz
        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File anno_vcf = "~{prefix}.vcf.gz"
        File anno_tbi = "~{prefix}.vcf.gz.tbi"
    }
}

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File]? vcfs_idx
        Boolean allow_overlaps = false
        Boolean naive = false
        Boolean sites_only = false
        Boolean sort_vcf_list = false
        String? outfile_prefix
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_override
    }

    # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
    # be held in memory or disk while working, potentially in a form that takes up more space)
    RuntimeAttr runtime_default = object {
        mem_gb: 3.75,
        disk_gb: ceil(10 + size(vcfs, "GB") * 2),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String outfile_name = outfile_prefix + ".vcf.gz"
    String allow_overlaps_flag = if allow_overlaps then "--allow-overlaps" else ""
    String naive_flag = if naive then "--naive" else ""
    String sites_only_command = if sites_only then "| bcftools view --no-version -G -Oz" else ""

    command <<<
        set -euo pipefail
        if ~{sort_vcf_list}; then
            VCFS=vcfs.list
            awk -F '/' '{print $NF"\t"$0}' ~{write_lines(vcfs)} | sort -k1,1V | awk '{print $2}' > $VCFS
        else
            VCFS="~{write_lines(vcfs)}"
        fi
        bcftools concat --no-version ~{allow_overlaps_flag} ~{naive_flag} -Oz --file-list $VCFS \
            ~{sites_only_command} > ~{outfile_name}
        tabix ~{outfile_name}
    >>>

    output {
      File concat_vcf = outfile_name
      File concat_vcf_idx = outfile_name + ".tbi"
    }
}

task MergeVcf {
    input {
        File annotated_vcf
        File annotated_tbi
        File? unannotated_vcf
        File? unannotated_tbi
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    String vcf_list = if (defined(unannotated_vcf))
        then "~{annotated_vcf} ~{unannotated_vcf}"
        else "~{annotated_vcf}"

    command <<<
        set -euxo pipefail
        bcftools concat -a -Oz -o ~{prefix}.merged.vcf.gz ~{vcf_list}
        tabix -p vcf ~{prefix}.merged.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.merged.vcf.gz"
        File merged_vcf_idx = "~{prefix}.merged.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr  runtime_default = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            ceil(10 + size(annotated_vcf, "GB")*2 + size(unannotated_vcf, "GB")*2),
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/ymostovoy/lr-utils-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }
}

task RevertSymbolicAlts {
    input {
        File annotated_vcf
        File annotated_tbi
        File original_vcf
        File original_tbi
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        
        python /revert_symbalts.py ~{annotated_vcf} ~{original_vcf} | bcftools view -Oz -o ~{prefix}.reverted.vcf.gz
        tabix -p vcf ~{prefix}.reverted.vcf.gz
    >>>

    output {
        File reverted_vcf = "~{prefix}.reverted.vcf.gz"
        File reverted_tbi = "~{prefix}.reverted.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr  runtime_default = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            ceil(10 + size(annotated_vcf, "GB")*2 + size(original_vcf, "GB")),
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/ymostovoy/lr-process-mendelian:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
