version 1.0

import "general/Structs.wdl"

workflow AnnotateSVAnnotate {
    input {
        File vcf
        File vcf_idx
        String prefix
        File coding_gtf
        File noncoding_bed
        Int min_svlen

        Array[String] contigs

        String pipeline_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_preprocess
        RuntimeAttr? runtime_attr_annotate_func
        RuntimeAttr? runtime_attr_concat_unannotated
        RuntimeAttr? runtime_attr_concat_annotated
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_postprocess
    }

    scatter (contig in contigs) {
        call SubsetVcf {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                prefix = "~{prefix}.~{contig}",
                locus = contig,
                min_svlen = min_svlen,
                runtime_attr_override = runtime_attr_subset_vcf
          }

        call PreprocessVcf {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_tbi,
                prefix = "~{prefix}.~{contig}.doubled",
                runtime_attr_override = runtime_attr_preprocess
        }

        call AnnotateFunctionalConsequences {
            input:
                vcf = PreprocessVcf.processed_vcf,
                vcf_idx = PreprocessVcf.processed_tbi,
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
            pipeline_docker=pipeline_docker,
            runtime_attr_override=runtime_attr_concat_unannotated
    }

    call ConcatVcfs as ConcatAnnotated {
        input:
            vcfs=AnnotateFunctionalConsequences.anno_vcf,
            vcfs_idx=AnnotateFunctionalConsequences.anno_tbi,
            allow_overlaps=true,
            outfile_prefix="~{prefix}.concat",
            pipeline_docker=pipeline_docker,
            runtime_attr_override=runtime_attr_concat_annotated
    }

    call MergeVcf {
        input:
            annotated_vcf = ConcatAnnotated.concat_vcf,
            annotated_tbi = ConcatAnnotated.concat_vcf_index,
            unannotated_vcf = ConcatUnannotated.concat_vcf,
            unannotated_tbi = ConcatUnannotated.concat_vcf_index,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_merge
    }

    call PostprocessVcf {
        input:
            annotated_vcf = MergeVcf.merged_vcf,
            annotated_tbi = MergeVcf.merged_vcf_index,
            original_vcf = vcf,
            original_tbi = vcf_idx,
            prefix = prefix,
            pipeline_docker = pipeline_docker,
            runtime_attr_override = runtime_attr_postprocess
    }

    output {
        File sv_annotated_vcf = PostprocessVcf.reverted_vcf
        File sv_annotated_vcf_index = PostprocessVcf.reverted_tbi
    }
}

task SubsetVcf {
    input {
        File vcf
        File vcf_idx
        String locus
        String prefix = "subset"
        Int min_svlen

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} --regions "~{locus}" --include "abs(INFO/SVLEN)>=~{min_svlen}" | bgzip > "~{prefix}.vcf.gz"
        tabix -p vcf "~{prefix}.vcf.gz"

        bcftools view ~{vcf} --regions "~{locus}" --include "abs(INFO/SVLEN)<~{min_svlen}" | bgzip > "~{prefix}.unannotated.vcf.gz"
        tabix -p vcf "~{prefix}.unannotated.vcf.gz"
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_tbi = "~{prefix}.vcf.gz.tbi"
        File unannotated_vcf = "~{prefix}.unannotated.vcf.gz"
        File unannotated_tbi = "~{prefix}.unannotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 4*ceil(size([vcf, vcf_idx], "GB")) + 2,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-utils-basic:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task PreprocessVcf {
    input {
        File vcf
        File vcf_idx
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # make all DEL svlengths positive
        # make all ALTs symbolic (required for SVAnnotate)

        python /opt/gnomad-lr/scripts/helpers/symbalts.py ~{vcf} | \
            python /opt/gnomad-lr/scripts/helpers/abs_svlen.py /dev/stdin | \
            bcftools view -Oz > ~{prefix}.vcf.gz

        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_tbi = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + size(vcf, "GB") * 2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-process-mendelian:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(10 + size(vcf, "GB") * 5),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int java_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 1000 * 0.7)

    command <<<
        set -euo pipefail

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

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-svannotate:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        String pipeline_docker

        RuntimeAttr? runtime_attr_override
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
      File concat_vcf_index = outfile_name + ".tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + size(vcfs, "GB") * 2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeVcf {
    input {
        File annotated_vcf
        File annotated_tbi
        File unannotated_vcf
        File unannotated_tbi
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        echo "~{annotated_vcf}" > vcf.list
        echo "~{unannotated_vcf}" >> vcf.list

        bcftools concat -a -f vcf.list -Oz -o ~{prefix}.merged.vcf.gz
        tabix -p vcf ~{prefix}.merged.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.merged.vcf.gz"
        File merged_vcf_index = "~{prefix}.merged.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + size(annotated_vcf, "GB")*2 + size(unannotated_vcf, "GB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task PostprocessVcf {
    input {
        File annotated_vcf
        File annotated_tbi
        File original_vcf
        File original_tbi
        String prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        python /opt/gnomad-lr/scripts/helpers/revert_symbalts.py ~{annotated_vcf} ~{original_vcf} | bcftools view -Oz -o ~{prefix}.reverted.vcf.gz
        tabix -p vcf ~{prefix}.reverted.vcf.gz
    >>>

    output {
        File reverted_vcf = "~{prefix}.reverted.vcf.gz"
        File reverted_tbi = "~{prefix}.reverted.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(10 + size(annotated_vcf, "GB")*2 + size(original_vcf, "GB")),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
