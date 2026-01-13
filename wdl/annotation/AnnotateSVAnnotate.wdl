version 1.0

import "../utils/Structs.wdl"

workflow AnnotateSVAnnotate {
    input {
        File vcf
        File vcf_idx

        String prefix
        Int min_svlen
        Array[String] contigs

        File coding_gtf
        File noncoding_bed

        String annotate_sv_annotate_docker
        String gatk_docker

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
                docker = annotate_sv_annotate_docker,
                runtime_attr_override = runtime_attr_subset_vcf
          }

        call PreprocessVcf {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.doubled",
                docker = annotate_sv_annotate_docker,
                runtime_attr_override = runtime_attr_preprocess
        }

        call AnnotateFunctionalConsequences {
            input:
                vcf = PreprocessVcf.processed_vcf,
                vcf_idx = PreprocessVcf.processed_vcf_idx,
                noncoding_bed = noncoding_bed,
                coding_gtf = coding_gtf,
                prefix = "~{prefix}.~{contig}.doubled.anno_func",
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_annotate_func
        }
    }

    call ConcatVcfs as ConcatUnannotated {
        input:
            vcfs = SubsetVcf.unannotated_vcf,
            vcfs_idx = SubsetVcf.unannotated_vcf_idx,
            allow_overlaps = true,
            prefix = "~{prefix}.unannotated.concat",
            docker = annotate_sv_annotate_docker,
            runtime_attr_override = runtime_attr_concat_unannotated
    }

    call ConcatVcfs as ConcatAnnotated {
        input:
            vcfs = AnnotateFunctionalConsequences.anno_vcf,
            vcfs_idx = AnnotateFunctionalConsequences.anno_vcf_idx,
            allow_overlaps = true,
            prefix = "~{prefix}.concat",
            docker = annotate_sv_annotate_docker,
            runtime_attr_override = runtime_attr_concat_annotated
    }

    call MergeVcf {
        input:
            annotated_vcf = ConcatAnnotated.concat_vcf,
            annotated_vcf_idx = ConcatAnnotated.concat_vcf_index,
            unannotated_vcf = ConcatUnannotated.concat_vcf,
            unannotated_vcf_idx = ConcatUnannotated.concat_vcf_index,
            prefix = prefix,
            docker = annotate_sv_annotate_docker,
            runtime_attr_override = runtime_attr_merge
    }

    call PostprocessVcf {
        input:
            annotated_vcf = MergeVcf.merged_vcf,
            annotated_vcf_idx = MergeVcf.merged_vcf_index,
            original_vcf = vcf,
            original_vcf_idx = vcf_idx,
            prefix = prefix,
            docker = annotate_sv_annotate_docker,
            runtime_attr_override = runtime_attr_postprocess
    }

    output {
        File sv_annotated_vcf = PostprocessVcf.reverted_vcf
        File sv_annotated_vcf_idx = PostprocessVcf.reverted_vcf_idx
    }
}

task SubsetVcf {
    input {
        File vcf
        File vcf_idx
        String locus
        String prefix = "subset"
        Int min_svlen
        String docker
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
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File unannotated_vcf = "~{prefix}.unannotated.vcf.gz"
        File unannotated_vcf_idx = "~{prefix}.unannotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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

task PreprocessVcf {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
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
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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

task AnnotateFunctionalConsequences {
    input {
        File vcf
        File vcf_idx
        File noncoding_bed
        File coding_gtf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int java_mem_mb = 1000 * 0.7 * ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]))

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
        File anno_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

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

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File]? vcfs_idx
        Boolean allow_overlaps = false
        Boolean naive = false
        Boolean sites_only = false
        Boolean sort_vcf_list = false
        String? prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String outfile_name = prefix + ".vcf.gz"
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
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
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

task MergeVcf {
    input {
        File annotated_vcf
        File annotated_vcf_idx
        File unannotated_vcf
        File unannotated_vcf_idx
        String prefix
        String docker
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
        mem_gb: 4,
        disk_gb: 2 * ceil(size(annotated_vcf, "GB") + size(unannotated_vcf, "GB")) + 10,
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

task PostprocessVcf {
    input {
        File annotated_vcf
        File annotated_vcf_idx
        File original_vcf
        File original_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        python /opt/gnomad-lr/scripts/helpers/revert_symbalts.py ~{annotated_vcf} ~{original_vcf} | bcftools view -Oz -o ~{prefix}.reverted.vcf.gz
        tabix -p vcf ~{prefix}.reverted.vcf.gz
    >>>

    output {
        File reverted_vcf = "~{prefix}.reverted.vcf.gz"
        File reverted_vcf_idx = "~{prefix}.reverted.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(annotated_vcf, "GB") + size(original_vcf, "GB")) + 10,
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
