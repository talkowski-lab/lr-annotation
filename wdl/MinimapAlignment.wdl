version 1.0

import "general/Structs.wdl"

workflow MinimapAlignment {
    input {
        File assembly_mat
        File assembly_pat
        String sample_id
        String minimap_flags = "-a -x asm20 --cs --eqx"
        Int minimap_threads = 32

        File ref_fasta
        File ref_fai

        String where_to_save

        String alignment_docker
        String finalize_docker

        RuntimeAttr? runtime_attr_override_align_asm2ref
        RuntimeAttr? runtime_attr_override_finalize
    }

    String workflow_name = "MinimapAlignment"
    String save_to_dir = where_to_save + "~{workflow_name}/~{sample_id}/"

    output {
        File minimap_assembled_bam_mat = save_to_dir + basename(AlignMat.bamOut)
        File minimap_assembled_bai_mat = save_to_dir + basename(AlignMat.baiOut)
        File minimap_assembled_paf_mat = save_to_dir + basename(AlignMat.pafOut)

        File minimap_assembled_bam_pat = save_to_dir + basename(AlignPat.bamOut)
        File minimap_assembled_bai_pat = save_to_dir + basename(AlignPat.baiOut)
        File minimap_assembled_paf_pat = save_to_dir + basename(AlignPat.pafOut)
    }


    call AlignAssembly as AlignMat { 
        input:
            assembly_fa = assembly_mat,
            sample_id = sample_id,
            flags = minimap_flags, 
            threads = minimap_threads,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            hap = 1,
            docker = alignment_docker,
            runtime_attr_override = runtime_attr_override_align_asm2ref
    }

    call AlignAssembly as AlignPat { 
        input:
            assembly_fa = assembly_pat,
            sample_id = sample_id,
            flags = minimap_flags,
            threads = minimap_threads,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            hap = 2,
            docker = alignment_docker,
            runtime_attr_override = runtime_attr_override_align_asm2ref
    }

    call FinalizeToDir as SaveBothHapsFiles { 
        input:
            files = [AlignMat.bamOut, AlignMat.pafOut, AlignMat.baiOut, AlignPat.bamOut, AlignPat.pafOut, AlignPat.baiOut],
            outdir = save_to_dir,
            docker = finalize_docker,
            runtime_attr_override = runtime_attr_override_finalize
    }
}

task AlignAssembly {
    input {
        File assembly_fa
        String sample_id
        File ref_fasta
        File ref_fai
        Int hap
        String flags

        Int threads
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String out_prefix = "~{sample_id}-asm_h~{hap}.minimap2"
    Int mm2_threads = threads - 4

    command <<<
        set -euo pipefail

        minimap2 \
            -t ~{mm2_threads} \
            ~{flags} \
            ~{ref_fasta} \
            ~{assembly_fa} \
        | samtools sort -@4 -o "~{out_prefix}.bam"

        samtools index -@3 "~{out_prefix}.bam"
        
        samtools view -h "~{out_prefix}.bam" \
        | k8 $(which paftools.js) sam2paf \
            -L \
            - \
        > "~{out_prefix}.paf"
    >>>

    output {
        File bamOut = "~{out_prefix}.bam"
        File baiOut = "~{out_prefix}.bam.bai"
        File pafOut = "~{out_prefix}.paf"
    }

    Int disk_size = 2*ceil(size(assembly_fa, "GB") + size(ref_fasta, "GB")) + 5

    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             4*threads,
        disk_gb:            disk_size,
        boot_disk_gb:       20,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
        runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FinalizeToDir {
    input {
        Array[File] files
        Array[String]? file_names
        String outdir
        File? keyfile
        
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")
    Boolean fail = if(defined(file_names)) then length(select_first([file_names])) != length(files) else false
    # this variable is defined because of meta-programing:
    # Cromwell generates the script to be executed at runtime (duing the run of the workflow),
    # but also at "compile time" when looked from the individual task perspective--the task is "compiled" right before it is run.
    # so optional variables, if not specified, cannot be used in the command section because at that "compile time", they are undefined
    # here we employ a hack:
    # if the optional input file_names isn't provided, it's not used anyway, so we don't worry about the literal correctness of
    # the variable's values--the variable used in generating the script--but only care that it is defined.
    Array[String] names_for_cromwell = select_first([file_names, ["correctness_doesnot_matter_here"]])

    command <<<
        set -euo pipefail

        if ~{fail}; then echo "input files and file_names don't have the same length!" && exit 1; fi

        if ~{defined(file_names)}; then
            paste \
                ~{write_lines(files)} \
                ~{write_lines(names_for_cromwell)} \
            > file_and_customname.tsv
            while IFS=$'\t' read -r ff nn; do
                gcloud storage cp \
                    "${ff}" \
                    "~{gcs_output_dir}"/"${nn}"
            done < file_and_customname.tsv
        else
            cat ~{write_lines(files)} | \
            gcloud storage cp -I "~{gcs_output_dir}"
        fi
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
