version 1.0

workflow GenotypeDepth {
  input {
    File vcf
    File vcf_index
    File training_intervals
    File median_coverage
    File rd
    File rd_index
    File reference_dict
    File ploidy_table
    File? depth_exclusion_intervals
    File? depth_exclusion_intervals_index

    String chr_x = "chrX"
    String chr_y = "chrY"

    String output_prefix
    String gatk_docker
  }

  call TrainSVGenotyping {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      median_coverage = median_coverage,
      rd = rd,
      rd_index = rd_index,
      reference_dict = reference_dict,
      ploidy_table = ploidy_table,
      chr_x = chr_x,
      chr_y = chr_y,
      depth_exclusion_intervals = depth_exclusion_intervals,
      depth_exclusion_intervals_index = depth_exclusion_intervals_index,
      output_prefix = output_prefix,
      gatk_docker = gatk_docker
  }

  call GenotypeSVs {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      output_prefix = output_prefix,
      median_coverage = median_coverage,
      rd = rd,
      rd_index = rd_index,
      reference_dict = reference_dict,
      ploidy_table = ploidy_table,
      depth_exclusion_intervals = depth_exclusion_intervals,
      depth_exclusion_intervals_index = depth_exclusion_intervals_index,
      rd_table = TrainSVGenotyping.rd_table,
      gatk_docker = gatk_docker
  }

  output {
    File genotyped_depth_vcf = GenotypeSVs.genotyped_vcf
    File genotyped_depth_vcf_index = GenotypeSVs.genotyped_vcf_index
  }
}

task TrainSVGenotyping {
  input {
    File vcf
    File vcf_index
    File training_intervals
    File median_coverage
    File rd
    File rd_index
    File reference_dict
    File ploidy_table
    String chr_x
    String chr_y
    File? depth_exclusion_intervals
    File? depth_exclusion_intervals_index
    String output_prefix

    String gatk_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  Int default_disk_gb = ceil(size([vcf, rd], "GB") * 2) + 32

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 16]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 10])
    docker: gatk_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 1])
  }

  Int command_mem_mb = ceil(select_first([mem_gib, 16]) * 800)

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx${command_mem_mb}" TrainSVGenotyping \
      -XL '~{chr_x}' -XL '~{chr_y}' \
      -V '~{vcf}' \
      --training_intervals '~{training_intervals}' \
      -O '~{output_prefix}.vcf.gz' \
      --median_coverage '~{median_coverage}' \
      --rd-file '~{rd}' \
      --sequence-dictionary '~{reference_dict}' \
      --ploidy-table '~{ploidy_table}' \
      ~{"--depth-exclusion-intervals '" +  depth_exclusion_intervals + "'"} \
      --output-dir ./ \
      --output-prefix '~{output_prefix}'

  >>>

  output {
    File rd_table = "~{output_prefix}.rd_geno_params.tsv"
  }
}

task GenotypeSVs {
  input {
    File vcf
    File vcf_index
    String output_prefix
    File median_coverage
    File rd
    File rd_index
    File reference_dict
    File ploidy_table
    File? depth_exclusion_intervals
    File? depth_exclusion_intervals_index
    File rd_table

    String gatk_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries

  }

  Int default_disk_gb = ceil(size([vcf, rd], "GB") * 2) + 32

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 8]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 10])
    docker: gatk_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 1])
  }

  Int command_mem_mb = ceil(select_first([mem_gib, 8]) * 800)

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{command_mem_mb}" GenotypeSVs \
      -V '~{vcf}' \
      -O '~{output_prefix}.vcf.gz' \
      --median-coverage '~{median_coverage}' \
      --rd-file '~{rd}' \
      --sequence-dictionary '~{reference_dict}' \
      --ploidy-table '~{ploidy_table}' \
      ~{"--depth-exclusion-intervals '" +  depth_exclusion_intervals + "'"} \
      --rd-table '~{rd_table}'
  >>>

  output {
    File genotyped_vcf = "~{output_prefix}.vcf.gz"
    File genotyped_vcf_index = "~{output_prefix}.vcf.gz.tbi"
  }
}
