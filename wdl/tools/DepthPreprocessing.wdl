version 1.0

# Defragment GATK-gCNV CNVs per sample and merge

workflow DepthPreprocessing {
  input {
    Array[String] samples
    Array[File] genotyped_segments_vcfs
    File contig_ploidy_calls_tar
    String batch
    String sv_base_mini_docker
    String sv_pipeline_docker
    Int gcnv_qs_cutoff
    Float? defragment_max_dist
  }

  scatter (i in range(length(samples))) {
    call GcnvVcfToBed {
      input:
        sample_id = samples[i],
        sample_index = i,
        vcf = genotyped_segments_vcfs[i],
        contig_ploidy_calls_tar = contig_ploidy_calls_tar,
        sv_pipeline_docker = sv_pipeline_docker,
        qs_cutoff = gcnv_qs_cutoff
    }
  }

  scatter (i in range(length(samples))) {
    call MergeSample as MergeSample_del {
      input:
        sample_id = samples[i],
        gcnv = GcnvVcfToBed.del_bed[i],
        max_dist = defragment_max_dist,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  scatter (i in range(length(samples))) {
    call MergeSample as MergeSample_dup {
      input:
        sample_id = samples[i],
        gcnv = GcnvVcfToBed.dup_bed[i],
        max_dist = defragment_max_dist,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  call MergeSet as MergeSet_del {
    input:
      beds = MergeSample_del.sample_bed,
      svtype = "DEL",
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker
  }

  call MergeSet as MergeSet_dup {
    input:
      beds = MergeSample_dup.sample_bed,
      svtype = "DUP",
      batch = batch,
      sv_base_mini_docker = sv_base_mini_docker
  }

  output {
    File del = MergeSet_del.out
    File del_index = MergeSet_del.out_idx
    File dup = MergeSet_dup.out
    File dup_index = MergeSet_dup.out_idx
  }
}

task GcnvVcfToBed {
  input {
    File vcf
    File contig_ploidy_calls_tar
    Int sample_index
    String sample_id
    Int qs_cutoff

    String sv_pipeline_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  Int default_disk_gb = ceil(size([vcf, contig_ploidy_calls_tar], "GB") * 2) + 50

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 4]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 10])
    docker: sv_pipeline_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 3])
    noAddress: true
  }

  command <<<
    set -euo pipefail

    tar xzf ~{contig_ploidy_calls_tar}
    # The tar file contains one directory per sample given to GATK DetermineGermlineContigPloidy,
    # with the naming scheme SAMPLE_0 to SAMPLE_N-1, presumably in the order the samples were given
    # to the tool.
    calls_dir='~{"SAMPLE_" + sample_index}'
    expected_sample_id='~{sample_id}'
    actual_sample_id="$(cat "${calls_dir}/sample_name.txt")"
    if [[ "${expected_sample_id}" != "${actual_sample_id}" ]]; then
      printf 'Expected sample ID does not match actual sample ID\n' >&2
      printf 'Expected: %s\n' "${expected_sample_id}" >&2
      printf 'Actual: %s\n' "${actual_sample_id}" >&2
      printf 'Likely that sample order for this task differs from order given to GATK DetermineGermlineContigPloidy\n' >&2
      exit 1
    fi
    cp "${calls_dir}/contig_ploidy.tsv" contig_ploidy.tsv

    tabix ~{vcf}
    python /opt/WGD/bin/convert_gcnv.py \
      --cutoff ~{qs_cutoff} \
      contig_ploidy.tsv \
      ~{vcf} \
      ~{sample_id}
  >>>

  output {
    File del_bed = "~{sample_id}.del.bed"
    File dup_bed = "~{sample_id}.dup.bed"
  }
}

task MergeSample {
  input {
    File gcnv
    String sample_id
    Float? max_dist

    String sv_pipeline_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  Int default_disk_gb = ceil(size(gcnv, "GB") * 2) + 50

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 4]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 10])
    docker: sv_pipeline_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 3])
    noAddress: true
  }

  command <<<
    set -euo pipefail

    sort ~{gcnv} -k1,1V -k2,2n > ~{sample_id}.bed
    bedtools merge -i ~{sample_id}.bed -d 0 -c 4,5,6,7 -o distinct > ~{sample_id}.merged.bed
    /opt/sv-pipeline/00_preprocessing/scripts/defragment_cnvs.py \
      --max-dist ~{if defined(max_dist) then max_dist else "0.25"} ~{sample_id}.merged.bed ~{sample_id}.merged.defrag.bed
    sort -k1,1V -k2,2n ~{sample_id}.merged.defrag.bed > ~{sample_id}.merged.defrag.sorted.bed
  >>>

  output {
    File sample_bed = "~{sample_id}.merged.defrag.sorted.bed"
  }
}

task MergeSet {
  input {
    Array[File] beds
    String svtype
    String batch

    String sv_base_mini_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  Int default_disk_gb = ceil(size(beds, "GB") * 2) + 50

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 4]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 10])
    docker: sv_base_mini_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 3])
    noAddress: true
  }

  command <<<
    set -euo pipefail

    cat ~{write_lines(beds)} \
      | xargs zcat \
      | sort -k1,1V -k2,2n \
      | awk -v OFS="\t" -v svtype=~{svtype} -v batch=~{batch} '{$4=batch"_"svtype"_"NR; print}' \
      | cat <(echo -e "#chr\\tstart\\tend\\tname\\tsample\\tsvtype\\tsources") - \
      | bgzip -c > ~{batch}.~{svtype}.bed.gz;
    tabix -p bed ~{batch}.~{svtype}.bed.gz
  >>>

  output {
    File out = "~{batch}.~{svtype}.bed.gz"
    File out_idx = "~{batch}.~{svtype}.bed.gz.tbi"
  }
}

