version 1.0

# Workflow: Split per-contig VCFs into SNV/indel/SV, run svtk vcf2bed, then
# concatenate across contigs for each variant class.
#
# VID classification (3rd column of VCF):
#   SNV    : last two fields separated by '-' are each a single nucleotide base
#             e.g. chr21-5010688-C-A
#   indel  : type token is DEL or INS (or other non-single-base token) AND size < 50 bp
#             e.g. chr21-5012503-DEL-1, chr21-5013525-DEL-2
#   SV     : type token is DEL, INS, DUP, INV, TRV or other non-SNV token AND size >= 50 bp
#             e.g. chr21-5264269-DUP-341, chr21-5249967-INS-13386

workflow SplitAndConvertVcfWorkflow {
  input {
    Array[File] input_vcfs         # one per contig
    Array[File] input_vcf_indices  # matched .tbi files

    String svtk_docker  = "us.gcr.io/broad-dsde-methods/talkowski-sv-pipeline:latest"
    String bcftools_docker = "staphb/bcftools:1.17"

    Int sv_size_threshold   = 50   # variants >= this size (bp) classified as SV

    Int additional_disk_gb  = 20
    Int machine_mem_gb      = 8
    Int preemptible_attempts = 1
  }

  # ── Task 1 ─────────────────────────────────────────────────────────────────
  # Split each per-contig VCF into SNV / indel / SV VCFs
  scatter (i in range(length(input_vcfs))) {
    call SplitVcfByVariantClass {
      input:
        input_vcf           = input_vcfs[i],
        input_vcf_index     = input_vcf_indices[i],
        sv_size_threshold   = sv_size_threshold,
        docker              = bcftools_docker,
        additional_disk_gb  = additional_disk_gb,
        machine_mem_gb      = machine_mem_gb,
        preemptible_attempts = preemptible_attempts

  # ── Task 2 ─────────────────────────────────────────────────────────────────
  # Run svtk vcf2bed independently for each class and each contig (scattered)
    call VcfToBedSnv {
      input:
        input_vcf            = SplitVcfByVariantClass.snv_vcf[i],
        docker               = svtk_docker,
        additional_disk_gb   = additional_disk_gb,
        machine_mem_gb       = machine_mem_gb,
        preemptible_attempts = preemptible_attempts
    }

    call VcfToBedIndel {
      input:
        input_vcf            = SplitVcfByVariantClass.indel_vcf[i],
        docker               = svtk_docker,
        additional_disk_gb   = additional_disk_gb,
        machine_mem_gb       = machine_mem_gb,
        preemptible_attempts = preemptible_attempts
    }

    call VcfToBedSv {
      input:
        input_vcf            = SplitVcfByVariantClass.sv_vcf[i],
        docker               = svtk_docker,
        additional_disk_gb   = additional_disk_gb,
        machine_mem_gb       = machine_mem_gb,
        preemptible_attempts = preemptible_attempts
    }
  }

  # ── Task 3 ─────────────────────────────────────────────────────────────────
  # Concatenate each BED class across all contigs — three independent tasks

  call ConcatBeds as ConcatSnvBeds {
    input:
      input_beds           = VcfToBedSnv.output_bed,
      output_prefix        = "merged.SNV",
      additional_disk_gb   = additional_disk_gb,
      machine_mem_gb       = machine_mem_gb,
      preemptible_attempts = preemptible_attempts
  }

  call ConcatBeds as ConcatIndelBeds {
    input:
      input_beds           = VcfToBedIndel.output_bed,
      output_prefix        = "merged.indel",
      additional_disk_gb   = additional_disk_gb,
      machine_mem_gb       = machine_mem_gb,
      preemptible_attempts = preemptible_attempts
  }

  call ConcatBeds as ConcatSvBeds {
    input:
      input_beds           = VcfToBedSv.output_bed,
      output_prefix        = "merged.SV",
      additional_disk_gb   = additional_disk_gb,
      machine_mem_gb       = machine_mem_gb,
      preemptible_attempts = preemptible_attempts
  }

  output {
    # Per-contig split VCFs
    Array[File] per_contig_snv_vcfs   = SplitVcfByVariantClass.snv_vcf
    Array[File] per_contig_indel_vcfs = SplitVcfByVariantClass.indel_vcf
    Array[File] per_contig_sv_vcfs    = SplitVcfByVariantClass.sv_vcf

    # Per-contig BED files
    Array[File] per_contig_snv_beds   = VcfToBedSnv.output_bed
    Array[File] per_contig_indel_beds = VcfToBedIndel.output_bed
    Array[File] per_contig_sv_beds    = VcfToBedSv.output_bed

    # Merged BEDs across all contigs
    File merged_snv_bed   = ConcatSnvBeds.merged_bed
    File merged_indel_bed = ConcatIndelBeds.merged_bed
    File merged_sv_bed    = ConcatSvBeds.merged_bed
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task 1: Split a single per-contig VCF into SNV / indel / SV
# ─────────────────────────────────────────────────────────────────────────────
task SplitVcfByVariantClass {
  input {
    File input_vcf
    File input_vcf_index
    Int  sv_size_threshold   = 50
    String docker
    Int additional_disk_gb   = 20
    Int machine_mem_gb       = 8
    Int preemptible_attempts = 1
  }

  String base = basename(input_vcf, ".vcf.gz")
  Int disk_gb = ceil(size(input_vcf, "GB") * 6) + additional_disk_gb

  command <<<
    set -euo pipefail

    VCF="~{input_vcf}"
    BASE="~{base}"
    SV_THRESH=~{sv_size_threshold}

    # Extract VCF header once
    bcftools view -h "$VCF" > header.vcf

    # Classify variants by VID (3rd column / ID field):
    #   SNV   : 4th field (type token) is a single nucleotide [ACGTN]
    #   indel : 4th field is NOT a single base AND size (4th field after split) < sv_size_threshold
    #   SV    : 4th field is NOT a single base AND size >= sv_size_threshold
    bcftools view -H "$VCF" | awk -v sv_thresh="$SV_THRESH" '
    BEGIN { OFS="\t" }
    {
      vid = $3
      n = split(vid, parts, "-")
      type_tok = parts[n-1]   # second-to-last field (e.g. "C", "DEL", "INS", "DUP")
      size_tok = parts[n]     # last field (e.g. "A", "1", "50", "341")

      # SNV: type_tok is a single nucleotide base
      if (type_tok ~ /^[ACGTNacgtn]$/) {
        print > "snv_records.vcf"
      } else {
        # indel or SV: size is numeric last token
        sz = size_tok + 0
        if (sz >= sv_thresh) {
          print > "sv_records.vcf"
        } else {
          print > "indel_records.vcf"
        }
      }
    }'

    # Build output VCFs from header + records (handle empty record files gracefully)
    for cls in snv indel sv; do
      out="${BASE}.${cls}.vcf.gz"
      if [ -s "${cls}_records.vcf" ]; then
        cat header.vcf "${cls}_records.vcf" | bcftools view -Oz -o "$out"
      else
        # Empty: write header-only VCF
        bcftools view -Oz -o "$out" header.vcf
      fi
      bcftools index -t "$out"
    done

    # Emit record counts for QC
    echo "SNV:   $([ -f snv_records.vcf ]   && wc -l < snv_records.vcf   || echo 0)"
    echo "Indel: $([ -f indel_records.vcf ] && wc -l < indel_records.vcf || echo 0)"
    echo "SV:    $([ -f sv_records.vcf ]    && wc -l < sv_records.vcf    || echo 0)"
  >>>

  output {
    File snv_vcf   = "~{base}.snv.vcf.gz"
    File indel_vcf = "~{base}.indel.vcf.gz"
    File sv_vcf    = "~{base}.sv.vcf.gz"
    File snv_vcf_index   = "~{base}.snv.vcf.gz.tbi"
    File indel_vcf_index = "~{base}.indel.vcf.gz.tbi"
    File sv_vcf_index    = "~{base}.sv.vcf.gz.tbi"
  }

  runtime {
    docker:      docker
    cpu:         2
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task 2 helpers: svtk vcf2bed for each variant class (separate tasks so
#                 SNV / indel / SV scatter arms run fully independently)
# ─────────────────────────────────────────────────────────────────────────────
task VcfToBedSnv {
  input {
    File input_vcf
    String docker
    Int additional_disk_gb   = 20
    Int machine_mem_gb       = 8
    Int preemptible_attempts = 1
  }

  String base = basename(input_vcf, ".vcf.gz")
  Int disk_gb = ceil(size(input_vcf, "GB") * 4) + additional_disk_gb

  command <<<
    set -euo pipefail
    svtk vcf2bed -i ALL --include-filters ~{input_vcf} ~{base}.bed
  >>>

  output {
    File output_bed = "~{base}.bed"
  }

  runtime {
    docker:      docker
    cpu:         2
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task VcfToBedIndel {
  input {
    File input_vcf
    String docker
    Int additional_disk_gb   = 20
    Int machine_mem_gb       = 8
    Int preemptible_attempts = 1
  }

  String base = basename(input_vcf, ".vcf.gz")
  Int disk_gb = ceil(size(input_vcf, "GB") * 4) + additional_disk_gb

  command <<<
    set -euo pipefail
    svtk vcf2bed -i ALL --include-filters ~{input_vcf} ~{base}.bed
  >>>

  output {
    File output_bed = "~{base}.bed"
  }

  runtime {
    docker:      docker
    cpu:         2
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

task VcfToBedSv {
  input {
    File input_vcf
    String docker
    Int additional_disk_gb   = 20
    Int machine_mem_gb       = 8
    Int preemptible_attempts = 1
  }

  String base = basename(input_vcf, ".vcf.gz")
  Int disk_gb = ceil(size(input_vcf, "GB") * 4) + additional_disk_gb

  command <<<
    set -euo pipefail
    svtk vcf2bed -i ALL --include-filters ~{input_vcf} ~{base}.bed
  >>>

  output {
    File output_bed = "~{base}.bed"
  }

  runtime {
    docker:      docker
    cpu:         2
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task 3: Concatenate BEDs across contigs (used independently per class)
# ─────────────────────────────────────────────────────────────────────────────
task ConcatBeds {
  input {
    Array[File] input_beds
    String output_prefix
    Int additional_disk_gb   = 20
    Int machine_mem_gb       = 8
    Int preemptible_attempts = 1
  }

  Int disk_gb = ceil(size(input_beds, "GB") * 4) + additional_disk_gb

  command <<<
    set -euo pipefail

    cat ~{sep=' ' input_beds} > ~{output_prefix}.bed
  >>>

  output {
    File merged_bed = "~{output_prefix}.bed"
  }

  runtime {
    cpu:         2
    memory:      machine_mem_gb + " GB"
    disks:       "local-disk " + disk_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries:  preemptible_attempts
  }
}
