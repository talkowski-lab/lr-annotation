# Long-Read Annotation


### [RenameInfoFields](wdl/annotation_utils/RenameInfoFields.wdl)
This utility renames INFO fields in a VCF, replacing each given field string and its header description with a new one, optionally sharding by record count. It outputs the VCF with renamed INFO fields.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF.
- `Array[String] current_info_strings`: INFO field strings to replace.
- `Array[String] replace_info_strings`: Replacement INFO field strings, aligned to `current_info_strings`.
- `Array[String] replace_info_descriptions`: Replacement header descriptions, aligned to `replace_info_strings`.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.

Outputs:
- `renamed_vcf`: VCF with renamed INFO fields.
- `renamed_vcf_idx`: Index for the renamed VCF.


### [ReplaceSampleCalls](wdl/annotation_utils/ReplaceSampleCalls.wdl)
This utility replaces the genotype calls of samples in a cohort VCF with the calls from a set of per-sample VCFs. It outputs the updated cohort VCF.

Inputs:
- `Array[File] sample_vcfs`: Per-sample VCFs providing the replacement calls.
- `Array[File] sample_vcf_idxs`: Indexes for `sample_vcfs`.
- `File cohort_vcf`: Cohort VCF whose calls are replaced.
- `File cohort_vcf_idx`: Index for the cohort VCF.

Outputs:
- `replaced_vcf`: Cohort VCF with replaced sample calls.
- `replaced_vcf_idx`: Index for the updated VCF.
