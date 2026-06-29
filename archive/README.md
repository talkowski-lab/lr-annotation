# Long-Read Annotation
## Annotation Utilities
### [MergeSites](wdl/annotation_utils/MergeSites.wdl)
This utility merges redundant records at the site level within a VCF by collapsing near-identical deletions and insertions. Deletions are collapsed using size-, reciprocal-overlap, sequence- and sample-similarity thresholds plus a breakpoint distance, insertions using size-, sequence- and sample-similarity plus a breakpoint distance, while all other variants pass through untouched. It outputs the merged VCF.

Inputs:
- `File vcf`: VCF to merge.
- `File vcf_idx`: Index for VCF.
- `Int del_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing deletions (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for collapsing deletions (default `0.0`).
- `Float del_sample_similarity`: Minimum sample similarity for collapsing deletions (default `0.5`).
- `Float del_sequence_similarity`: Minimum sequence similarity for collapsing deletions (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for collapsing deletions (default `0.7`).
- `Int del_size_max`: Maximum deletion size to collapse, or `-1` for no maximum (default `-1`).
- `Int del_size_min`: Minimum deletion size to collapse (default `0`).
- `Int ins_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing insertions (default `200`).
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for collapsing insertions (default `0.0`).
- `Float ins_sample_similarity`: Minimum sample similarity for collapsing insertions (default `0.5`).
- `Float ins_sequence_similarity`: Minimum sequence similarity for collapsing insertions (default `0.7`).
- `Float ins_size_similarity`: Minimum size similarity for collapsing insertions (default `0.7`).
- `Int ins_size_max`: Maximum insertion size to collapse, or `-1` for no maximum (default `-1`).
- `Int ins_size_min`: Minimum insertion size to collapse (default `0`).

Outputs:
- `merged_vcf`: Site-merged VCF.
- `merged_vcf_idx`: Index for the merged VCF.


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
