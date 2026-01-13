### [AgglovarMerge](wdl/AgglovarMerge.wdl)
This workflow leverages [Agglovar](https://github.com/BeckLaboratory/agglovar) in order to merge calls across multiple input VCFs in a custom-defined manner.

Inputs:
- `vcfs`: Series of VCFs to merge.
- `vcf_idxs`: Indexes of VCFs to merge.
- `match_ref`: Determines whether to require a match on the REF tag.
- `match_alt`: Determines whether to require a match on the ALT tag.
- `ro_min`: Minimum required reciprocal overlap percentage.
- `size_ro_min`: Minimum required reciprocal length proportion.
- `offset_match`: Maximum allowable breakpoint distance.
- `offset_prop_max`: Maximum allowable breakpoint distance as a proportion of variant length.
- `match_prop_min`: Minimum required match on sequence context.


###  [AnnotateSTRs](archive/wdl/AnnotateSTRs.wdl)
This workflow is based on a [script](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/filter_vcf_to_STR_variants.py) developed by Ben Weisburd to annotate STRs based on their sequence context. It involves reading the nucleotide content of each variant, and comparing this to its surrounding reference genome context.

References:
- `ref_fa`.


### [BenchmarkSTRs](wdl/BenchmarkSTRs.wdl)
TODO


### [ExtractRegionFromBAM](wdl/ExtractRegionFromBAM.wdl)
TODO


### [FlagSingletonReads](wdl/FlagSingletonReads.wdl)
TODO


### [PAV](wdl/PAV.wdl)
This workflow leverages [PAV](https://github.com/BeckLaboratory/pav) in order to call variants using assemblies. 

Inputs:
- `mat_haplotypes`: Series of unaligned maternal haplotypes.
- `pat_haplotypes`: Series of unaligned paternal haplotypes.


### [TruvariMerge](wdl/TruvariMerge.wdl)
This workflow leverages [Truvari](https://github.com/ACEnglish/truvari) in order to merge VCFs whilst collapsing redundant calls. It runs `collapse` with its default parameters of 70% sequence similarity, 70% size similarity, reciprocal overlap of 0% and breakpoint distance of 500 bp.

Inputs:
- `vcfs`: Series of VCFs to merge.
- `vcf_idxs`: Indexes of VCFs to merge.
- `ref_fa`.
- `ref_fai`.