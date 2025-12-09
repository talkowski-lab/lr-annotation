###  [AnnotateSTRs](archive/wdl/AnnotateSTRs.wdl)
This workflow is based on a [script](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/filter_vcf_to_STR_variants.py) developed by Ben Weisburd to annotate STRs based on their sequence context. It involves reading the nucleotide content of each variant, and comparing this to its surrounding reference genome context.

References:
- `ref_fa`.


### [Merge AF Annotated VCFs](archive/scripts/merge/merge_af_annotated_vcfs.py)
This script takes in multiple VCFs containing various AF-related annotations, as produced by the AF annotation workflow, and synthesizes them into an integrated VCF in which each record contains the annotations for that record across all the input VCFs.

Below illustrates an example of its use:
```
python ./scripts/merge/merge_af_annotated_vcfs.py \
    ./data/annotated_af/annotated.ancestry_sex.vcf.gz \
    ./data/annotated_af/annotated.ancestry.vcf.gz \
    ./data/annotated_af/annotated.sex.vcf.gz \
    -o ./data/annotated_af/annotated_af.vcf.gz
```


### [Merge Functionally Annotated VCFs](archive/scripts/merge/merge_functionally_annotated_vcfs.py)
This script takes in VCFs produced by the functional annotation workflows, and synthesizes them into an integrated VCF in which any < 50 bp contains only the VEP annotations while any variant ≥ 50 bp contains only SVAnnotate annotations.

Below illustrates an example of its use:
```
python ./scripts/merge/merge_functionally_annotated_vcfs.py \
    --vep-vcf ./data/functionally_annotated/vep.vcf.gz \
    --svannotate-vcf ./data/functionally_annotated/svannotate.vcf.gz \
    -o ./data/functionally_annotated/functionally_annotated.vcf.gz
```
