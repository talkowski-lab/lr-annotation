# Long-Read Annotation
Pipeline for long-read callset annotation.



## Pipeline Steps
- [AnnotateSVAN.wdl](#annotatesvanwdl).
- [AnnotatePALMER.wdl](#annotatepalmerwdl).
- [AnnotateL1MEAIDFilter.wdl](#annotatel1meaidfilterwdl).
- [AnnotateVEPHail.wdl](#annotatevephailwdl).
- [AnnotateSVAnnotate.wdl](#annotatesvannotatewdl).
- [AnnotateAF.wdl](#annotateafwdl).



## Cohort
- [HGSVC3]: 65 samples.
- [HPRC2]: 231 samples.
- Overlapping: 5 samples (HG002, HG00733, HG02818, NA19036, NA19240).
- Missing: CHM13 (only in HPRC VCF in place of HG00272), HG03492 (only in HPRC2 Github).
- Total: 291 samples.



## Annotation Workflows
### [AnnotateSVAN.wdl](wdl/AnnotateSVAN.wdl)
This workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It involves running  Tandem Repeat Finder (TRF) on the inserted or deleted sequence for each SV in the input VCF.

References:
- `SVAN Package`: [Current version](https://github.com/REPBIO-LAB/SVAN) from the Github repository, which is built directly in the [docker](dockerfiles/Dockerfile.AnnotateSVAN).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.
- `vntr_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/VNTR_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `exons_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/EXONS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `repeats_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/REPEATS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `mei_fasta`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/CONSENSUS.fa) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.


## [AnnotatePALMER.wdl]
TODO


## [AnnotateL1MEAIDFilter.wdl]
TODO


### [AnnotateVEPHail.wdl](wdl/AnnotateVEPHail.wdl)
This workflow leverages the [Ensembl Variant Effect Predictor (Ensembl VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) to annotate predicted functional effects based on site-level information. It requires numerous reference files that provide context to these annotations, and uses Hail in order to run this annotation process in a more efficient and scalable manner.

References:
- `VEP Package`: [v105](https://github.com/REPBIO-LAB/SVAN) from the VEP Dockerhub repository, which is used as the base image in the [docker](dockerfiles/Dockerfile.AnnotateVEPHail).
- `Hail Package`: [Current version](https://github.com/REPBIO-LAB/SVAN) from the Github repository, which is built directly using `pip` in the [docker](dockerfiles/Dockerfile.AnnotateVEPHail).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.
- `top_level_fa`: [hg38 Release 76](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/GRCh38.dna.toplevel.chr.fa.gz) from [Ensemble](https://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz).
- `ref_vep_cache`: [v105](gs://gcp-public-data--gnomad/resources/vep/v105/homo_sapiens_merged_vep_105_GRCh38.tar.gz) from [VEP archives](https://ftp.ensembl.org/pub/release-105/variation/).
- Additional References: Based on the pre-defined references used in the installed version of VEP, with the most up-to-date list of these found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache). These additional references include the MANE protein coding GTF, GENCODE gene list, ClinVar annotation set etc.


### [AnnotateSVAnnotate.wdl](wdl/AnnotateSVAnnotate.wdl)
This workflow employs the GATK tool [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) to annotate predicted functional effects for structural variants. It conditionally only runs structural variants through this workflow, ignoring all SNVs and InDels.

References:
- `GATK Package`: [v4.5.0.0](https://github.com/broadinstitute/gatk), which is built in the base image `quay.io/ymostovoy/lr-utils-basic:2.0` that is used in the [docker](dockerfiles/Dockerfile.AnnotateSVAnnotate).
- `noncoding_bed`: [Panel for hg38](gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed) from the GATK-SV featured workspace.
- `coding_gtf`: [GENCODE v39](gs://talkowski-sv-gnomad-output/zero/RerunAnno/genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf) from the gnomAD workspace.


### [AnnotateAF.wdl](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)
This workflow, which is located in the `kj_project_gnomad_lr` branch of the GATK-SV repository, annotates the internal allele frequencies (AFs) based on  sample ancestries and sexes. It runs all variants in the input VCF through this workflow, including structural variants. It is based off the `AnnotateVcf.wdl` workflow.

References:
- `GATK-SV Package`: [sv-pipeline-docker v0.28.3](https://github.com/broadinstitute/gatk-sv), which is used in the docker created by [build_docker.py](https://github.com/broadinstitute/gatk-sv/blob/5c4e659ba3747b1053b860ead5c0d7ff82768ea9/scripts/docker/build_docker.py).
- `par_bed`: [Panel for hg38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed) from the GATK-SV featured workspace.



## Additional Workflows
### [BenchmarkAnnotations.wdl](wdl/BenchmarkAnnotations.wdl)
TODO


### [FlagSingletonReads](wdl/FlagSingletonReads.wdl)
TODO


### TODO: Merge N Annotated VCFs
This workflow should take in a VCF from each of the above annotation steps, and synthesize the annotations across each, applying logic to merge SNV vs SV information. You could even possibly specify an ordering or specific steps to run (or not run). Each of the input VCFs should be the same size though, which means that each of the above workflows should simply annotate rather than also subset.

For now, we do the following to merge our functionally annotated VCFs with our AF annotated VCFs:
```
python ./scripts/merge/merge_af_annotated_vcfs.py \
	./data/annotated_af/annotated_af.vcf.gz \
	./data/functionally_annotated/functionally_annotated.vcf.gz \
	-o ./data/annotated_merged/af_functionally_annotated.vcf.gz
```


### [MinimapAlignment.wdl](wdl/MinimapAlignment.wdl)
TODO


### [PALMER.wdl](wdl/PALMER.wdl)
TODO


### [PALMERToVcf.wdl](wdl/PALMERToVcf.wdl)
TODO


### [RepeatMasker](wdl/RepeatMasker.wdl)
TODO


### [TRGT](wdl/TRGT.wdl)
This workflow is based on a [tool](https://github.com/PacificBiosciences/trgt) developed by PacBio to annotate STRs from long-reads. It requires a reference context file that indicates tandem-repeat regions, and genotypes these to identify STRs.

References:
- `TRGT Package`: [v3.0.0](https://github.com/PacificBiosciences/trgt), which is the current version in the repository, and is included directly in the manually-listed docker in the [TRGT workflow](wdl/TRGT.wdl).
- `reference_fasta`: [hg38 with no alts](gs://fc-8c3900db-633f-477f-96b3-fb31ae265c44/resources/references/grch38_noalt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa) from the All of Us workspace.
- `repeat_catalog`: [Panel for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz) from [Ben's repository](https://github.com/broadinstitute/tandem-repeat-catalog/releases) as used in All of Us.



## Archived Workflows
###  [AnnotateSTRs.wdl](archive/wdl/AnnotateSTRs.wdl)
This workflow is based on a [script](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/filter_vcf_to_STR_variants.py) developed by Ben Weisburd to annotate STRs based on their sequence context. It involves reading the nucleotide content of each variant, and comparing this to its surrounding reference genome context.

References:
- `STR Analysis Package`: [Current version](https://github.com/broadinstitute/str-analysis/tree/main) from the Github repository, which is built directly in the [docker](dockerfiles/Dockerfile.AnnotateSTRs).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.


### [Merge AF Annotated VCFs](archive/scripts/merge/merge_af_annotated_vcfs.py)
This script takes in multiple VCFs containing various AF-related annotations, as produced by the AF annotation workflow, and synthesizes them into a single VCF in which each record contains the annotations for that record across all the input VCFs.

Below illustrates an example of its use:
```
python ./scripts/merge/merge_af_annotated_vcfs.py \
	./data/annotated_af/annotated.ancestry_sex.vcf.gz \
	./data/annotated_af/annotated.ancestry.vcf.gz \
	./data/annotated_af/annotated.sex.vcf.gz \
	-o ./data/annotated_af/annotated_af.vcf.gz
```


### [Merge Functionally Annotated VCFs](archive/scripts/merge/merge_functionally_annotated_vcfs.py)
This script takes in VCFs produced by the functional annotation workflows, and synthesizes them into a single VCF in which any < 50 bp contains only the VEP annotations while any variant ≥ 50 bp contains only SVAnnotate annotations.

Below illustrates an example of its use:
```
python ./scripts/merge/merge_functionally_annotated_vcfs.py \
	--vep-vcf ./data/functionally_annotated/vep.vcf.gz \
	--svannotate-vcf ./data/functionally_annotated/svannotate.vcf.gz \
	-o ./data/functionally_annotated/functionally_annotated.vcf.gz
```
