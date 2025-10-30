# Long-Read Annotation
Pipeline for long-read callset annotation.



## Pipeline Steps
- [AnnotateAF.wdl](#annotateafwdl).
- [AnnotateL1MEAIDFilter.wdl](#annotatel1meaidfilterwdl).
- [AnnotatePALMER.wdl](#annotatepalmerwdl).
- [AnnotateSVAN.wdl](#annotatesvanwdl).
- [AnnotateSVAnnotate.wdl](#annotatesvannotatewdl).
- [AnnotateVEPHail.wdl](#annotatevephailwdl).



## Cohort
- [HGSVC3]: 65 samples.
	- Data:
		- 32x aligned reads, produced after downsampling by Fabio.
		- High coverage assemblies, derived directly from HGSVC.
	- Links:
		- [VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/).
		- [Samples](https://www.internationalgenome.org/data-portal/data-collection/hgsvc3).
	- Metadata file has 67 samples, though it misses NA24385 (HG002 in Terra) from its VCF, renames NA21487 from its VCF to GM21487 and additionally includes GM19320, GM20355 & GM19129.
- [HPRC2]: 232 samples.
	- Data:
		- High coverage aligned reads, produced directly by Fabio.
		- High coverage assemblies, derived directly from HPRC.
	- Links:
		- [VCF](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/release2/minigraph-cactus/).
		- [Samples](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv).
	- Metadata file has 234 samples, as it also includes GRC38 and CHM13.
	- VCF has 232 samples, though it misses HG00272 and instead includes CHM13.
	- The HPRC Y2 workspace includes an additional sample HG03492.
- Overlapping: 5 samples (HG002, HG00733, HG02818, NA19036, NA19240).
- Mismatched: NA24385 in HGSVC3 is named HG002 in HPRC.
- Total: 292 samples.
- The 'modified' pedigree files were made to match the callset VCFs.
- The sample IDs in the Terra workspace were made to match that of the HPRC Y2 workspace.



## Annotation Workflows
### [AnnotateAF.wdl](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)
This workflow annotates internal allele frequencies based on sex and sample ancestrie. It runs on all variants in the input VCF, including structural variants. It is based off the `AnnotateVcf.wdl` workflow.

References:
- `GATK-SV Package`: [sv-pipeline-docker v0.28.3](https://github.com/broadinstitute/gatk-sv), which is used in the docker created by [build_docker.py](https://github.com/broadinstitute/gatk-sv/blob/5c4e659ba3747b1053b860ead5c0d7ff82768ea9/scripts/docker/build_docker.py).
- `par_bed`: [Panel for hg38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed) from the GATK-SV featured workspace.

Additional Inputs:
- `sample_pop_assignments`: Two column file containing sample IDs in the first column and ancestry labels in the second column.
- `ped_file`: Cohort PED file.


### [AnnotateL1MEAIDFilter.wdl]
This workflow runs [L1ME-AID](https://github.com/Markloftus/L1ME-AID) and then [INTACT_MEI](https://github.com/xzhuo/INTACT_MEI) to annotate and then filter MEIs in the input VCF. It outputs a filtered version of a _RepeatMasker_ file.

Additional Inputs:
- `rm_fa`: Fasta file output from _RepeatMasker_.
- `rm_out`: Output file from _RepeatMasker_.


### [AnnotatePALMER.wdl]
TODO


### [AnnotateSVAN.wdl](wdl/AnnotateSVAN.wdl)
This workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It involves running  Tandem Repeat Finder (TRF) on the inserted or deleted sequence for each SV in the input VCF.

References:
- `SVAN Package`: [Current version](https://github.com/REPBIO-LAB/SVAN) from the Github repository, which is built directly in the [docker](dockerfiles/Dockerfile.AnnotateSVAN).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.
- `vntr_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/VNTR_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `exons_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/EXONS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `repeats_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/REPEATS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `mei_fasta`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/CONSENSUS.fa) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.


### [AnnotateSVAnnotate.wdl](wdl/AnnotateSVAnnotate.wdl)
This workflow employs the GATK tool [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) to annotate predicted functional effects for structural variants. It conditionally only runs structural variants through this workflow, ignoring all SNVs and InDels.

References:
- `GATK Package`: [v4.5.0.0](https://github.com/broadinstitute/gatk), which is built in the base image `quay.io/ymostovoy/lr-utils-basic:2.0` that is used in the [docker](dockerfiles/Dockerfile.AnnotateSVAnnotate).
- `noncoding_bed`: [Panel for hg38](gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed) from the GATK-SV featured workspace.
- `coding_gtf`: [GENCODE v39](gs://talkowski-sv-gnomad-output/zero/RerunAnno/genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf) from the gnomAD workspace.


### [AnnotateVEPHail.wdl](wdl/AnnotateVEPHail.wdl)
This workflow leverages the [Ensembl Variant Effect Predictor (Ensembl VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) to annotate predicted functional effects based on site-level information. It requires numerous reference files that provide context to these annotations, and uses Hail in order to run this annotation process in a more efficient and scalable manner.

References:
- `VEP Package`: [v105](https://github.com/REPBIO-LAB/SVAN) from the VEP Dockerhub repository, which is used as the base image in the [docker](dockerfiles/Dockerfile.AnnotateVEPHail).
- `Hail Package`: [Current version](https://github.com/REPBIO-LAB/SVAN) from the Github repository, which is built directly using `pip` in the [docker](dockerfiles/Dockerfile.AnnotateVEPHail).
- `reference_fasta`: [hg38](gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta) from the GATK-SV featured workspace.
- `top_level_fa`: [hg38 Release 76](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/GRCh38.dna.toplevel.chr.fa.gz) from [Ensemble](https://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz).
- `ref_vep_cache`: [v105](gs://gcp-public-data--gnomad/resources/vep/v105/homo_sapiens_merged_vep_105_GRCh38.tar.gz) from [VEP archives](https://ftp.ensembl.org/pub/release-105/variation/).
- Additional References: Based on the pre-defined references used in the installed version of VEP, with the most up-to-date list of these found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache). These additional references include the MANE protein coding GTF, GENCODE gene list, ClinVar annotation set etc.



## Downstream Workflows
### [BenchmarkAnnotations.wdl](wdl/BenchmarkAnnotations.wdl)
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



## Additional Workflows
### [AgglovarMerge](wdl/AgglovarMerge.wdl)
TODO


### [FlagSingletonReads](wdl/FlagSingletonReads.wdl)
TODO


### [MinimapAlignment.wdl](wdl/MinimapAlignment.wdl)
TODO


### [PALMER.wdl](wdl/PALMER.wdl)
TODO


### [PALMERToVcf.wdl](wdl/PALMERToVcf.wdl)
TODO


### [PAV](wdl/PAV.wdl)
TODO


### [RepeatMasker](wdl/RepeatMasker.wdl)
TODO


### [TRGT](wdl/TRGT.wdl)
This workflow is based on a [tool](https://github.com/PacificBiosciences/trgt) developed by PacBio to annotate STRs from long-reads. It requires a reference context file that indicates tandem-repeat regions, and genotypes these to identify STRs.

References:
- `TRGT Package`: [v3.0.0](https://github.com/PacificBiosciences/trgt), as seen in the docker used in the [TRGT workflow](wdl/TRGT.wdl), as it was used for All of Us Phase 2. An earlier version of TRGT yielded per-sample costs of nearly 4x this version.
- `reference_fasta`: [hg38 with no alts](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/) from the HGSVC3 consortia.
- `repeat_catalog`: [Panel for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz) from [Ben's repository](https://github.com/broadinstitute/tandem-repeat-catalog/releases) as used in All of Us.


### [TruvariMerge](wdl/TruvariMerge.wdl)
TODO



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



## Code Conventions
### WDL
- Workflows should be structured in the following order, with each of the below separated by a blank line:
	1. Imports.
	2. Inputs.
	3. Definition of variables dynamically generated in the workflow itself - unless these require outputs from other tasks, in which case they can go after the tasks they are dependent on.
	4. Calls to tasks.
	5. Outputs.
- Tasks should be structured in the following order, with each of the below separated by a blank line:
	1. Inputs.
	2. Definition of variables dynamically generated in the task itself.
	3. Command.
	4. Outputs.
	5. Runtime settings - which should first be its default runtime settings, followed by a select first with the runtime override, then the actual runtime block.
- Workflow inputs should be structured in the following order, with each of the below separated by a blank lines:
	1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their index file if applicable).
	2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PED files, sample metadata files etc.
	3. Reference files - e.g. reference fasta, their index files, catalogs used for annotations, etc.
	4. Runtime information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
	5. All RuntimeAttr? inputs - there should be one per task called, with its name reflective of the task's function.
- Workflows should take in an input `prefix` that is passed to every task that creates output files, which should be used in conjunction with a descriptive suffix when creating output files.
- Workflow imports should not be renamed using the `as` operator.
- Workflows should never contain any blank comments - e.g. `#########################`.
- Workflows should never contain be any consecutive blank lines - i.e. they should have a maximum of one blank line at a time.
- Inputs passed to a task should not have blank lines between inputs.
- The order of inputs passed to a task should reflect their order in the inputs on the workflow level.
- Inputs passed to a task should have a space on either side of the `=` character.
- The inputs section of a task should not have blank lines between inputs.
- Tasks should always have input fields `docker` and `runtime_attr_override` defined, though what is passed to each one of these when calling the task should be explicitly named - e.g. `gatk_docker` and `runtime_attr_override_svannotate` respectively.
- Every command block within a task should begin with `set -euo pipefail` followed by a blank line.
- The default `disk_size` for a task should be calculated dynamically based on the largest sized input file - or multiple if there are several large inputs, like multiple reference fasta files or input catalogs. It should be defined in-line in the default runtime attributes section, unless it is a complicated function in which it can have a dedicated variable `disk_size`.
- The default `mem_gb`, `boot_disk_gb` and `compute_cores` for a task should be explicitly defined rather than based on an input file - it should be set based on the intensity of compute needed by that task.
- The default `preemptible_tries` for a task should always be 1.
- The default `max_retries` for a task should always be 0.
- The names of workflows and tasks should never include a `_` character within them - rather, they should always be in Pascal case.
- The names of inputs, variables and outputs should include a `_` to separate words, and be entirely lowercase unless they refer to a noun that is capitalized (e.g. PALMER or L1MEAID) - i.e. they should always be in snake case.
- There should never be any additional indentation in order to better align parts of the code to the length or horizontal/vertical spacing of other components in its section - indentation should only be applied at the start of a line.
- All mentions of `fasta` should instead use `fa` - e.g. `ref_fa` instead of `ref_fasta`.
- All mentions of `fasta_index`, `fasta_fai` or `fa_fai`  should instead use `fai` - e.g. `ref_fai` instead of `ref_fasta_index`, `ref_fasta_fai` or `ref_fa_fai`.
- All mentions of `vcf_index` or `vcf_tbi` should instead use `vcf_idx`.
- All VCFs should have suffix `_vcf`, and be coupled with a VCF index file that has a suffix `_vcf_idx`.
- General tasks that can be generalized and used across workflows should live in `Helpers.wdl` and be imported by consumer workflows rather than explicitly defined in the workflow file itself.
- Workflow file names must always match the workflow defined within them.
- All code written should use 4-space tabs for indentation.


### Python
- All code should be formatted in-line with black's formatting, which can be applied via `black .`.
- Running `flake8` should yield no errors.
- All code written should use 2-space tabs for indentation.


### Dockerfiles
TODO
