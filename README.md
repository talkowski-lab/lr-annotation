# Long-Read Annotation
This repository serves as a home for all scripts, workflows and processes for annotating long-read callsets.



## Cohort
- HGSVC3.
	- Data:
		- 65 total samples.
		- 32x aligned reads, produced after downsampling by Fabio.
		- High coverage assemblies, derived directly from HGSVC.
		- [VCF](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/).
		- [Samples](https://www.internationalgenome.org/data-portal/data-collection/hgsvc3).
	- Notes:
		- Metadata file has 67 samples, though it misses NA24385 (HG002 in Terra) from its VCF, renames NA21487 from its VCF to GM21487 and additionally includes GM19320, GM20355 & GM19129.
- HPRC2.
	- Data:
		- 232 total samples.
		- High coverage aligned reads, produced directly by Fabio.
		- High coverage assemblies, derived directly from HPRC.
		- [VCF](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/release2/minigraph-cactus/).
		- [Samples](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv).
	- Notes: 
		- Metadata file has 234 samples, as it also includes GRC38 and CHM13.
		- VCF has 232 samples, though it misses HG00272 and instead includes CHM13.
		- The HPRC Y2 workspace includes an additional sample HG03492.
- Additional Notes:
	- Overlapping Samples: 5 (HG002, HG00733, HG02818, NA19036, NA19240).
	- Mismatched Samples: NA24385 in HGSVC3 is named HG002 in HPRC.
	- The pedigrees in the Terra workspace were made to match the callset VCFs.
	- The sample IDs in the Terra workspace were made to match that of the HPRC Y2 workspace.



## Pipeline-Wide References
- `ref_fa`: Reference sequence, derived from what was used for [PAV](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/).
- `ref_fai`: Index for the reference sequence.



## Annotation Workflows
### [AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)
This workflow leverages [AnnotateVcf](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling/workflows/broad-firecloud-dsde-methods/20-AnnotateVcf) from the GATK-SV pipeline in order to annotate internal allele frequencies based on sample sexes and ancestries. It runs on all variants in the input VCF, including SVs.

Inputs:
- `sample_pop_assignments`: Two column file containing sample IDs in the first column and ancestry labels in the second column.
- `ped_file`: Six column file containing the cohort pedigree, with specifications described in [this article](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).

References:
- `par_bed`: [Panel for hg38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed) from the GATK-SV featured workspace.


### [AnnotateL1MEAIDFilter](wdl/AnnotateL1MEAIDFilter.wdl)
This workflow leverages [L1ME-AID](https://github.com/Markloftus/L1ME-AID) and [INTACT_MEI](https://github.com/xzhuo/INTACT_MEI) in order to annotate and then filter MEIs in the input VCF. It outputs a filtered version of the output from _RepeatMasker_.

Inputs:
- `rm_fa`: Output by _RepeatMasker_.
- `rm_out`: Output by _RepeatMasker_.


### [AnnotatePALMER](wdl/AnnotatePALMER.wdl)
This workflow leverages [PALMER](https://github.com/WeichenZhou/PALMER) in order to annotate MEI calls for a cohort in a given cohort VCF. It retains the genotypes present in the VCF, simply adding an INFO field `ME_TYPE` to insertions whose characteristics match those of the PALMER calls.

Inputs:
- `rm_fa`: Output by _RepeatMasker_.
- `rm_out`: Output by _RepeatMasker_.

References:
- `ref_fai`.


### [AnnotateSVAN](wdl/AnnotateSVAN.wdl)
This workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) in order to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It involves running  Tandem Repeat Finder (TRF) on the inserted or deleted sequence for each SV in the input VCF.

References:
- `ref_fa`.
- `vntr_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/VNTR_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `exons_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/EXONS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `repeats_bed`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/REPEATS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `mei_fasta`: [hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/CONSENSUS.fa) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.


### [AnnotateSVAnnotate](wdl/AnnotateSVAnnotate.wdl)
This workflow leverages [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) in order to annotate predicted functional effects for SVs. It conditionally only runs SV through this workflow, ignoring all SNVs and InDels.

References:
- `noncoding_bed`: [Panel for hg38](gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed) from the GATK-SV featured workspace.
- `coding_gtf`: [GENCODE v39](gs://talkowski-sv-gnomad-output/zero/RerunAnno/genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf) from the gnomAD workspace.


### [AnnotateVEPHail](wdl/AnnotateVEPHail.wdl)
This workflow leverages [the Ensembl Variant Effect Predictor (VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) in order to annotate predicted functional effects based on site-level information. It requires numerous references that provide context to these annotations, and uses Hail in order to run this annotation process in a more efficient and scalable manner.

References:
- `ref_fa`.
- `top_level_fa`: [hg38 Release 76](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/GRCh38.dna.toplevel.chr.fa.gz) from [Ensemble](https://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz).
- `ref_vep_cache`: [v105](gs://gcp-public-data--gnomad/resources/vep/v105/homo_sapiens_merged_vep_105_GRCh38.tar.gz) from [VEP archives](https://ftp.ensembl.org/pub/release-105/variation/). This also contains a series of additional references, which include the MANE protein coding GTF, GENCODE gene list, ClinVar annotation set etc, with the most up-to-date list of these found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).



## Downstream Workflows
### [BenchmarkAnnotations](wdl/BenchmarkAnnotations.wdl)
This workflow ingests two VCFs and finds matching variants across them in order to compare the AF & VEP annotations of these matched pairs. This serves as a degree of benchmarking, as it ensures that annotations applied to a larger cohort (e.g. gnomAD) are in line with those we annotate. It also enables the identification of variants that are outliers relative to exiting cohorts by pulling out those with a large amount of discordance in their annotation across the callsets.

The workflow undergoes multiple rounds of variant matching in order to determine matched pairs:
1. Exact match across CHROM, POS, REF and ALT.
2. Truvari match with overlap percentages of 90%, 70% and 50%.
3. Matching based on `bedtools closest`, finetuned for SVs.

Inputs:
- `vcf_truth`: VCF containing SNV & indels to evaluate against.
- `vcf_sv_truth`: VCF containing SVs to evaluate against

References:
- `ref_fa`.
- `ref_fai`.


### [FlagSingletonReads](wdl/FlagSingletonReads.wdl)
This workflow looks for variants that are supported by just 1 read, and adds the `SINGLE_READ_SUPPORT` filter status to them. Such variants are defined as those that are called in just 1 sample, and have an allele depth (AD) of just 1 in that sample.


### TODO: Merge N Annotated VCFs
This workflow will take in a VCF from each of the above annotation steps, and synthesize the annotations across each, applying logic to merge SNV vs SV information. You could even possibly specify an ordering or specific steps to run (or not run). Each of the input VCFs should be the same size though, which means that each of the above workflows should simply annotate rather than also subset.

For now, we do the following to merge our functionally annotated VCFs with our AF annotated VCFs:
```
python ./scripts/merge/merge_af_annotated_vcfs.py \
	./data/annotated_af/annotated_af.vcf.gz \
	./data/functionally_annotated/functionally_annotated.vcf.gz \
	-o ./data/annotated_merged/af_functionally_annotated.vcf.gz
```



## Additional Workflows
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


### [MinimapAlignment](wdl/MinimapAlignment.wdl)
This workflow leverages [Minimap2](https://github.com/BeckLaboratory/agglovar) in order to align assemblies to a reference.

Inputs:
- `assembly_mat`: Maternal assembly.
- `assembly_pat`: Paternal assembly.
- `minimap_flags`: Parameters to use when running Minimap2.

References:
- `ref_fa`.
- `ref_fai`.


### [PALMER](wdl/PALMER.wdl)
This workflow runs PALMER on a pair of aligned assembly haplotypes in order to generate MEI calls. It then convets the raw PALMER calls generated into a VCF, merges calls across the haplotypes to create a diploid VCF per haplotype then finally integrates these into a final VCF containing multiple MEI types.

Inputs:
- `bam_pat`: Aligned assembly for paternal haplotype.
- `bai_pat`: Index for aligned assembly for paternal haplotype.
- `bam_mat`: Aligned assembly for maternal haplotype.
- `bai_mat`: Index for aligned assembly for maternal haplotype.
- `mei_types`: Series of MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `contigs`: Series of contigs to make PALMER calls on.
- `override_palmer_calls_pat`: Optional PALMER calls for paternal haplotype to use if overriding running of PALMER.
- `override_palmer_tsd_reads_pat`: Optional PALMER TSD reads for paternal haplotype to use if overriding running of PALMER.
- `override_palmer_calls_mat`: Optional PALMER calls for maternal haplotype to use if overriding running of PALMER.
- `override_palmer_tsd_reads_mat`: Optional PALMER TSD reads for maternal haplotype to use if overriding running of PALMER.

References:
- `ref_fa`.
- `ref_fai`.


### [PAV](wdl/PAV.wdl)
This workflow leverages [PAV](https://github.com/BeckLaboratory/pav) in order to call variants using assemblies. 

Inputs:
- `mat_haplotypes`: Series of unaligned maternal haplotypes.
- `pat_haplotypes`: Series of unaligned paternal haplotypes.


### [RepeatMasker](wdl/RepeatMasker.wdl)
This workflow leverages [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) in order to annotate repeated regions in an input VCF.


### [TRGT](wdl/TRGT.wdl)
This workflow leverages [TRGT](https://github.com/PacificBiosciences/trgt) in order to genotype short-tandem repeats. 

Inputs:
- `bam`: Aligned reads.
- `bai`: Index for aligned reads.
- `sex`: Sex of sample (one of `M` or `F`).

References:
- `ref_fa`.
- `ref_fai`.
- `repeat_catalog`: [Panel for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz) from [Ben's repository](https://github.com/broadinstitute/tandem-repeat-catalog/releases) as used in All of Us Phase 2.


### [TruvariMerge](wdl/TruvariMerge.wdl)
This workflow leverages [Truvari](https://github.com/ACEnglish/truvari) in order to merge VCFs whilst collapsing redundant calls. It runs `collapse` with its default parameters of 70% sequence similarity, 70% size similarity, reciprocal overlap of 0% and breakpoint distance of 500 bp.

Inputs:
- `vcfs`: Series of VCFs to merge.
- `vcf_idxs`: Indexes of VCFs to merge.

References:
- `ref_fa`.
- `ref_fai`.


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
	1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their indexes if applicable).
	2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PEDs, metadata files etc.
	3. References - e.g. reference fasta, their indexes, catalogs used for annotations, etc.
	4. Runtime information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
	5. All RuntimeAttr? inputs - there should be one per task called, with its name reflective of the task's function.
- Workflows should take in an input `prefix` that is passed to every task that creates output files, which should be used in conjunction with a descriptive suffix when creating outputs.
- Workflow imports should not be renamed using the `as` operator.
- Workflows should never contain any blank comments - e.g. `#########################`.
- Workflows should never contain be any consecutive blank lines - i.e. they should have a maximum of one blank line at a time.
- Inputs passed to a task should not have blank lines between inputs.
- The order of inputs passed to a task should reflect their order in the inputs on the workflow level.
- Inputs passed to a task should have a space on either side of the `=` character.
- The inputs section of a task should not have blank lines between inputs.
- Tasks should always have input fields `docker` and `runtime_attr_override` defined, though what is passed to each one of these when calling the task should be explicitly named - e.g. `gatk_docker` and `runtime_attr_override_svannotate` respectively.
- Every command block within a task should begin with `set -euo pipefail` followed by a blank line.
- The default `disk_size` for a task should be calculated dynamically based on the largest sized input file - or multiple if there are several large inputs, like multiple reference fastas or input catalogs. It should be defined in-line in the default runtime attributes section, unless it is a complicated function in which it can have a dedicated variable `disk_size`.
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
- Tsks that can be generalized and used across workflows should live in `Helpers.wdl` and be imported by consumer workflows, rather than explicitly defined in a standalone workflow itself.
- Workflow file names must always match the workflow defined within them.
- All code written should use 4-space tabs for indentation.


### Python
- All code should be formatted in-line with black's formatting, which can be applied via `black .`.
- Running `flake8` should yield no errors.
- All code written should use 2-space tabs for indentation.



## Archived
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
