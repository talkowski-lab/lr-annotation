# Long-Read Annotation
This repository serves as a home for all scripts, workflows and processes for annotating long-read callsets.



## Cohort
- HPRC.
	- 232 total samples.
	- [Metadata File](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv).
		- 234 total samples.
		- Additionally includes GRC38 and CHM13.
		- Fabio's table and the DeepVariant callset also include HG03492.
- HGSVC.
	- 65 total samples.
	- [Metadata File](https://www.internationalgenome.org/data-portal/data-collection/hgsvc3).
		- 67 total samples.
		- Renames NA21487 to GM21487.
		- Duplicates NA19129 (also includes GM19129) and NA20355 (also includes GM20355).
		- Additionally includes GM19320.
		- Misses NA24385 (HG002 in Terra).
- HPRC & HGSVC Overlapping Samples: HG002, HG00733, HG02818, NA19036, NA19240.
- All of Us Phase 1.
	- 1027 total samples.
	- Metadata file extracted from VCF - all samples are unrelated and of African ancestry..


## References
- `coding_gtf`: [GENCODE v39](gs://talkowski-sv-gnomad-output/zero/RerunAnno/genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf) from the gnomAD workspace.
- `exons_bed`: [Loci for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/EXONS_hg38.bed) from the [references](https://zenodo.org/records/15229020/files/hg38.tar.gz) listed in the SVAN repository.
- `fix_variant_collisions_java`: [Script for phasing](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/FixVariantCollisions-9daf9b2.java) from the DSP Long-Read SV team.
- `genetic_maps_tsv`: [Maps for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/genetic_map_b38.tsv) from the [GLIMPSE references](https://github.com/odelaneau/GLIMPSE/tree/master/maps/genetic_maps.b38).
- `mei_catalog`: [Loci for hg38](gs://todo) is a combination of ALU and LINE loci from [RepeatMasker](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema) and SVA loci from [van Bree et al](https://genome.cshlp.org/content/suppl/2022/03/24/gr.275515.121.DC1).
- `mei_fa`: [Sequences for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/CONSENSUS.fa) from the [SVAN references](https://zenodo.org/records/15229020/files/hg38.tar.gz).
- `mei_fa_amb`: Index for `mei_fa`.
- `mei_fa_ann`: Index for `mei_fa`.
- `mei_fa_bwt`: Index for `mei_fa`.
- `mei_fa_mmi`: Index for `mei_fa`.
- `mei_fa_pac`: Index for `mei_fa`.
- `mei_fa_sa`: Index for `mei_fa`.
- `noncoding_bed`: [Panel for hg38](gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/noncoding.sort.hg38.bed) from the [GATK-SV references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `par_bed`: [Panel for hg38](gs://gatk-sv-resources-public/hg38/v0/sv-resources/resources/v1/hg38.par.bed) from the [GATK-SV references](https://broadinstitute.github.io/gatk-sv/docs/resources).
- `ref_fa`: [Sequences for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/hg38.no_alt.fa) from the [PAV references](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/).
- `ref_fai`: Index for `ref_fa`.
- `ref_vep_cache`: [Cache for v105](gs://gcp-public-data--gnomad/resources/vep/v105/homo_sapiens_merged_vep_105_GRCh38.tar.gz) from the [VEP archives](https://ftp.ensembl.org/pub/release-105/variation/), with the most up-to-date list of these found [here](https://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache).
- `repeats_bed`: [Loci for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/REPEATS_hg38.bed) from the [SVAN references](https://zenodo.org/records/15229020/files/hg38.tar.gz).
- `repeat_catalog_trgt`: [Catalog for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz) from [the TR Catalog references](https://github.com/broadinstitute/tandem-repeat-catalog/releases), as used in All of Us Phase.
- `top_level_fa`: [Release 76 sequences for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/GRCh38.dna.toplevel.chr.fa.gz) from [Ensemble](https://ftp.ensembl.org/pub/release-76/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz).
- `vntr_bed`: [Loci for hg38](gs://fc-107e0442-e00c-4bb9-9810-bbe370bda6e5/files_kj/references/VNTR_hg38.bed) from the [SVAN references](https://zenodo.org/records/15229020/files/hg38.tar.gz).



## Annotation Workflows
### [AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)
This workflow leverages [AnnotateVcf](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling/workflows/broad-firecloud-dsde-methods/20-AnnotateVcf) from the GATK-SV pipeline in order to annotate internal allele frequencies based on sample sexes and ancestries. It runs on all variants in the input VCF, including SVs.

Inputs:
- `sample_pop_assignments`: Two column file containing sample IDs in the first column and ancestry labels in the second column.
- `ped_file`: Six column file containing the cohort pedigree, with specifications described in [this article](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).
- `contigs`.
- `par_bed`.


### [AnnotateL1MEAID](wdl/annotation/AnnotateL1MEAIDFilter.wdl)
This workflow first runs _RepeatMasker_ on an input VCF. It then uses its output to run [L1ME-AID](https://github.com/Markloftus/L1ME-AID) and [INTACT_MEI](https://github.com/xzhuo/INTACT_MEI) in order to annotate and filter MEI calls.

TODO


### [AnnotateMEDs](wdl/annotation/AnnotateMEDs.wdl)
TODO


### [AnnotateMEIs](wdl/annotation/AnnotateMEIs.wdl)
TODO


### [AnnotatePALMER](wdl/annotation/AnnotatePALMER.wdl)
This workflow leverages [PALMER](https://github.com/WeichenZhou/PALMER) in order to annotate MEI calls for a cohort in a given cohort VCF. It retains the genotypes present in the VCF, simply adding an INFO field `ME_TYPE` to insertions whose characteristics match those of the PALMER calls.

Inputs:
- `mei_types`: MEI types to run on - must be a subset of [`ALU`, `SVA`, `LINE` or `HERVK`].
- `rm_fa`: Output by _RepeatMasker_.
- `rm_out`: Output by _RepeatMasker_.
- `contigs`.
- `ref_fai`.


### [AnnotateSVAN](wdl/annotation/AnnotateSVAN.wdl)
This workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) in order to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It involves running  Tandem Repeat Finder (TRF) on the inserted or deleted sequence for each SV in the input VCF.

Inputs:
- `contigs`.
- `exons_bed`.
- `mei_fasta`.
- `ref_fa`.
- `repeats_bed`.
- `vntr_bed`.


### [AnnotateSVAnnotate](wdl/annotation/AnnotateSVAnnotate.wdl)
This workflow leverages [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) in order to annotate predicted functional effects for SVs. It conditionally only runs SV through this workflow, ignoring all SNVs and InDels.

Inputs:
- `coding_gtf`.
- `contigs`.
- `noncoding_bed`.


### [AnnotateTRs](wdl/annotation/AnnotateTRs.wdl)
TODO


### [AnnotateVEPHail](wdl/annotation/AnnotateVEPHail.wdl)
This workflow leverages [the Ensembl Variant Effect Predictor (VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) in order to annotate predicted functional effects based on site-level information. It requires numerous references that provide context to these annotations, and uses Hail in order to run this annotation process in a more efficient and scalable manner.

Inputs:
- `ref_fa`.
- `ref_vep_cache`.
- `top_level_fa`.


### [BenchmarkAnnotations](wdl/annotation/BenchmarkAnnotations.wdl)
This workflow ingests two VCFs and finds matching variants across them in order to compare the AF & VEP annotations of these matched pairs. This serves as a degree of benchmarking, as it ensures that annotations applied to a larger cohort (e.g. gnomAD) are in line with those we annotate. It also enables the identification of variants that are outliers relative to exiting cohorts by pulling out those with a large amount of discordance in their annotation across the callsets.

The workflow undergoes multiple rounds of variant matching in order to determine matched pairs:
1. Exact match across CHROM, POS, REF and ALT.
2. Truvari match with overlap percentages of 90%, 70% and 50%.
3. Matching based on `bedtools closest`, finetuned for SVs.

Inputs:
- `vcf_truth`: VCF containing SNV & indels to evaluate against.
- `vcf_sv_truth`: VCF containing SVs to evaluate against.
- `ref_fa`.
- `ref_fai`.



## Annotation Utilities
### [AnnotateAlleleType](wdl/annotation_utils/AnnotateAlleleType.wdl)
TODO


### [AnnotateVcf](wdl/annotation_utils/AnnotateVcf.wdl)
TODO


### [CreatePEDAncestry](wdl/annotation_utils/CreatePEDAncestry.wdl)
TODO


### [ExtractSampleVcfs](wdl/annotation_utils/ExtractSampleVcfs.wdl)
TODO


### [FillPhasedGenotypes](wdl/annotation_utils/FillPhasedGenotypes.wdl)
TODO


### [IntegrateVcfs](wdl/annotation_utils/IntegrateHGSVC.wdl)
TODO


### [SubsetTRGTToCatalog](wdl/annotation_utils/SubsetTRGTToCatalog.wdl)
TODO


### [SubsetTsvToColumns](wdl/annotation_utils/SubsetTsvToColumns.wdl)
TODO


### [SubsetVcfToContigs](wdl/annotation_utils/SubsetVcfToContigs.wdl)
TODO



## Tools
### [Automop](wdl/tools/Automop.wdl)
TODO


### [HiPhase](wdl/tools/HiPhase.wdl)
TODO


### [HiPhaseMerge](wdl/tools/HiPhaseMerge.wdl)
TODO


### [MinimapAlignment](wdl/tools/MinimapAlignment.wdl)
This workflow leverages [Minimap2](https://github.com/BeckLaboratory/agglovar) in order to align assemblies to a reference.

Inputs:
- `assembly_mat`: Maternal assembly.
- `assembly_pat`: Paternal assembly.
- `minimap_flags`: Parameters to use when running Minimap2.
- `ref_fa`.
- `ref_fai`.


### [PALMER](wdl/tools/PALMER.wdl)
This workflow runs PALMER on a pair of aligned assembly haplotypes in order to generate MEI calls. It then convets the raw PALMER calls generated into a VCF, merges calls across the haplotypes to create a diploid VCF per haplotype and then finally integrates these into a final VCF containing multiple MEI types.

Inputs:
- `bam_pat`: Aligned assembly for paternal haplotype.
- `bai_pat`: Index for aligned assembly for paternal haplotype.
- `bam_mat`: Aligned assembly for maternal haplotype.
- `bai_mat`: Index for aligned assembly for maternal haplotype.
- `mei_types`: Series of MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `truvari_collapse_params`: Truvari parameters to use when merging across haplotypes.
- `override_palmer_calls_pat`: Optional PALMER calls for paternal haplotype, causing the workflow to bypass its execution. 
- `override_palmer_tsd_reads_pat`: Optional PALMER TSD reads for paternal haplotype, causing the workflow to bypass its execution. 
- `override_palmer_calls_mat`: Optional PALMER calls for maternal haplotype, causing the workflow to bypass its execution. 
- `override_palmer_tsd_reads_mat`: Optional PALMER TSD reads for maternal haplotype, causing the workflow to bypass its execution. 
- `contigs`.
- `ref_fa`.
- `ref_fai`.


### [PALMERDiploid](wdl/tools/PALMERDiploid.wdl)
TODO


### [PALMERMerge](wdl/tools/PALMERMerge.wdl)
TODO


### [RepeatMasker](wdl/tools/RepeatMasker.wdl)
This workflow leverages [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) in order to annotate repeated regions in an input VCF.

TODO


### [SHAPEITPhase](wdl/tools/SHAPEITPhase.wdl)
TODO


### [TRGT](wdl/tools/TRGT.wdl)
This workflow leverages [TRGT](https://github.com/PacificBiosciences/trgt) in order to genotype short-tandem repeats. 

Inputs:
- `bam`: Aligned reads.
- `bai`: Index for aligned reads.
- `sex`: Sex of sample (one of `M` or `F`).
- `ref_fa`.
- `ref_fai`.
- `repeat_catalog_trgt`.


### [TRGTLPS](wdl/tools/TRGTLPS.wdl)
TODO


### [TRGTMerge](wdl/tools/TRGTMerge.wdl)
TODO


### [VcfDist](wdl/tools/VcfDist.wdl)
TODO



## Output Schema
- `INFO/allele_length`: Allele length - positive for insertions, negative for deletions and 0 for SNVs.
- `INFO/allele_type`: Allele type, which is one of the below.
	- `snv`: Single nucleotide variant.
	- `ins`: Insertion.
	- `del`: Deletion.
	- `dup`: Tandem duplication.
	- `dup_interspersed`: Interspersed duplication.
	- `complex_dup`: Complex duplication.
	- `inv_dup`: Inverted duplication.
	- `numt`: Nuclear-mitochondrial segment.
	- `trv`: Tandem repeat.
	- `alu_ins`: ALU insertion.
	- `line_ins`: LINE insertion.
	- `sva_ins`: SVA insertion.
	- `alu_del`: ALU deletion.
	- `line_del`: LINE deletion.
	- `sva_del`: SVA deletion.
- `SOURCE`: Source of call, which is one of the below.
	- `DeepVariant`: SNV or indel call made by the DeepVariant pipeline.
	- `HPRC_SV_Integration`: Structural variant call made by the HPRC SV Integration pipeline.
	- `TRGT`: Tandem repeat call made by TRGT.
- `ORIGIN`: Origin of duplicated sequence for duplications and NUMTs.
- `SUB_FAMILY`: Sub-family for MEI calls.
- `TRID`: TR identifier for TR calls, as well as non-TR calls that that are completely enveloped by a TR call. 
- Functional Annotations.
	- `vep`: Annotations from the Variant Effect Predictor (VEP).
	- `PREDICTED_`: Annotations from SVAnnotate, which are all prefixed by `PREDICTED_`.
- gnomAD_V4 Benchmarking.
	- `gnomAD_V4_match_type`: Method for generating match, which is one of the below.
		- `EXACT_MATCH`: Exact match across CHROM, POS, REF and ALT.
		- `TRUVARI_{X}`: Truvari match requiring X% sequence similarity. 
		- `BEDTOOLS_CLOSEST`: Bedtools closest match finetuned for SVs.
	- `gnomAD_V4_match_ID`: Variant ID of matched variant.
	- `gnomAD_V4_match_source`: Source of matched variant, which is one of the below.
		- `SNV_indel`: SNV & indel callset.
		- `SV`: SV callset.
- Allele Frequencies.
	- `AN`: Count of alleles genotyped.
	- `AC`: Count of non-reference alleles.
	- `AF`: Proportion of alleles that are non-reference.
	- `NCR`: Proportion of alleles that don't have a genotype call.
	- `AP_allele`: Allele purity per-allele (multiallelic sites only).
	- `MC_allele`: Motif count per-allele (multiallelic sites only).
	- `LPS_allele`: Longest polymer sequence per-allele (multiallelic sites only).
- Filters.
	- `LARGE_SNV_INDEL`: Variant with `SOURCE = "DeepVariant"` that has  `INFO/allele_length ≥ 50`.
	- `SMALL_SV`: Variant with `SOURCE = "HPRC_SV_Integration"` that has `INFO/allele_length < 50`.
	- `TRGT_OVERLAPPED`: Variant with `SOURCE != "TRGT"` that is completely enveloped by a call with `SOURCE = "TRGT"`. 



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
- Inputs should be structured in the following order, with each of the below separated by a blank lines:
	1. Core input files that will be run through the workflow - e.g. VCFs being annotated, BAMs being analyzed etc (as well as their indexes if applicable). Also the contigs to be run on as well as the prefix.
	2. Parameters that govern how the file will be processed - e.g. prefixes, modes, input arguments to tools being called, PEDs, metadata files etc.
	3. Reference files - e.g. reference fasta, their indexes, catalogs used for annotations, etc.
	4. Runtime-related information that are not of type RuntimeAttr - e.g. docker paths, cores if applicable, sharding information if applicable.
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
- Tasks should also have a prefix input defined, which is passed and set at the workflow level - the outputs from the task should simply use the prefix along with the file type.
- Every command block within a task should begin with `set -euo pipefail` followed by a blank line.
- The default `disk_size` for a task should be calculated dynamically based on the largest sized input file - or multiple if there are several large inputs, like multiple reference fastas or input catalogs. It should be defined in-line in the default runtime attributes section, unless it is a complicated function in which it can have a dedicated variable `disk_size`.
- The default `mem_gb`, `boot_disk_gb` and `compute_cores` for a task should be explicitly defined rather than based on an input file - it should be set based on the intensity of compute needed by that task.
- The default `preemptible_tries` for a task should always be 2.
- The default `max_retries` for a task should always be 0.
- The names of workflows and tasks should never include a `_` character within them - rather, they should always be in Pascal case.
- The names of inputs, variables and outputs should include a `_` to separate words, and be entirely lowercase unless they refer to a noun that is capitalized (e.g. PALMER or L1MEAID) - i.e. they should always be in snake case.
- There should never be any additional indentation in order to better align parts of the code to the length or horizontal/vertical spacing of other components in its section - indentation should only be applied at the start of a line.
- All mentions of `fasta` should instead use `fa` - e.g. `ref_fa` instead of `ref_fasta`.
- All mentions of `fasta_index`, `fasta_fai` or `fa_fai`  should instead use `fai` - e.g. `ref_fai` instead of `ref_fasta_index`, `ref_fasta_fai` or `ref_fa_fai`.
- All mentions of `vcf_index` or `vcf_tbi` should instead use `vcf_idx`.
- All VCFs should have suffix `_vcf`, and be coupled with a VCF index file that has a suffix `_vcf_idx`.
- Tasks that can be generalized and used across workflows should live in `Helpers.wdl` and be imported by consumer workflows, rather than explicitly defined in a standalone workflow itself.
- Workflow file names must always match the workflow defined within them.


### Python
- All code should be formatted in-line with black's formatting, which can be applied via `black .`.
- All code should be compliant with `flake8`


### Codebase
- Workflows in `wdl/annotation/` should begin with _Annotate_.
- Workflows directly run in the pipeline should be in either `wdl/annotation/`, `wdl/annotation_utils/` or `wdl/tools/`.
- Workflows in `wdl/annotation/`, `wdl/annotation_utils/` and `wdl/tools/` should have an entry in `dockstore.yml`.
- Workflows in `wdl/annotation/`, `wdl/annotation_utils/` and `wdl/tools/` should be documented.


### Workspace
- All references should be passed in via workspace data.
- All dockers should be passed in via workspace data.
