# LR Annotation
This repository contains code for processing cohorts sequenced with long-read genetic sequencing methods - from variant calling and callset integration, all the way through to phasing and annotation. The pipeline was run on two cohorts as part of its initial release - a combined set of 292 samples from the HPRC & HGSVC consortia and 1027 samples from Phase 1 of the All of Us initiative.

The pipeline covers three broad stages. 
1. Callset Generation: Per-sample variant calls from DeepVariant (for SNVs/indels), a series of SV callers and TRGT (for tandem repeats) are integrated into cohort-level VCFs and filtered.
2. Callset Processing: A series of preprocessing steps convert these VCFs into the desired format. Then, variants are physically phased per sample using HiPhase, and these phase blocks are then linked via backbone phasing. 
3. Callset Annotation: A suite of annotation workflows characterizes each variant - calling mobile element insertions (MEIs) and deletions (MEDs) using PALMER, SVAN and L1ME-AID, identifying duplications and NUMTs, annotating functional effects via VEP and SVAnnotate, computing in-silico predictor scores, assigning allele frequencies, linking variants to external databases (dbSNP, dbVaR, gnomAD) and more.

**Stack:** All workflows are written in WDL and executed via Cromwell on Terra, which is built on top of GCP. The code logic is a mixture of Python, R, Bash and Hail. 


## Documentation
- [Annotations](docs/annotations.md) - VCF INFO fields, FORMAT fields and filter definitions.
- [Cohort](docs/cohort.md) - sample cohorts, sizes and metadata sources.
- [Conventions](docs/conventions.md) - WDL and Python style conventions.
- [Pipeline](docs/pipeline.md) - end-to-end pipeline description covering callset generation, preprocessing and annotation.
- [References](docs/references.md) - all reference files and their GCS locations.
- [Repository Structure](docs/repository-structure.md) - directory layout, Dockerfile build/push conventions, scripts overview and CI/CD.
- [Workflows](docs/workflows.md) - annotation workflows, annotation utilities and tool wrappers with inputs and outputs.
