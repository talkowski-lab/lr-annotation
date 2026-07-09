# LR Annotation
This repository contains all code for processing long-read sequencing cohorts end-to-end — from variant calling and callset integration through phasing and annotation. It was run on two cohorts: a combined HPRC/HGSVC set (297 samples) and All of Us Phase 1 (1027 samples).

The pipeline covers three broad stages. First, per-sample variant calls from DeepVariant (SNVs/indels), multiple SV callers and TRGT (tandem repeats) are integrated into cohort-level VCFs and filtered. Second, variants are phased per sample using HiPhase and extended into longer phase blocks via backbone phasing. Third, a suite of annotation workflows characterizes each variant — calling mobile element insertions (MEIs) and deletions (MEDs) using PALMER, SVAN and L1ME-AID, identifying duplications and NUMTs, annotating functional effects via VEP and SVAnnotate, computing in-silico predictor scores, assigning allele frequencies, and linking variants to external databases (dbSNP, dbVar, gnomAD).

**Stack:** WDL 1.0 (executed via Cromwell on Terra), Python 3.8+, Bash, Hail. Cloud platform is GCP (GCS URIs `gs://`). Reference genome is GRCh38.


## Repository Layout
```
wdl/
  annotation/        # Main annotation workflows (prefix: Annotate*)
  annotation_utils/  # VCF manipulation and utility workflows
  tools/             # Individual bioinformatics tool wrappers
  utils/             # Shared structs (Structs.wdl) and helper tasks (Helpers.wdl)
scripts/
  vep/               # VEP and Hail-based annotation scripts
  mei/               # MEI analysis scripts
  benchmark/         # Benchmarking scripts
  miscellaneous/     # Other utility scripts
dockerfiles/         # Dockerfile.<ToolName> for each container
data/
  references/        # Reference genomes, catalogs, BED files
  base_vcfs/         # Test/base VCF files
  metadata/          # Sample metadata, pedigrees, ancestry
  qc_annotations/    # QC annotation outputs
docs/                # Extended documentation
archive/             # Deprecated workflows, scripts, and dockerfiles (do not modify)
```


## Documentation

- [Cohort](docs/cohort.md) — sample cohorts, sizes and metadata sources.
- [References](docs/references.md) — all reference files and their GCS locations.
- [Workflows](docs/workflows.md) — annotation workflows, annotation utilities and tool wrappers with inputs and outputs.
- [Annotations](docs/annotations.md) — VCF INFO fields, FORMAT fields and filter definitions.
- [Conventions](docs/conventions.md) — WDL, Python and codebase conventions.
- [Pipeline](docs/pipeline.md) — end-to-end pipeline description covering callset generation, preprocessing and annotation.
- [Scratch](docs/scratch.md) — processing notes and ad-hoc run logs.
