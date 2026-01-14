# Genomic characterization of a suspected hereditary cancer cohort

This repository contains the scripts developed for my Master's Thesis (TFM), focused on the genomic characterization of a cohort with suspected hereditary cancer.

The project integrates clinical data curation, somatic variant processing, Tumor Mutational Burden (TMB) calculation, visualization of mutational features and statistical analyses to describe the cohort at a genomic level.

## Analysis workflow

Somatic variant data are processed following this general workflow:

VCF --> MAF --> TMB calculation --> Visualization and statistical analysis

![workflow](./images/workflow.png)

No variant calling was performed as part of this project.

## Repository structure

#### clinical_data/

Scripts for cleaning, modifying and enriching the clinical data table.

#### vcf2maf/

Pipeline to convert somatic VCF files into MAF format using `vcf2maf`

#### TMB/

Scripts related to Tumor Mutational Burden analysis:
- Filtering of non-synonymous mutations
- TMB calculation pipelines (maftools and mutscape)
- Manual corrections and annotations
- Statistical analysis of TMB values

#### plots_scripts/

Scripts for graphical representations, including:
- Cohort overview plots
- Tumor type and subtype distribution
- Age-related plots
- TMB boxplots per method and tumor type
- Visualization of provided mutational signatures
- Oncoplots (global and focused on specific gene sets)
- Somatic interaction plots
- Variant Allele Frequency (VAF) plots

#### environments/

- environment.yml: minimal environment
- environment_full.yml: full environment with all dependencies

Conda environments are provided to ensure reproducibility of the analyses.
