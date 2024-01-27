# Genotype refinement

This pipeline uses pedigree information to refine genotype qualities following the [GATK's Genotype Refinement workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531432-Genotype-Refinement-workflow-for-germline-short-variants).
First, it applies the GATK v4 CalculateGenotypePosteriors tool to recompute genotype qualities and probabilites using the provided family information (e.g. trios), which accounts for Mendel's law.
Second, it applies the GATK v4 VariantFiltration tool to annotate individual genotypes with GQ < 20 (commonly accepted threshold).
Lastly, it applies the GATK v4 VariantAnnotator tool to annotate possible de novo mutations.

> [!CAUTION]  
> Run this pipeline after the **Genotype** sup-pipeline. Otherwise, make sure that your input VCF file(s) was/were run through the GATK's Variant Quality Score Recalibration (VQSR).

## Parallelization
If you provide chunked VCF files (e.g. by chromosome or regions), then they will be processed in parallel.
The preceding **Genotype** sub-pipeline provides already chunked VCF files.

## Dependencies/pre-requisites
* Nextflow
* apptainer
* GATK v4 apptainer/singularity image build from [GATK's docker repository](https://hub.docker.com/r/broadinstitute/gatk/) 
