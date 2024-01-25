# Single-sample variant calling

This step runs GATK v4 HaplotypeCaller on each individual (i.e. BAM/CRAM file) and saves the called individual genotypes in the GATK v4 GenomicsDB data store, which can be used for subsequent joint genotyping.

> [!TIP]
> You can update an existing GenomicsDB by adding additional individuals as long as the existing GenomicsDB was created using the same pipeline.

## Assumptions

We assume, that the inividual BAM/CRAM files were already passed through the tools for marking duplicate reads and recallibrating base qualities.

## Incrementally adding new sequenced individuals

If you are performing variant calling on a first sequencing batch of individuals in your study or processing all individuals at once, then leave the ```genomicsDB``` configuration parameter empty.

If you want to add the next sequencing batch of individuals to the existing variant calls, then set the ```genomicsDB``` configuration parameter to the GenomicsDB directory from your previous run. To generate a new multi-sample VCF, run the **Genotype** sub-pipeline on the updated GenomicsDB directory.

> [!CAUTION]
> The incrementally adding new sequenced individuals will only work if the existing GenomicsDB directory was generated using the same **HaplotypeCaller** pipeline and ```intervals```. 

## Parallelization

The pipeline distributes computations by splitting each individual genome into the [GATK's pre-defined intervals](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) and analyzing each interval in-parallel. This technique requires allocation of many small machines (1CPU, ~16GB RAM), which may be more advantegeous than allocating large machines on some systems. It also provides finer granularity for resuming the computations. 

## Dependencies/pre-requisites
* Nextflow
* apptainer
* GATK v4 apptainer/singularity image build from [GATK's docker repository](https://hub.docker.com/r/broadinstitute/gatk/) 
* [GATK v4 resource bundle and scattered calling intervals](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
* Human genome reference file (i.e. *.fa and associated indices)
