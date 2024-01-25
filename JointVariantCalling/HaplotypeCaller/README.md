# Single-sample variant calling

This step runs GATK v4 HaplotypeCaller on each individual (i.e. BAM/CRAM file) and saves the called individual genotypes in the GATK v4 GenomicsDB data store, which can be used for subsequent joint genotyping.

## Assumptions

We assume, that the inividual BAM/CRAM files were already passed through the tools for marking duplicate reads and recallibrating base qualities.

## Parallelization

The pipeline distributes computations by splitting each individual genome into the [GATK's pre-defined intervals](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) and analyzing each interval in-parallel. This technique requires allocation of many small machines (1CPU, ~16GB RAM), which may be more advantegeous than allocating large machines on some systems. It also provides finer granularity for resuming the computations. 

## Dependencies/pre-requisites
* Nextflow
* apptainer
* GATK 4 apptainer/singularity image build from [GATK's docker repository](https://hub.docker.com/r/broadinstitute/gatk/) 
* [GATK 4 resource bundle and scattered calling intervals](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
* Human genome reference file (i.e. *.fa and associated indices)
