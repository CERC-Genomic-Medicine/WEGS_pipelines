#!/bin/bash
#SBATCH --job-name=variantCalling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

module load nextflow
module load singularity
module load bcftools
module load samtools
nextflow run /path/to/Pipeline.nf -w /path/to/working/directory