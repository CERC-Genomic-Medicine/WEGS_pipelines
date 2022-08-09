#!/bin/bash
#SBATCH --job-name=jvc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4096
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

module load nextflow
module load singularity
module load bcftools
module load samtools
nextflow run /path/to/Pipeline.nf