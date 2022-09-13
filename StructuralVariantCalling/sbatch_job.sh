#!/bin/bash
#SBATCH --job-name=sv_call
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4096
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

module load nextflow
nextflow run /path/to/Pipeline.nf -w /path/to/working/directory