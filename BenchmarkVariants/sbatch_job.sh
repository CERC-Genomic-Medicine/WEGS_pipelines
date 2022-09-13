#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00

module load nextflow
nextflow run /path/to/Pipeline.nf -w /path/to/working/directory