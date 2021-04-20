#!/bin/bash
#SBATCH --job-name=benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

module load nextflow
nextflow run /home/praveen/projects/def-vmooser/praveen/ExomePlus/BenchmarkVariants/Pipeline.nf -w /home/praveen/scratch/wes+wgs