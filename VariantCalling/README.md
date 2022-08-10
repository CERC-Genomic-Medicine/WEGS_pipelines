# GATK Variant Calling pipeline

## 1. Description

Pipeline to perfform GATK variant calling for individual samples. 

## 2. Execute

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/ExomePlus.git
   ```
2. Traverse to the folder `VariantCalling`


3. Modify `nextflow.config` configuration file:
    * `params.inputFiles` -- path to the input files for Pipeline.nf (can use wildcard expressions to select files)
    * `params.inputFileType` -- input file type (bam or cram).
    * `params.intervals` -- path to the gatk scattered calling intervals files based on the reference genome build (for parallel processing).
    * `params.referenceDir` -- path to the folder containing reference genome fasta file.
    * `params.referenceGenome` -- name of the reference genome fasta without extension.
    * `params.gatkContainer` -- path to the GATK singularity image file (.sif).
    * `params.bundle` -- path to the gatk bundle folder. 
    * `params.resultFolder` -- path to result folder.
    * `params.doBQSR` = set 'true' to enable recalibration


4. Run pipeline (Interactive SLURM job):
    ```
    salloc --time=12:00:00 --ntasks=1 --mem-per-cpu=16G
    module load nextflow
    module load singularity
    module load bcftools
    module load samtools
    nextflow run Pipeline.nf
    ```
8. Run pipeline (sbatch job):
    ```
    sbatch sbatch_job.sh
    ```
## 3. Troubleshoot

1. `Failed to submit process to grid scheduler for execution`
   ```
   nextflow run Pipeline.nf -resume
   ```
   
2. `libnet.so: failed to map segment from shared object`, then try to increase the amount of memory in your `salloc` or `sbatch` job.

