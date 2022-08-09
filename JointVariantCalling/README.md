# GATK Joint Variant Calling pipeline

## 1. Description

Pipeline to perfform GATK joint variant calling for larger cohorts. 
RecalibrationPipeline.nf can be used to run BQSR (Base Quality Score Recalibrator) if required.
Pipeline.nf is the main pipline that performs joint cariant calling. The pipleine can be run in parts (HaplotypeCaller, GenomicDBImport, GenotypeGVCFs, VariantFilter) to avoid submitting jobs spanning weeks. 

## 2. Execute

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/ExomePlus.git
   ```
2. Traverse to the folder `JointVariantCalling`


3. Modify `nextflow.config` configuration file:
    * `params.inputFilesToRecal` -- path to the input files for reaclibration, if running RecalibrationPipeline.nf (can use wildcard expressions to select files)
    * `params.recalFolder` -- path to the folder to save the recalibrated files.
    * `params.inputFiles` -- path to the input files for Pipeline.nf (joint variant calling, can use wildcard expressions to select files)
    * `params.inputFileType` -- input file type (bam or cram).
    * `params.regions` -- path to the file having list of regions to split the genomic DB (e.g., chr1:10001-207666)
    * `params.chroms` -- path to the file having list of chromosomes to generate chromosome wise final VCF files.
    * `params.label` -- cohort name.
    * `params.intervals` -- path to the gatk scattered calling intervals files based on the reference genome build (for parallel processing).
    * `params.referenceDir` -- path to the folder containing reference genome fasta file.
    * `params.referenceGenome` -- name of the reference genome fasta without extension.
    * `params.gatkContainer` -- path to the GATK singularity image file (.sif).
    * `params.gatkBundle` -- path to the gatk bundle folder. 
    * `params.genomicsDBImportFolder` -- path to save the genomicDB.
    * `params.genomicsDB` -- path to the existing genomicDB folder in case of running genomicDBImport in "Add" mode.
   
    * `params.runHaplotypeCaller` -- set true to run HaplotypeCaller.
    * `params.runGenomicDBImport` -- set true to run GenomicDBImport.
    * `params.genomicDBImportMode` -- set the genomicDBImport mode (options "Intialize" or "Add")
    * `params.runGenotypeGVCFs` -- set true to run GenotypeGVCFs.
    * `params.runVariantFilter` -- set true to run VariantFilter.


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

