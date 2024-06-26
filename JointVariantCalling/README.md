# Joint variant calling for SNVs and short indels

This pipeline performs joint variant calling of SNVs (single-nucleotide variations) and short indels (insertions/deletions) in large cohorts using GATK v4 and following the [GATK's Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels). For execution tracking convenience, the pipeline is split into multiple sub-pipelines, namely: **Recalibration**, **HaplotypeCaller**, **Genotype**, and **Refine**.

> [!TIP]
> The pipeline supports appending additionally sequenced individuals. See ```genomicsDB``` configuration parameter inside the **Genotype** sub-pipeline.

> [!TIP]
> The pipeline can be applied to small studies of related individuals (e.g. trios). See ```pedigree``` configuration parameter inside the **Genotype** sub-pipeline.


<p align="center">
<img src="JointVariantCalling.png" width="60%"/>
</p>

## Assumptions

We assume, that the inividual BAM/CRAM files were already passed through the tools for marking duplicate reads. If your BAM/CRAM files were also passed through the tools for base quality recalibration, then you can skip the **Recalibration** sub-pipeline.

> [!TIP]
> It is always a good idea to run your BAM/CRAM files through the [VerifyBamID pipeline](https://github.com/CERC-Genomic-Medicine/verifybam_pipeline) before the variant calling, which verifies sequencing depth, checks for possible contamination, checks for missing chromosomes, estimates genetic ancestry, and confirms ploidy of sex chromosomes against the reported sex values.

## Steps

> [!IMPORTANT]  
> If the pipeline failed, and you want to resume it from where it stopped (e.g. after fixing the issue), add **-resume** flag to the Nextflow command.

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/WEGS_pipelines.git
   ```

2. (Optional) **Recalibration**. If the base qualities in your BAM/CRAM files where not yet recalibrated (check your header information ```samtools view -H <BAM/CRAM>``` or ask a person who generated the files), then you need to run this step before variant calling.
   1. Go to ```WEGS_pipelines/JointVariantCalling/Recalibration``` directory.
   2. Modify the ```nextflow.config``` configuration file as needed. See README.md inside this directory for more details.
   3. Run ```nextflow run Pipeline.nf```.

3. **HaplotypeCaller**. In this step, you will generate genotype calls for each individual in your dataset separately.
   1. Go to ```WEGS_pipelines/JointVariantCalling/HaplotypeCaller``` directory.
   2. Modify the ```nextflow.config``` configuration file as needed. If you did **Recalibration** step, then you should use BAM/CRAM files produced by it as an input to this step. See README.md inside this directory for more details.
   3. Run ```nextflow run Pipeline.nf```.
  
4. **Genotype**. The previous step (i.e. **HaplotypeCaller**) produced the GenomicsDB data store, which stores all genotypes for all individuals. In this step, you will use GenomicsDB to perform joint genotyping, which will refine individual genotype calls based on observed genotypes across all individuals.
   1. Go to ```WEGS_pipelines/JointVariantCalling/Genotype``` directory.
   2. Modify the ```nextflow.config``` configuration file as needed. The input for this step is the GenomicsDB output from **HaplotypeCaller**. See README.md inside this directory for more details.
   3. Run ```nextflow run Pipeline.nf```.
   4. The output will be inside the ```VCFs``` folder, chunked into intervals (chunking is usefull if you want to parallelize any downstream analyses). If you want to merge all the results into a single large VCF file, then use:
      ```
      find VCFs/ -name "*.vcf.gz" | sort -V > files.txt
      bcftools concat -f files.txt -Oz -o all.vcf.gz
      ```
5. (Optional) **Refine**. If you have a small study of related individuals (e.g. trios), then you may want to identify Mendelian errors and flag de novo mutations.
   1. Go to ```WEGS_pipelines/JointVariantCalling/Refine``` directory.
   2. Modify the ```nextflow.config`` configuration file as needed.
   3. Run ```nextflow run Pipeline.nf```.
   4. The output will be inside the ```Refined_VCFs``` folder, chunked into intervals. If you want to merge all the results into a single large VCF file, then use:
      ```
      find Refined_VCFs/ -name "*.vcf.gz" | sort -V > files.txt
      bcftools concat -f files.txt -Oz -o all.vcf.gz
      ```
   
## Dependencies/pre-requisites

Within each sub-pipeline directory, you will find a list of dependencies/pre-requisites specific for that step. Here, we provide an aggregated list:
* Nextflow
* apptainer
* GATK v4 apptainer/singularity image build from [GATK's docker repository](https://hub.docker.com/r/broadinstitute/gatk/) 
* [GATK v4 resource bundle and scattered calling intervals](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
* Human genome reference file (i.e. *.fa and associated indices)
* bcftools
* tabix
* samtools

## Job scheduler

The pipeline was designed assuming that you are using the SLURM job scheduler.
Since the execution of the some sub-pipelines may take days (depending on the size of your data), we recommend launching them as SLURM jobs, e.g.
```
sbatch --time=5-0 --mem=8G --cpus-per-task=1 --wrap="nextflow run Pipeline.nf" -o pipeline_run.log
```

Depending on your system, you may need to pre-load the required tools first, e.g.:
```
module load nextflow
module load apptainer
module load bcftools
module load samtools
```

## Troubleshoot
* `Failed to submit process to grid scheduler for execution`
   ```
   nextflow run Pipeline.nf -resume
   ```
   
* `libnet.so: failed to map segment from shared object`, then try to increase the amount of memory in your `salloc` or `sbatch` job.

