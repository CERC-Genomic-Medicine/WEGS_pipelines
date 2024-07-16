# Pipeline to merge WES and WGS BAM/CRAM files

This pipeline uses samtools to merge sequencing reads generated for the same sample in WES and WGS experiments.
It expects WES reads and WGS reads to be stored in two separate single sample BAM/CRAM files.
Before merging these two BAM/CRAM files, the pipeline will replace all `SM` tag values with the provided sample name.

> [!CAUTION]
> The `SM` tag values, which store the sample name, can be different in the input WES and WGS BAM/CRAM files. 
> The pipeline doesn't do any checks to see if reads really belong to the same individual.
> **Thus, it is very important, that you correctly pair WES and WGS files based on your manifest file from sequencing experiments.**

> [!IMPORTANT]
> The pipeline needs to produce two temporary BAM files - one for WES and another for WGS data - with the updated headers. Thus, depending on your storage capacity, you may want to process your data in batches.

## Steps

1. Clone this repository to the directory where you will run the pipeline:
   ```
   git clone https://github.com/CERC-Genomic-Medicine/WEGS_pipelines.git
   ```
   
2. Go to `WEGS_pipelines/WES+WGS` directory.
3. Create a file with the list of paired WES and WGS files to merge. The file must have no header and must contain one line per sample like in the following example:
   ```
   sample1 /path/to/wes/sample1.wes.bam /path/to/wgs/sample1.wgs.bam
   sample2 /path/to/wes/sample2.wes.bam /path/to/wgs/sample2.wgs.bam
   ...
   sampleN /path/to/wes/sampleN.wes.bam /path/to/wgs/sampleN.wgs.bam
   ```
> [!IMPORTANT]
> The BAM/CRAM filenames for WES and WGS data can be identical as long as their full paths are different. Or, WES and WGS data can be located on the same path as long as their BAM/CRAM filenames are different.
> The specified sample name will be used to name the final CRAM files with the merged reads and in the corresponding `SM` tags.
   
4. Modify the `nextflow.config` configuration file as needed.
5. Run `nextflow run Pipeline.nf`.

## Dependencies/pre-requisites
1. samtools
2. Nextflow
3. Human genome reference file (i.e. *.fa and associated indices)
