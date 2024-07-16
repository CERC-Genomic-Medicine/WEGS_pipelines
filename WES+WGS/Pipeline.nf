#!/usr/bin/env nextflow

/*
* AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2022
*/

// To make sure that sample names are consistent across the read groups and also between the WES and WGS files.
// Current `samtools merge` version can't process CRAM files, so we need to make sure to write temporary files in the BAM format.
process HarmonizeSM {
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "2h"
   errorStrategy 'retry'
   maxRetries 1
   
   input:
      tuple val(sample), path(wes, stageAs: 'WES.bam'), path(wgs, stageAs: 'WGS.bam')
   output:
   tuple val(sample), path("${wes.getBaseName()}.newheader.bam"), path("${wgs.getBaseName()}.newheader.bam")

   """
   samtools view -H ${wes} | sed "s/SM:[^\t]*/SM:${sample}/g" | samtools reheader - ${wes} | samtools view -b -o ${wes.getBaseName()}.newheader.bam
   samtools view -H ${wgs} | sed "s/SM:[^\t]*/SM:${sample}/g" | samtools reheader - ${wgs} | samtools view -b -o ${wgs.getBaseName()}.newheader.bam
   """
}

// Merge the WES and WGS files
process Merge {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"
   errorStrategy 'retry'
   maxRetries 1

   input:
   tuple val(sample), path(wes), path(wgs)

   output:
   tuple val(sample), path("${sample}.WEGS.cram"), path("${sample}.WEGS.cram.crai")

   publishDir "${params.resultFolder}", pattern: "${sample}.WEGS.cram*", mode: "copy"

   """
   samtools merge -u -o - ${wes} ${wgs} | samtools view -C -T ${params.refGenome} -o ${sample}.WEGS.cram
   samtools index ${sample}.WEGS.cram ${sample}.WEGS.cram.crai
   """
}

workflow {
	bams = Channel.from(file(params.bamsListPath).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1]), file(fields[2]) ] }
	Merge(HarmonizeSM(bams))
}
