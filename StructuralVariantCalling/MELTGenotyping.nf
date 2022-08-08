#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

inputFiles = Channel.fromPath(params.bamsFolder).map { file -> [file, file + ".bai"] }.into { melt_alu; melt_hervk; melt_line1; melt_sva }

process MeltALUGroup{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "72h"
  
   output:
   val true into meltALUGroup

   """
   module load bowtie2
   module load java
   module load samtools
   

   java -Xmx12G -jar ${params.melt}/MELT.jar GroupAnalysis \
    -discoverydir ${params.result}/melt/ALU/ \
    -w ${params.result}/melt/ALU/ \
    -t ${params.mei}/ALU_MELT.zip \
    -h ${params.referenceGenome} \
    -n ${params.melt}/add_bed_files/Hg38/Hg38.genes.bed
   """

}

process MeltHERVKGroup{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "72h"

   output:
   val true into meltHERVKGroup

   """
   module load bowtie2
   module load java
   module load samtools
   

   java -Xmx12G -jar ${params.melt}/MELT.jar GroupAnalysis \
    -discoverydir ${params.result}/melt/HERVK/ \
    -w ${params.result}/melt/HERVK/ \
    -t ${params.mei}/HERVK_MELT.zip \
    -h ${params.referenceGenome} \
    -n ${params.melt}/add_bed_files/Hg38/Hg38.genes.bed
   """

}

process MeltLINE1Group{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "72h"

   output:
   val true into meltLINE1Group

   """
   module load bowtie2
   module load java
   module load samtools
   

   java -Xmx12G -jar ${params.melt}/MELT.jar GroupAnalysis \
    -discoverydir ${params.result}/melt/LINE1/ \
    -w ${params.result}/melt/LINE1/ \
    -t ${params.mei}/LINE1_MELT.zip \
    -h ${params.referenceGenome} \
    -n ${params.melt}/add_bed_files/Hg38/Hg38.genes.bed
   """

}


process MeltSVAGroup{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "72h"

   output:
   val true into meltSVAGroup

   """
   module load bowtie2
   module load java
   module load samtools
   

   java -Xmx12G -jar ${params.melt}/MELT.jar GroupAnalysis \
    -discoverydir ${params.result}/melt/SVA/ \
    -w ${params.result}/melt/SVA/ \
    -t ${params.mei}/SVA_MELT.zip \
    -h ${params.referenceGenome} \
    -n ${params.melt}/add_bed_files/Hg38/Hg38.genes.bed
   """

}

process MeltALU{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:
   val(flag) from meltALUGroup
   tuple file(inputFile), file(index) from melt_alu

   output:
   val true into makeVcfALU

   """
   module load bowtie2
   module load java
   module load samtools

   java -Xmx12G -jar ${params.melt}/MELT.jar Genotype \
    -bamfile ${inputFile} \
    -t ${params.mei}/ALU_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/ALU/ \
    -p ${params.result}/melt/ALU/
   """
}

process MeltHERVK{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:
   val(flag) from meltHERVKGroup   
   tuple file(inputFile), file(index) from melt_hervk

   output:
   val true into makeVcfHERVK

   """
   module load bowtie2
   module load java
   module load samtools

   java -Xmx12G -jar ${params.melt}/MELT.jar Genotype \
    -bamfile ${inputFile} \
    -t ${params.mei}/HERVK_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/HERVK/ \
    -p ${params.result}/melt/HERVK/
   """
}

process MeltLINE1{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:   
   val(flag) from meltLINE1Group
   tuple file(inputFile), file(index) from melt_line1

   output:
   val true into makeVcfLINE1

   """
   module load bowtie2
   module load java
   module load samtools

   java -Xmx12G -jar ${params.melt}/MELT.jar Genotype \
    -bamfile ${inputFile} \
    -t ${params.mei}/LINE1_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/LINE1/ \
    -p ${params.result}/melt/LINE1/
   """
}

process MeltSVA{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:
   val(flag) from meltSVAGroup   
   tuple file(inputFile), file(index) from melt_sva

   output:
   val true into makeVcfSVA

   """
   module load bowtie2
   module load java
   module load samtools

   java -Xmx12G -jar ${params.melt}/MELT.jar Genotype \
    -bamfile ${inputFile} \
    -t ${params.mei}/SVA_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/SVA/ \
    -p ${params.result}/melt/SVA/
   """
}


process MeltALUVcf{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:
   val(flag) from makeVcfALU.collect() 

   """
   module load bowtie2
   module load java
   module load samtools
   
   java -Xmx12G -jar ${params.melt}/MELT.jar MakeVCF \
    -genotypingdir ${params.result}/melt/ALU/ \
    -t ${params.mei}/ALU_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/ALU/ \
    -p ${params.result}/melt/ALU/ \
    -o ${params.result}/melt/ALU/
   """
}

process MeltHERVKVcf{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:
   val(flag) from makeVcfHERVK.collect() 

   """
   module load bowtie2
   module load java
   module load samtools

   java -Xmx12G -jar ${params.melt}/MELT.jar MakeVCF \
    -genotypingdir ${params.result}/melt/HERVK/ \
    -t ${params.mei}/HERVK_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/HERVK/ \
    -p ${params.result}/melt/HERVK/ \
    -o ${params.result}/melt/HERVK/
   """
}

process MeltLINE1Vcf{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:
   val(flag) from makeVcfLINE1.collect() 

   """
   module load bowtie2
   module load java
   module load samtools

   java -Xmx12G -jar ${params.melt}/MELT.jar MakeVCF \
    -genotypingdir ${params.result}/melt/LINE1/ \
    -t ${params.mei}/LINE1_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/LINE1/ \
    -p ${params.result}/melt/LINE1/ \
    -o ${params.result}/melt/LINE1/
   """
}

process MeltSVAVcf{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   input:
   val(flag) from makeVcfSVA.collect() 

   """
   module load bowtie2
   module load java
   module load samtools

   java -Xmx12G -jar ${params.melt}/MELT.jar MakeVCF \
    -genotypingdir ${params.result}/melt/SVA/ \
    -t ${params.mei}/SVA_MELT.zip \
    -h ${params.referenceGenome} \
    -w ${params.result}/melt/SVA/ \
    -p ${params.result}/melt/SVA/ \
    -o ${params.result}/melt/SVA/
   """
}