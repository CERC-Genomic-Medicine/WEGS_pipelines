#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/


path = params.bamsFolder + "/*bam"
Channel.fromPath(path).map { file -> [file.getSimpleName(), file, file + ".bai"] }.into { del; dup; inv }

process Del {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:
   tuple val(id), file(bam), file(idx) from del

   """
   module load samtools
   module load bcftools
   module load python
   export PATH="\${PATH}:/home/praveen/projects/rrg-vmooser/praveen/tools/SV/lumpy-sv/bin:/home/praveen/projects/rrg-vmooser/praveen/tools/smoove"
   ${params.smoove} genotype -p 1 --name $id --outdir ${params.result}/genotype/DELs/$id --fasta ${params.referenceGenome} --duphold --vcf ${params.result}/genotype/samples_merged_DEL.vcf $bam
   """
}

process Dup {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:
   tuple val(id), file(bam), file(idx) from dup

   """
   module load samtools
   module load bcftools
   module load python
   export PATH="\${PATH}:/home/praveen/projects/rrg-vmooser/praveen/tools/SV/lumpy-sv/bin:/home/praveen/projects/rrg-vmooser/praveen/tools/smoove"
   ${params.smoove} genotype -p 1 --name $id --outdir ${params.result}/genotype/DUPs/$id --fasta ${params.referenceGenome} --duphold --vcf ${params.result}/genotype/samples_merged_DUP.vcf $bam
   """
}

process Inv {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:
   tuple val(id), file(bam), file(idx) from inv

   """
   module load samtools
   module load bcftools
   module load python
   export PATH="\${PATH}:/home/praveen/projects/rrg-vmooser/praveen/tools/SV/lumpy-sv/bin:/home/praveen/projects/rrg-vmooser/praveen/tools/smoove"
   ${params.smoove} genotype -p 1 --name $id --outdir ${params.result}/genotype/INVs/$id --fasta ${params.referenceGenome} --duphold --vcf ${params.result}/genotype/samples_merged_INV.vcf $bam
   """
}

