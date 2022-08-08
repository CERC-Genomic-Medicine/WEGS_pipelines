#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

process MakeResultDirs{
    label "PreProcessing"
    cache "lenient"
    executor "local"
    cpus 1

    output:
    val 1 into del, dup, inv

    """
    mkdir -p ${params.result}/genotype/Results
    """
}

process Del {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:
   val(flag) from del

   """
   module load samtools
   module load bcftools
   module load python
   module load r
   export PATH="\${PATH}:/home/praveen/projects/rrg-vmooser/praveen/tools/SV/lumpy-sv/bin:/home/praveen/projects/rrg-vmooser/praveen/tools/smoove"
   ${params.smoove} paste --name DEL --outdir ${params.result}/genotype/Results ${params.result}/genotype/DELs/*/*.genotyped.vcf.gz
   """
}

process Dup {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:
   val(flag) from dup

   """
   module load samtools
   module load bcftools
   module load python
   module load r
   export PATH="\${PATH}:/home/praveen/projects/rrg-vmooser/praveen/tools/SV/lumpy-sv/bin:/home/praveen/projects/rrg-vmooser/praveen/tools/smoove"
   ${params.smoove} paste --name DUP --outdir ${params.result}/genotype/Results ${params.result}/genotype/DUPs/*/*.genotyped.vcf.gz
   """
}

process Inv {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:
   val(flag) from inv

   """
   module load samtools
   module load bcftools
   module load python
   module load r
   export PATH="\${PATH}:/home/praveen/projects/rrg-vmooser/praveen/tools/SV/lumpy-sv/bin:/home/praveen/projects/rrg-vmooser/praveen/tools/smoove"
   ${params.smoove} paste --name INV --outdir ${params.result}/genotype/Results ${params.result}/genotype/INVs/*/*.genotyped.vcf.gz
   """
}





