#!/usr/bin/env nextflow

/*
* AUTHOR: Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 1.0
* YEAR: 2024
*/

process Refine {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1 
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"

   input:
   tuple val(prefix), path(vcf), path(vcf_index)
   each path(pedigree)

   output:
   tuple path("${prefix}.gtrefined.vcf.gz"), path("${prefix}.gtrefined.vcf.gz.tbi") 

   storeDir "Refined_VCFs/"

   """
   # we set 'Xmx' argument to maximal possibe and limit the memory use through the process's 'memory' setting
   gatk --java-options "-Xmx16G" CalculateGenotypePosteriors -V ${vcf} -ped ${pedigree} --skip-population-priors -O temp.refined.vcf.gz
   gatk --java-options "-Xmx16G" VariantFiltration -V temp.refined.vcf.gz  -O temp.refined.lowGQ.vcf.gz -G-filter "GQ < 20" -G-filter-name lowGQ
   gatk --java-options "-Xmx16G" VariantAnnotator -OVI true -V temp.refined.lowGQ.vcf.gz -A PossibleDeNovo -ped ${pedigree} -O ${prefix}.gtrefined.vcf.gz
   """
}


workflow {
   // We make sure that index file (i.e. tbi or csi) is the last item in the list
   input_vcfs = Channel.fromFilePairs(params.inputVcfs, flat: true, size: 2).map(it -> ((it[1].getExtension() == "gz") || (it[1].getExtension() == "bcf")) ? [it[0], it[1], it[2]] : [it[0], it[2], it[1]])
   Refine(input_vcfs, Channel.fromPath(params.pedigree))
}
