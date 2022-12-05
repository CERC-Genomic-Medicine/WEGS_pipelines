#!/usr/bin/env nextflow

/*
* AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 3.0
* YEAR: 2021
*/


process HappyExome {  
   cpus 1
   memory "8 GB"
   time "4h"
   errorStrategy 'retry'
   maxRetries 3

   beforeScript "source ${params.virtualenv}"

   input:
   tuple val(name), path(study_vcf), path(study_target), path(truth_vcf), path(truth_high_confidence_regions)
   env HGREF

   output:
   path "*.happy_benchmark.*"

   publishDir "${params.resultFolder}", pattern: "*.happy_benchmark.*", mode: "copy"

   """
   python ${params.happy} --threads 1 ${truth_vcf} ${study_vcf} -f ${truth_high_confidence_regions} -T ${study_target} -r ${params.reference} -o ${study_vcf.getName().toString().replace(".vcf.gz", "")}.happy_benchmark
   """
}


process HappyGenome {  
   cpus 1
   memory "8 GB"
   time "4h"
   errorStrategy 'retry'
   maxRetries 3

   beforeScript "source ${params.virtualenv}"

   input:
   tuple val(name), path(study_vcf), path(truth_vcf), path(truth_high_confidence_regions)
   env HGREF

   output:
   path "*.happy_benchmark.*"

   publishDir "${params.resultFolder}", pattern: "*.happy_benchmark.*", mode: "copy"

   """
   python ${params.happy} --threads 1 ${truth_vcf} ${study_vcf} -f ${truth_high_confidence_regions} -r ${params.reference} -o ${study_vcf.getName().toString().replace(".vcf.gz", "")}.happy_benchmark
   """
}


workflow {
	truth = Channel.from(file(params.truthFiles).readLines()).map { line1 -> fields1 = line1.split(); [ fields1[0], file(fields1[1]), file(fields1[2]) ] }
	
	if (params.exome == true) {
		input = Channel.from(file(params.inputFiles).readLines()).map { line2 -> fields2 = line2.split(); [ fields2[0], file(fields2[1]), file(fields2[2]) ] }
		HappyExome(input.combine(truth, by: 0), params.reference)
	} else {
		input = Channel.from(file(params.inputFiles).readLines()).map { line2 -> fields2 = line2.split(); [ fields2[0], file(fields2[1]) ] }
		HappyGenome(input.combine(truth, by: 0), params.reference)
	}
}
