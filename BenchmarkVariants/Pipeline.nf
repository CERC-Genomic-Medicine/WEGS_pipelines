truth = Channel.from(file(params.truth_files).readLines()).map { line1 -> fields1 = line1.split(); [ fields1[0], file(fields1[1]), file(fields1[2]) ] }
input = Channel.from(file(params.input_files).readLines()).map { line2 -> fields2 = line2.split(); [ fields2[0], file(fields2[1]), file(fields2[2]) ] }

process Happy {  
   cpus 1
   memory "8 GB"
   time "2h"

   beforeScript "source ${params.virtualenv}"

   input:
   tuple val(name), file(study_vcf), file(study_target), file(truth_vcf), file(truth_high_confidence_regions) from input.combine(truth, by: 0)
   env HGREF from params.reference

   output:
   file "${study_vcf.getSimpleName()}.*"

   publishDir "${params.result_folder}", pattern: "${study_vcf.getSimpleName()}.*" 

   """
   #source ${params.virtualenv}
   python ${params.happy} --threads 1 ${truth_vcf} ${study_vcf} -f ${truth_high_confidence_regions} -T ${study_target} -r ${params.reference} -o ${study_vcf.getSimpleName()}
   """
}