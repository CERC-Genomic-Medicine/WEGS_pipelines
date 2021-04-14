beds = params.beds_list_path ? Channel.from(file(params.beds_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1])] } : Channel.empty()

process MergeAdjacentCCDSRegions {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(label), file(bed) from beds

   output:
   tuple val(label), file("${label}.adjmerge.bed") into beds_for_flanking
   file "${label}.adjmerge.bed" into beds_adj_merged

   publishDir "${params.result_folder}", pattern: "${label}.adjmerge.bed"

   """
   bedtools merge -i ${bed} > ${label}.adjmerge.bed
   """
}

process FlankCCDSRegions {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(label), file(mergedBed) from beds_for_flanking

   output:
   tuple val(label), file("${label}.adjmerge.bed"), file("${label}.adjmerge.flanking.bed") into beds_for_stats
   file "${label}.adjmerge.flanking.bed" into beds_flanking

   publishDir "${params.result_folder}", pattern: "${label}.adjmerge.flanking.bed"

   """
   bedtools flank -i ${mergedBed} -g ${params.bed_genome_path} -b 10 > ${label}.adjmerge.flanking.bed
   """
}


