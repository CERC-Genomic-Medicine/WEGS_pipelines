bams = Channel.from(file(params.bams_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1]), file(fields[2])] }

process CheckSampleNameMismatch {
   cache "lenient"
   cpus 1
   memory "2 GB"
   time "1h"
   
   input:
   tuple val(sample), file('bam1.bam'),file('bam2.bam') from bams
   
   output:
   tuple val(sample), file('bam1.bam'),file('bam2.bam') into bams_to_merge

   """
   arr=(\$(samtools view -H bam1.bam | grep '^@RG' | sed "s/.*SM:\\([^\t]*\\).*/\1/g"))
   arr+=(\$(samtools view -H bam2.bam | grep '^@RG' | sed "s/.*SM:\\([^\t]*\\).*/\1/g"))
   arrU=(\$(printf "%s\n" "\${arr[@]}" | sort -u))
   if [[ \${#arrU[@]} -gt 1 ]]; then echo "Sample names SM are not unique, please update the bam headers and retry";exit 1;fi
   """
}

process MergeAndMarkDuplicates {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"
   errorStrategy "finish"
   

   input:
   tuple val(sample), file(bam1), file(bam2) from bams_to_merge

   output:
   tuple file("${sample}.dedup.bam"),file("${sample}.dedup.bai"),file("${sample}.marked_dup_metrics.txt") into bams_merged_dedup
   tuple val(sample), file(bam1), file(bam2), file("${sample}.dedup.bam") into bam_for_stats
   tuple val(sample), file("${sample}.dedup.bam") into bam_to_cram


   publishDir "${params.result_folder}", pattern: "${sample}.marked_dup_metrics.txt"

   """
   java -jar -Xmx8G -XX:ParallelGCThreads=1 $EBROOTPICARD/picard.jar MarkDuplicates CREATE_INDEX=true I=${bam1} I=${bam2} O=${sample}.dedup.bam M=${sample}.marked_dup_metrics.txt SORTING_COLLECTION_SIZE_RATIO=0.05 MAX_RECORDS_IN_RAM=100000 COMPRESSION_LEVEL=3 VALIDATION_STRINGENCY=STRICT
   """
}


process GenerateCramFromBam {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "4h"
   errorStrategy "finish"

   input:
   tuple val(sample),file(mergedBam) from bam_to_cram

   output:
   tuple val(sample),file("${sample}.dedup.cram") into cram

   publishDir "${params.result_folder}", pattern: "${sample}.dedup.cram"

   """
   samtools view -C -T ${params.ref_genome} ${mergedBam} > ${sample}.dedup.cram
   """
}

bam_for_stats.into { wes_bam_ccds; wes_bam_flanking; wgs_bam_autosomal; merged_bam_ccds; merged_bam_flanking; wes_before_merge_stats; wgs_before_merge_stats; after_merge_stats}

process GenerateStatsBeforeMergeWes {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from wes_before_merge_stats

   output:
   file "${sample}-wes.dedup.bam.stats" into wes_bef_merge_stats

   publishDir "${params.result_folder}", pattern: "${sample}-wes.dedup.bam.stats"

   """
   samtools stats ${bam1} > ${sample}-wes.dedup.bam.stats
   """
}

process GenerateStatsBeforeMergeWgs {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from wgs_before_merge_stats

   output:
   file "${sample}-wgs.dedup.bam.stats" into wgs_bef_merge_stats

   publishDir "${params.result_folder}", pattern: "${sample}-wgs.dedup.bam.stats"

   """
   samtools stats ${bam2} > ${sample}-wgs.dedup.bam.stats
   """
}

process GenerateStatsAfterMerge {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from after_merge_stats

   output:
   file "${sample}-merged.dedup.bam.stats" into aft_merge_stats

   publishDir "${params.result_folder}", pattern: "${sample}-merged.dedup.bam.stats"

   """
   samtools stats ${mergedBam} > ${sample}-merged.dedup.bam.stats
   """
}

process GenerateCustomStatsFlankingWes {
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from wes_bam_flanking

   output:
   file "${sample}-wes.flanking.depth" into wes_flanking_depth

   publishDir "${params.result_folder}", pattern: "${sample}-wes.flanking.depth"

   """
   python ${params.project_dir}/customStats.py -i ${bam1} -b ${params.ccds_flanking_file} -o ${sample}-wes.flanking.depth
   """
}

process GenerateCustomStatsCCDSWes {
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from wes_bam_ccds

   output:
   file "${sample}-wes.ccds.depth" into wes_ccds_depth

   publishDir "${params.result_folder}", pattern: "${sample}-wes.ccds.depth"

   """
   python ${params.project_dir}/customStats.py -i ${bam1} -b ${params.ccds_adjmerge_file} -o ${sample}-wes.ccds.depth
   """
}
process GenerateCustomStatsFlankingMerged {
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from merged_bam_flanking

   output:
   file "${sample}-merged.flanking.depth" into merged_flanking_depth

   publishDir "${params.result_folder}", pattern: "${sample}-merged.flanking.depth"

   """
   python ${params.project_dir}/customStats.py -i ${mergedBam} -b ${params.ccds_flanking_file} -o ${sample}-merged.flanking.depth
   """
}

process GenerateCustomStatsCCDSMerged {
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from merged_bam_ccds

   output:
   file "${sample}-merged.ccds.depth" into merged_ccds_depth

   publishDir "${params.result_folder}", pattern: "${sample}-merged.ccds.depth"


   """
   python ${params.project_dir}/customStats.py -i ${mergedBam} -b ${params.ccds_adjmerge_file} -o ${sample}-merged.ccds.depth
   """
}

process GenerateCustomStatsWgs {
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), file(bam2), file(mergedBam) from wgs_bam_autosomal

   output:
   file "${sample}-wgs.autosomal.depth" into wgs_autosomal_depth

   publishDir "${params.result_folder}", pattern: "${sample}-wgs.autosomal.depth"


   """
   python ${params.project_dir}/customStats.py -a -i ${bam2} -o ${sample}-wgs.autosomal.depth
   """
}
