inputFiles = Channel.empty()

if( params.inputFileType == "bam" ) {
inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getBaseName(), file, file.getBaseName() + ".bai"] }
}
else if( params.inputFileType == "cram" ) {
inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getBaseName(), file, file + ".crai"] }
}

intervals = Channel.fromPath(params.intervals).map { file -> [file.getParent().toString().split('/').last(), file]}

process HaplotypeCaller {
   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "2h"

   container "${params.gatk_container}"
   containerOptions "-B ${params.reference_dir}:/ref"

   input:   
   tuple val(interval_label), file(interval_list), val(input_label), file(input), file(index) from intervals.combine(inputFiles)

   output:
   tuple val(input_label), file("${interval_label}.${input_label}.vcf.gz"), file("${interval_label}.${input_label}.vcf.gz.tbi") into interval_vcfs

   """
   gatk --java-options -Xmx4G HaplotypeCaller --native-pair-hmm-threads 1 -G StandardAnnotation -G StandardHCAnnotation -L ${interval_list} -R /ref/Homo_sapiens.GRCh38.fa -I ${input} -O ${interval_label}.${input_label}.vcf.gz
   """
}


process Merge {
   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "1h"

   input:
   tuple val(label), file(vcfs), file(vcf_indices) from interval_vcfs.groupTuple(by: 0)

   output:
   tuple val(label), file("${label}.vcf.gz"), file("${label}.vcf.gz.tbi") into vcfs

   publishDir "${params.result_folder}", pattern: "${label}.vcf.gz*"

   """
   find . -name "*.vcf.gz" | sort > files.txt
   bcftools concat -f files.txt -Oz -o ${label}.vcf.gz
   tabix ${label}.vcf.gz
   """
}

process Filter {

   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "1h"

   input:
   tuple val(label),file(vcf), file(vcf_index) from vcfs

   output:
   tuple val(label),file("${label}.hard_filter.vcf.gz"), file("${label}.hard_filter.vcf.gz.tbi") into hard_filtered_vcfs

   publishDir "${params.result_folder}", pattern: "${label}.hard_filter.vcf.gz*"

   """
   bcftools filter -s FAIL -i '(TYPE="snp" & INFO/QD>=2 & INFO/FS<=60 & INFO/MQ>=40 & INFO/SOR<=3 & (INFO/MQRankSum="." | INFO/MQRankSum>=-12.4) & (INFO/ReadPosRankSum="." | INFO/ReadPosRankSum>=-8.0)) | (TYPE="indel" & INFO/QD>=2 & INFO/FS<=200 & INFO/SOR<=10 & (INFO/ReadPosRankSum="." | INFO/ReadPosRankSum>=-20))' ${vcf} -Oz -o ${label}.hard_filter.vcf.gz
   tabix ${label}.hard_filter.vcf.gz
   """ 

}

