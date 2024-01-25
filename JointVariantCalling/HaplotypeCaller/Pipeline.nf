#!/usr/bin/env nextflow

/*
* AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2022
*/


process HaplotypeCaller {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input:   
   tuple path(interval_list), path(input), path(index) 

   output:
   path "*.${input.getSimpleName()}.g.vcf.gz"
   
   """
   gatk --java-options -Xmx14G HaplotypeCaller --ERC GVCF --native-pair-hmm-threads 1 -G StandardAnnotation -G StandardHCAnnotation -L ${interval_list} -R /ref/${params.referenceGenome} -I ${input} -O temp.g.vcf.gz

   chunk_size=10000000
   while read -r chr interval_start interval_stop; do
      for chunk_start in `seq \${interval_start} \${chunk_size} \${interval_stop}`; do  
         chunk_stop=\$((chunk_start + chunk_size - 1))
         if [ "\${chunk_stop}" -gt "\${interval_stop}" ]; then
	    chunk_stop=\${interval_stop}
         fi
         tabix -h temp.g.vcf.gz \${chr}:\${chunk_start}-\${chunk_stop} | awk -v r=\${chunk_start} -F "\t" '/^#/ { print } !/^#/ {if(\$2 >= r) {print}}' | bgzip -c > \${chr}_\${chunk_start}_\${chunk_stop}.${input.getSimpleName()}.g.vcf.gz
      done    
   done < <(grep -v "^@" ${interval_list} | cut -f1-3)
   """
}


process GenomicsDBImport {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16G"
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"
   
   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input:
   tuple val(region), path(gvcfs)

   output:
   path "${region.split('_')[1]}_${region.split('_')[2]}"

   storeDir "GenomicsDB/${region.split('_')[0]}/"

   """
   for f in *.g.vcf.gz; do id=`echo \${f} | cut -f2 -d.`; echo -e "\${id}\t\${f}"; done > sampleMap.txt
   for f in *.g.vcf.gz; do tabix \${f}; done
   region=${region.split("_")[0]}:${region.split("_")[1]}-${region.split("_")[2]}
   path=${region.split("_")[1]}_${region.split("_")[2]}
   gatk --java-options "-Xms4g -Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport -L \${region} --genomicsdb-workspace-path \${path} --batch-size 200 --sample-name-map sampleMap.txt
   """
}


process GenomicsDBUpdate {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16G"
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"
   
   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input:
   tuple val(region), path(gvcfs), path(genomicsdb)

   output:
   path("${region.split('_')[1]}_${region.split('_')[2]}"), includeInputs: true

   storeDir "GenomicsDB/${region.split('_')[0]}/"

   """
   for f in *.g.vcf.gz; do id=`echo \${f} | cut -f2 -d.`; echo -e "\${id}\t\${f}"; done > sampleMap.txt
   for f in *.g.vcf.gz; do tabix \${f}; done
   region=${region.split("_")[0]}:${region.split("_")[1]}-${region.split("_")[2]}
   path=${region.split("_")[1]}_${region.split("_")[2]}
   gatk --java-options "-Xms4g -Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport -L \${region} --genomicsdb-update-workspace-path \${path} --batch-size 200 --sample-name-map sampleMap.txt
   """
}


workflow {
   intervals = Channel.fromPath(params.intervals)
	
   // Pair *.bam and *.bai (or *.cram and *.crai) files by their prefix, and make sure that *.bam (or *.cram) file is first in the emmited pair 
   bams = Channel.fromFilePairs(params.inputFiles, flat: true, size: 2).map(it -> ((it[1].getExtension() == "bam") || (it[1].getExtension() == "cram")) ? [it[1], it[2]] : [it[2], it[1]])
   chunked_gvcfs = HaplotypeCaller(intervals.combine(bams)).flatten().map{ file -> [file.getSimpleName(), file]}.groupTuple(by: 0, sort: "hash")
   if (params.genomicsDB) {
      chunked_genomicsDB = Channel.fromPath(params.genomicsDB, type: "dir").map{ file-> [ file.getParent().toString().split('/').last() + "_" + file.getName(), file] }
      GenomicsDBUpdate(chunked_gvcfs.join(chunked_genomicsDB))
   } else {
	GenomicsDBImport(chunked_gvcfs)
   }	
}
