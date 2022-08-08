#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

inputFiles = Channel.fromPath(params.inputFiles).map { file -> file.getSimpleName() }

process MakeWorkDir{
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "1h"

    input:   
    val input_label from inputFiles

    output:
    val input_label into delly, lumpy, manta, breakdancer, breakseq, cnvnator

    """
    mkdir -p ${params.result}/prepareVCFsToMerge/${input_label}
    """



}

process Delly {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"

   input:   
   val(input_label) from delly

   output:
   val input_label into delly_pop
   
   """
   module load bcftools
   module load r
   
   bcftools view ${params.result}/delly/${input_label}.vcf.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.delly.vcf.01
   ${params.survivor} filter ${input_label}.delly.vcf.01 NA 50 10000000 0 -1 ${input_label}.delly.vcf.02
   awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.delly.vcf.02 > ${input_label}.delly.vcf.03
   Rscript ${params.scripts}/changeSampleName.R ${input_label}.delly.vcf.03 DELLY_${input_label} ${input_label}.delly.vcf.04
   Rscript ${params.scripts}/renameRefAlt.R ${input_label}.delly.vcf.04 ${input_label}.delly.vcf.05.gz
   gunzip ${input_label}.delly.vcf.05.gz
   cp ${input_label}.delly.vcf.05 ${input_label}.delly.vcf.filt
   bgzip -c ${input_label}.delly.vcf.filt > ${input_label}.delly.vcf.filt.gz
   tabix -p vcf ${input_label}.delly.vcf.filt.gz
   bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${input_label}.delly.vcf.filt.gz > ${input_label}.delly.vcf.filt.genotypes
   sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.delly.vcf.filt | grep -E '^#|SVTYPE=INS' > ${input_label}.delly.vcf.filt.DUPtoINS
   grep -E '^#|SVTYPE=DEL' ${input_label}.delly.vcf.filt > ${input_label}.delly.vcf.filt.DEL
   grep -E '^#|SVTYPE=DUP' ${input_label}.delly.vcf.filt > ${input_label}.delly.vcf.filt.DUP
   grep -E '^#|SVTYPE=INS' ${input_label}.delly.vcf.filt > ${input_label}.delly.vcf.filt.INS
   grep -E '^#|SVTYPE=INV' ${input_label}.delly.vcf.filt > ${input_label}.delly.vcf.filt.INV
   grep -E '^#|SVTYPE=BND' ${input_label}.delly.vcf.filt > ${input_label}.delly.vcf.filt.BND
   sed 's/SVTYPE=BND/SVTYPE=TRA/' ${input_label}.delly.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${input_label}.delly.vcf.filt.BNDtoTRA
   cp ${input_label}.delly.vcf.filt*  ${params.result}/prepareVCFsToMerge/${input_label}

   """
}

process Manta {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"

   input:   
   val(input_label) from manta

   output:
   val input_label into manta_pop
   
   """
   module load bcftools
   module load r
   
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   bcftools view ${params.result}/manta/${input_label}.vcf.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.manta.vcf.01
   ${params.survivor} filter ${input_label}.manta.vcf.01 NA 50 -1 0 -1 ${input_label}.manta.vcf.02
   awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.manta.vcf.02 > ${input_label}.manta.vcf.03
   Rscript ${params.scripts}/changeSampleName.R ${input_label}.manta.vcf.03 MANTA_${input_label} ${input_label}.manta.vcf.04
   Rscript ${params.scripts}/renameRefAlt.R ${input_label}.manta.vcf.04 ${input_label}.manta.vcf.05
   Rscript ${params.scripts}/renameSVIds.R ${input_label}.manta.vcf.05 ${input_label}.manta.vcf.06.gz
   gunzip ${input_label}.manta.vcf.06.gz
   cp ${input_label}.manta.vcf.06 ${input_label}.manta.vcf.filt
   bgzip -c ${input_label}.manta.vcf.filt > ${input_label}.manta.vcf.filt.gz
   tabix -p vcf ${input_label}.manta.vcf.filt.gz
   bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${input_label}.manta.vcf.filt.gz > ${input_label}.manta.vcf.filt.genotypes
   sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.manta.vcf.filt | grep -E '^#|SVTYPE=INS' > ${input_label}.manta.vcf.filt.DUPtoINS
   grep -E '^#|SVTYPE=DEL' ${input_label}.manta.vcf.filt > ${input_label}.manta.vcf.filt.DEL
   grep -E '^#|SVTYPE=DUP' ${input_label}.manta.vcf.filt > ${input_label}.manta.vcf.filt.DUP
   grep -E '^#|SVTYPE=INS' ${input_label}.manta.vcf.filt > ${input_label}.manta.vcf.filt.INS
   grep -E '^#|SVTYPE=INV' ${input_label}.manta.vcf.filt > ${input_label}.manta.vcf.filt.INV
   grep -E '^#|SVTYPE=BND' ${input_label}.manta.vcf.filt > ${input_label}.manta.vcf.filt.BND
   sed 's/SVTYPE=BND/SVTYPE=TRA/' ${input_label}.manta.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${input_label}.manta.vcf.filt.BNDtoTRA
   cp ${input_label}.manta.vcf.filt*  ${params.result}/prepareVCFsToMerge/${input_label}

   """
}

process Lumpy {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"

   input:   
   val(input_label) from lumpy

   output:
   val input_label into lumpy_pop
   
   """
   module load bcftools
   module load r
   
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   cp ${params.result}/lumpy/${input_label}.vcf.gz  ${input_label}.vcf.gz
   gunzip ${input_label}.vcf.gz
   awk -F '\t' '{if(\$0 ~ "#") print; else {if(\$7 == ".") \$7="PASS"; print} }' OFS='\t' ${input_label}.vcf > ${input_label}.vcf.PASS
   mv ${input_label}.vcf.PASS ${input_label}.vcf
   bgzip -c ${input_label}.vcf > ${input_label}.vcf.gz
   tabix -p vcf ${input_label}.vcf.gz   
   bcftools view ${input_label}.vcf.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.lumpy.vcf.01
   ${params.survivor} filter ${input_label}.lumpy.vcf.01 NA 50 -1 10000000 -1 ${input_label}.lumpy.vcf.02
   Rscript ${params.scripts}/changeSampleName.R ${input_label}.lumpy.vcf.02 LUMPY_${input_label} ${input_label}.lumpy.vcf.03
   Rscript ${params.scripts}/renameRefAlt.R ${input_label}.lumpy.vcf.03 ${input_label}.lumpy.vcf.04.gz
   gunzip ${input_label}.lumpy.vcf.04.gz
   cp ${input_label}.lumpy.vcf.04 ${input_label}.lumpy.vcf.filt
   bgzip -c ${input_label}.lumpy.vcf.filt > ${input_label}.lumpy.vcf.filt.gz
   tabix -p vcf ${input_label}.lumpy.vcf.filt.gz
   bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${input_label}.lumpy.vcf.filt.gz > ${input_label}.lumpy.vcf.filt.genotypes
   sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.lumpy.vcf.filt | grep -E '^#|SVTYPE=INS' > ${input_label}.lumpy.vcf.filt.DUPtoINS
   grep -E '^#|SVTYPE=DEL' ${input_label}.lumpy.vcf.filt > ${input_label}.lumpy.vcf.filt.DEL
   grep -E '^#|SVTYPE=DUP' ${input_label}.lumpy.vcf.filt > ${input_label}.lumpy.vcf.filt.DUP
   grep -E '^#|SVTYPE=INV' ${input_label}.lumpy.vcf.filt > ${input_label}.lumpy.vcf.filt.INV
   grep -E '^#|SVTYPE=BND' ${input_label}.lumpy.vcf.filt > ${input_label}.lumpy.vcf.filt.BND
   sed 's/SVTYPE=BND/SVTYPE=TRA/' ${input_label}.lumpy.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${input_label}.lumpy.vcf.filt.BNDtoTRA
   cp ${input_label}.lumpy.vcf.filt*  ${params.result}/prepareVCFsToMerge/${input_label}

   """
}


process Breakseq {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"

   input:   
   val(input_label) from breakseq

   output:
   val input_label into breakseq_pop
   
   """
   module load bcftools
   module load r
   
   bcftools view ${params.result}/breakseq/${input_label}.vcf.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.breakseq.vcf.01
   ${params.survivor} filter ${input_label}.breakseq.vcf.01 NA 50 -1 10000000 -1 ${input_label}.breakseq.vcf.02
   awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.breakseq.vcf.02 > ${input_label}.breakseq.vcf.03
   Rscript ${params.scripts}/changeSampleName.R ${input_label}.breakseq.vcf.03 BREAKSEQ_${input_label} ${input_label}.breakseq.vcf.04
   Rscript ${params.scripts}/renameRefAlt.R ${input_label}.breakseq.vcf.04 ${input_label}.breakseq.vcf.05.gz
   gunzip ${input_label}.breakseq.vcf.05.gz
   cp ${input_label}.breakseq.vcf.05 ${input_label}.breakseq.vcf.filt
   bgzip -c ${input_label}.breakseq.vcf.filt > ${input_label}.breakseq.vcf.filt.gz
   tabix -p vcf ${input_label}.breakseq.vcf.filt.gz
   bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${input_label}.breakseq.vcf.filt.gz > ${input_label}.breakseq.vcf.filt.genotypes
   sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.breakseq.vcf.filt | grep -E '^#|SVTYPE=INS' > ${input_label}.breakseq.vcf.filt.DUPtoINS
   grep -E '^#|SVTYPE=DEL' ${input_label}.breakseq.vcf.filt > ${input_label}.breakseq.vcf.filt.DEL
   grep -E '^#|SVTYPE=DUP' ${input_label}.breakseq.vcf.filt > ${input_label}.breakseq.vcf.filt.DUP
   grep -E '^#|SVTYPE=INS' ${input_label}.breakseq.vcf.filt > ${input_label}.breakseq.vcf.filt.INS
   grep -E '^#|SVTYPE=INV' ${input_label}.breakseq.vcf.filt > ${input_label}.breakseq.vcf.filt.INV
   grep -E '^#|SVTYPE=BND' ${input_label}.breakseq.vcf.filt > ${input_label}.breakseq.vcf.filt.BND
   grep -E '^#|SVTYPE=TRA' ${input_label}.breakseq.vcf.filt > ${input_label}.breakseq.vcf.filt.TRA
   sed 's/SVTYPE=BND/SVTYPE=TRA/' ${input_label}.breakseq.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${input_label}.breakseq.vcf.filt.BNDtoTRA
   cp ${input_label}.breakseq.vcf.filt*  ${params.result}/prepareVCFsToMerge/${input_label}

   """
}


process CNVnator {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"

   input:   
   val(input_label) from cnvnator

   output:
   val input_label into cnvnator_pop
   
   """
   module load bcftools
   module load r
   
   bcftools view ${params.result}/cnvnator/${input_label}.vcf.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.cnvnator.vcf.01
   ${params.survivor} filter ${input_label}.cnvnator.vcf.01 NA 50 -1 10000000 -1 ${input_label}.cnvnator.vcf.02
   awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.cnvnator.vcf.02 > ${input_label}.cnvnator.vcf.03
   Rscript ${params.scripts}/changeSampleName.R ${input_label}.cnvnator.vcf.03 CNVNATOR_${input_label} ${input_label}.cnvnator.vcf.04
   Rscript ${params.scripts}/renameRefAlt.R ${input_label}.cnvnator.vcf.04 ${input_label}.cnvnator.vcf.05.gz
   gunzip ${input_label}.cnvnator.vcf.05.gz
   cp ${input_label}.cnvnator.vcf.05 ${input_label}.cnvnator.vcf.filt
   bgzip -c ${input_label}.cnvnator.vcf.filt > ${input_label}.cnvnator.vcf.filt.gz
   tabix -p vcf ${input_label}.cnvnator.vcf.filt.gz
   bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${input_label}.cnvnator.vcf.filt.gz > ${input_label}.cnvnator.vcf.filt.genotypes
   sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.cnvnator.vcf.filt | grep -E '^#|SVTYPE=INS' > ${input_label}.cnvnator.vcf.filt.DUPtoINS
   grep -E '^#|SVTYPE=DEL' ${input_label}.cnvnator.vcf.filt > ${input_label}.cnvnator.vcf.filt.DEL
   grep -E '^#|SVTYPE=DUP' ${input_label}.cnvnator.vcf.filt > ${input_label}.cnvnator.vcf.filt.DUP
   grep -E '^#|SVTYPE=INS' ${input_label}.cnvnator.vcf.filt > ${input_label}.cnvnator.vcf.filt.INS
   grep -E '^#|SVTYPE=INV' ${input_label}.cnvnator.vcf.filt > ${input_label}.cnvnator.vcf.filt.INV
   grep -E '^#|SVTYPE=BND' ${input_label}.cnvnator.vcf.filt > ${input_label}.cnvnator.vcf.filt.BND
   grep -E '^#|SVTYPE=TRA' ${input_label}.cnvnator.vcf.filt > ${input_label}.cnvnator.vcf.filt.TRA
   sed 's/SVTYPE=BND/SVTYPE=TRA/' ${input_label}.cnvnator.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${input_label}.cnvnator.vcf.filt.BNDtoTRA
   cp ${input_label}.cnvnator.vcf.filt*  ${params.result}/prepareVCFsToMerge/${input_label}

   """
}

process Breakdancer {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"

   input:   
   val(input_label) from breakdancer

   output:
   val input_label into breakdancer_pop
   
   """
   module load bcftools
   module load r

   cp ${params.result}/breakdancer/${input_label}.vcf.gz  ${input_label}.vcf.gz
   gunzip ${input_label}.vcf.gz
   cat /dev/null > ${input_label}.fixColumn
   grep "^##" ${input_label}.vcf >> ${input_label}.fixColumn
   echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' >> ${input_label}.fixColumn
   echo '##FORMAT=<ID=CN,Number=1,Type=String,Description="Copy number genotype for imprecise events">' >> ${input_label}.fixColumn
   grep "^#CHROM" ${input_label}.vcf | awk -v newdata="FORMAT\t${input_label}" 'BEGIN{FS=OFS="\t"} {print \$0 OFS newdata}' >> ${input_label}.fixColumn
   grep -v "^#" ${input_label}.vcf | awk -v newdata="GT:CN\t./.:." 'BEGIN{FS=OFS="\t"} {print \$0 OFS newdata}' >> ${input_label}.fixColumn
   bgzip -c ${input_label}.fixColumn > ${input_label}.fixColumn.gz
   tabix -p vcf ${input_label}.fixColumn.gz
   bcftools view ${input_label}.fixColumn.gz -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.breakdancer.vcf.01
   ${params.survivor} filter ${input_label}.breakdancer.vcf.01 NA 50 10000000 0 -1 ${input_label}.breakdancer.vcf.02
   awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.breakdancer.vcf.02 > ${input_label}.breakdancer.vcf.03
   Rscript ${params.scripts}/changeSampleName.R ${input_label}.breakdancer.vcf.03 BREAKDANCER_${input_label} ${input_label}.breakdancer.vcf.04
   Rscript ${params.scripts}/renameRefAlt.R ${input_label}.breakdancer.vcf.04 ${input_label}.breakdancer.vcf.05.gz
   gunzip ${input_label}.breakdancer.vcf.05.gz
   cp ${input_label}.breakdancer.vcf.05 ${input_label}.breakdancer.vcf.filt
   bgzip -c ${input_label}.breakdancer.vcf.filt > ${input_label}.breakdancer.vcf.filt.gz
   tabix -p vcf ${input_label}.breakdancer.vcf.filt.gz
   bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${input_label}.breakdancer.vcf.filt.gz > ${input_label}.breakdancer.vcf.filt.genotypes
   sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.breakdancer.vcf.filt | grep -E '^#|SVTYPE=INS' > ${input_label}.breakdancer.vcf.filt.DUPtoINS
   grep -E '^#|SVTYPE=DEL' ${input_label}.breakdancer.vcf.filt > ${input_label}.breakdancer.vcf.filt.DEL
   grep -E '^#|SVTYPE=DUP' ${input_label}.breakdancer.vcf.filt > ${input_label}.breakdancer.vcf.filt.DUP
   grep -E '^#|SVTYPE=INS' ${input_label}.breakdancer.vcf.filt > ${input_label}.breakdancer.vcf.filt.INS
   grep -E '^#|SVTYPE=INV' ${input_label}.breakdancer.vcf.filt > ${input_label}.breakdancer.vcf.filt.INV
   grep -E '^#|SVTYPE=BND' ${input_label}.breakdancer.vcf.filt > ${input_label}.breakdancer.vcf.filt.BND
   grep -E '^#|SVTYPE=TRA' ${input_label}.breakdancer.vcf.filt > ${input_label}.breakdancer.vcf.filt.TRA
   sed 's/SVTYPE=BND/SVTYPE=TRA/' ${input_label}.breakdancer.vcf.filt | grep -E '^#|SVTYPE=TRA' > ${input_label}.breakdancer.vcf.filt.BNDtoTRA
   cp ${input_label}.breakdancer.vcf.filt*  ${params.result}/prepareVCFsToMerge/${input_label}

   """
}

process DellyPop {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:   
   val(input_label) from delly_pop.collect()

   """
   module load bcftools
   module load vcftools
   find ${params.result}/prepareVCFsToMerge/ -name "*delly.vcf.filt" -type f > vcfCallFiles_ALL.list
   ${params.survivor} merge vcfCallFiles_ALL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_ALL.vcf
   vcf-sort -c callsMerged_ALL.vcf > callsMerged_ALL.sorted.vcf
   bgzip -c callsMerged_ALL.sorted.vcf > callsMerged_ALL.sorted.vcf.gz
   tabix -p vcf callsMerged_ALL.sorted.vcf.gz
   bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_ALL.sorted.vcf.gz > callsMerged_ALL.sorted.clean.vcf
   bgzip -c callsMerged_ALL.sorted.clean.vcf > delly.callsMerged_ALL.sorted.clean.vcf.gz
   tabix -p vcf delly.callsMerged_ALL.sorted.clean.vcf.gz

   find ${params.result}/prepareVCFsToMerge/ -name "*delly.vcf.filt.INS" -type f > vcfCallFiles_INS.list
   ${params.survivor} merge vcfCallFiles_INS.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} delly.callsMerged_INS.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*delly.vcf.filt.DEL" -type f > vcfCallFiles_DEL.list
   ${params.survivor} merge vcfCallFiles_DEL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} delly.callsMerged_DEL.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*delly.vcf.filt.DUP" -type f > vcfCallFiles_DUP.list
   ${params.survivor} merge vcfCallFiles_DUP.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} delly.callsMerged_DUP.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*delly.vcf.filt.INV" -type f > vcfCallFiles_INV.list
   ${params.survivor} merge vcfCallFiles_INV.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} delly.callsMerged_INV.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*delly.vcf.filt.BND" -type f > vcfCallFiles_BND.list
   ${params.survivor} merge vcfCallFiles_BND.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} delly.callsMerged_BND.vcf

   cp delly.callsMerged* ${params.result}/prepareVCFsToMerge/SURVIVOR 
   """
}

process MantaPop {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:   
   val(input_label) from manta_pop.collect()

   """
   module load bcftools
   module load vcftools
   find ${params.result}/prepareVCFsToMerge/ -name "*manta.vcf.filt" -type f > vcfCallFiles_ALL.list
   ${params.survivor} merge vcfCallFiles_ALL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_ALL.vcf
   vcf-sort -c callsMerged_ALL.vcf > callsMerged_ALL.sorted.vcf
   bgzip -c callsMerged_ALL.sorted.vcf > callsMerged_ALL.sorted.vcf.gz
   tabix -p vcf callsMerged_ALL.sorted.vcf.gz
   bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_ALL.sorted.vcf.gz > callsMerged_ALL.sorted.clean.vcf
   bgzip -c callsMerged_ALL.sorted.clean.vcf > manta.callsMerged_ALL.sorted.clean.vcf.gz
   tabix -p vcf manta.callsMerged_ALL.sorted.clean.vcf.gz

   find ${params.result}/prepareVCFsToMerge/ -name "*manta.vcf.filt.INS" -type f > vcfCallFiles_INS.list
   ${params.survivor} merge vcfCallFiles_INS.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} manta.callsMerged_INS.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*manta.vcf.filt.DEL" -type f > vcfCallFiles_DEL.list
   ${params.survivor} merge vcfCallFiles_DEL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} manta.callsMerged_DEL.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*manta.vcf.filt.DUP" -type f > vcfCallFiles_DUP.list
   ${params.survivor} merge vcfCallFiles_DUP.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} manta.callsMerged_DUP.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*manta.vcf.filt.INV" -type f > vcfCallFiles_INV.list
   ${params.survivor} merge vcfCallFiles_INV.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} manta.callsMerged_INV.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*manta.vcf.filt.BND" -type f > vcfCallFiles_BND.list
   ${params.survivor} merge vcfCallFiles_BND.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} manta.callsMerged_BND.vcf

   cp manta.callsMerged* ${params.result}/prepareVCFsToMerge/SURVIVOR 
   """
}

process LumpyPop {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:   
   val(input_label) from lumpy_pop.collect()

   """
   module load bcftools
   module load vcftools
   find ${params.result}/prepareVCFsToMerge/ -name "*lumpy.vcf.filt" -type f > vcfCallFiles_ALL.list
   ${params.survivor} merge vcfCallFiles_ALL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_ALL.vcf
   vcf-sort -c callsMerged_ALL.vcf > callsMerged_ALL.sorted.vcf
   bgzip -c callsMerged_ALL.sorted.vcf > callsMerged_ALL.sorted.vcf.gz
   tabix -p vcf callsMerged_ALL.sorted.vcf.gz
   bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_ALL.sorted.vcf.gz > callsMerged_ALL.sorted.clean.vcf
   bgzip -c callsMerged_ALL.sorted.clean.vcf > lumpy.callsMerged_ALL.sorted.clean.vcf.gz
   tabix -p vcf lumpy.callsMerged_ALL.sorted.clean.vcf.gz

   find ${params.result}/prepareVCFsToMerge/ -name "*lumpy.vcf.filt.DEL" -type f > vcfCallFiles_DEL.list
   ${params.survivor} merge vcfCallFiles_DEL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} lumpy.callsMerged_DEL.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*lumpy.vcf.filt.DUP" -type f > vcfCallFiles_DUP.list
   ${params.survivor} merge vcfCallFiles_DUP.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} lumpy.callsMerged_DUP.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*lumpy.vcf.filt.INV" -type f > vcfCallFiles_INV.list
   ${params.survivor} merge vcfCallFiles_INV.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} lumpy.callsMerged_INV.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*lumpy.vcf.filt.BND" -type f > vcfCallFiles_BND.list
   ${params.survivor} merge vcfCallFiles_BND.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} lumpy.callsMerged_BND.vcf

   cp lumpy.callsMerged* ${params.result}/prepareVCFsToMerge/SURVIVOR 
   """
}

process BreakdancerPop {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "64 GB"
   time "24h"

   input:   
   val(input_label) from breakdancer_pop.collect()

   """
   module load bcftools
   module load vcftools
   module load htslib/1.11
   find ${params.result}/prepareVCFsToMerge/ -name "*breakdancer.vcf.filt" -type f > vcfCallFiles_ALL.list
   ${params.survivor} merge vcfCallFiles_ALL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_ALL.vcf
   vcf-sort -c callsMerged_ALL.vcf > callsMerged_ALL.sorted.vcf
   bgzip -c callsMerged_ALL.sorted.vcf > callsMerged_ALL.sorted.vcf.gz
   tabix -p vcf callsMerged_ALL.sorted.vcf.gz
   bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_ALL.sorted.vcf.gz > callsMerged_ALL.sorted.clean.vcf
   bgzip -c callsMerged_ALL.sorted.clean.vcf > breakdancer.callsMerged_ALL.sorted.clean.vcf.gz
   tabix -p vcf breakdancer.callsMerged_ALL.sorted.clean.vcf.gz

   find ${params.result}/prepareVCFsToMerge/ -name "*breakdancer.vcf.filt.INS" -type f > vcfCallFiles_INS.list
   ${params.survivor} merge vcfCallFiles_INS.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakdancer.callsMerged_INS.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakdancer.vcf.filt.DEL" -type f > vcfCallFiles_DEL.list
   ${params.survivor} merge vcfCallFiles_DEL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakdancer.callsMerged_DEL.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakdancer.vcf.filt.DUP" -type f > vcfCallFiles_DUP.list
   ${params.survivor} merge vcfCallFiles_DUP.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakdancer.callsMerged_DUP.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakdancer.vcf.filt.INV" -type f > vcfCallFiles_INV.list
   ${params.survivor} merge vcfCallFiles_INV.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakdancer.callsMerged_INV.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakdancer.vcf.filt.BND" -type f > vcfCallFiles_BND.list
   ${params.survivor} merge vcfCallFiles_BND.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakdancer.callsMerged_BND.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakdancer.vcf.filt.TRA" -type f > vcfCallFiles_TRA.list
   ${params.survivor} merge vcfCallFiles_TRA.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakdancer.callsMerged_TRA.vcf

   cp breakdancer.callsMerged* ${params.result}/prepareVCFsToMerge/SURVIVOR 
   """
}

process BreakseqPop {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:   
   val(input_label) from breakseq_pop.collect()

   """
   module load bcftools
   module load vcftools
   find ${params.result}/prepareVCFsToMerge/ -name "*breakseq.vcf.filt" -type f > vcfCallFiles_ALL.list
   ${params.survivor} merge vcfCallFiles_ALL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_ALL.vcf
   vcf-sort -c callsMerged_ALL.vcf > callsMerged_ALL.sorted.vcf
   bgzip -c callsMerged_ALL.sorted.vcf > callsMerged_ALL.sorted.vcf.gz
   tabix -p vcf callsMerged_ALL.sorted.vcf.gz
   bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_ALL.sorted.vcf.gz > callsMerged_ALL.sorted.clean.vcf
   bgzip -c callsMerged_ALL.sorted.clean.vcf > breakseq.callsMerged_ALL.sorted.clean.vcf.gz
   tabix -p vcf breakseq.callsMerged_ALL.sorted.clean.vcf.gz

   find ${params.result}/prepareVCFsToMerge/ -name "*breakseq.vcf.filt.INS" -type f > vcfCallFiles_INS.list
   ${params.survivor} merge vcfCallFiles_INS.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakseq.callsMerged_INS.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakseq.vcf.filt.DEL" -type f > vcfCallFiles_DEL.list
   ${params.survivor} merge vcfCallFiles_DEL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakseq.callsMerged_DEL.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakseq.vcf.filt.DUP" -type f > vcfCallFiles_DUP.list
   ${params.survivor} merge vcfCallFiles_DUP.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakseq.callsMerged_DUP.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakseq.vcf.filt.INV" -type f > vcfCallFiles_INV.list
   ${params.survivor} merge vcfCallFiles_INV.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakseq.callsMerged_INV.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakseq.vcf.filt.BND" -type f > vcfCallFiles_BND.list
   ${params.survivor} merge vcfCallFiles_BND.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakseq.callsMerged_BND.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*breakseq.vcf.filt.TRA" -type f > vcfCallFiles_TRA.list
   ${params.survivor} merge vcfCallFiles_TRA.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} breakseq.callsMerged_TRA.vcf

   cp breakseq.callsMerged* ${params.result}/prepareVCFsToMerge/SURVIVOR 
   """
}

process CNVnatorPop {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:   
   val(input_label) from cnvnator_pop.collect()

   """
   module load bcftools
   module load vcftools
   find ${params.result}/prepareVCFsToMerge/ -name "*cnvnator.vcf.filt" -type f > vcfCallFiles_ALL.list
   ${params.survivor} merge vcfCallFiles_ALL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_ALL.vcf
   vcf-sort -c callsMerged_ALL.vcf > callsMerged_ALL.sorted.vcf
   bgzip -c callsMerged_ALL.sorted.vcf > callsMerged_ALL.sorted.vcf.gz
   tabix -p vcf callsMerged_ALL.sorted.vcf.gz
   bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_ALL.sorted.vcf.gz > callsMerged_ALL.sorted.clean.vcf
   bgzip -c callsMerged_ALL.sorted.clean.vcf > cnvnator.callsMerged_ALL.sorted.clean.vcf.gz
   tabix -p vcf cnvnator.callsMerged_ALL.sorted.clean.vcf.gz

   find ${params.result}/prepareVCFsToMerge/ -name "*cnvnator.vcf.filt.INS" -type f > vcfCallFiles_INS.list
   ${params.survivor} merge vcfCallFiles_INS.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} cnvnator.callsMerged_INS.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*cnvnator.vcf.filt.DEL" -type f > vcfCallFiles_DEL.list
   ${params.survivor} merge vcfCallFiles_DEL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} cnvnator.callsMerged_DEL.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*cnvnator.vcf.filt.DUP" -type f > vcfCallFiles_DUP.list
   ${params.survivor} merge vcfCallFiles_DUP.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} cnvnator.callsMerged_DUP.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*cnvnator.vcf.filt.INV" -type f > vcfCallFiles_INV.list
   ${params.survivor} merge vcfCallFiles_INV.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} cnvnator.callsMerged_INV.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*cnvnator.vcf.filt.BND" -type f > vcfCallFiles_BND.list
   ${params.survivor} merge vcfCallFiles_BND.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} cnvnator.callsMerged_BND.vcf

   find ${params.result}/prepareVCFsToMerge/ -name "*cnvnator.vcf.filt.TRA" -type f > vcfCallFiles_TRA.list
   ${params.survivor} merge vcfCallFiles_TRA.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} cnvnator.callsMerged_TRA.vcf

   cp cnvnator.callsMerged* ${params.result}/prepareVCFsToMerge/SURVIVOR 
   """
}