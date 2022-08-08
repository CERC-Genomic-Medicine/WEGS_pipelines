#!/usr/bin/env Rscript

#################################################################################

# Copyright (c) 2019 Ricardo Vialle and Raj Laboratory
# Contact: Ricardo Vialle <ricardovialle@gmail.com>

#################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Usage: Rscript harmonize_insertions cohort_root_dir output_vcf", call.=FALSE)
} 

work_dir = args[1]
output_vcf = args[2]

library(vcfR)
library(GenomicRanges)
library(regioneR)

getoverlap <- function(ref,test){
  refGR <- toGRanges(ref[,1:5])
  testGR <- toGRanges(test[,1:5])
  
  hits = findOverlaps(refGR, testGR)
  overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
  refPercentOverlap <- width(overlaps) / width(refGR[queryHits(hits)])
  testPercentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
  
  hits <- hits[(refPercentOverlap > 0.7) & (testPercentOverlap > 0.7)]
  return(hits)
}

## INS from other callers
survivor_dir = paste0(work_dir, "/mergeSamples")
ins = read.vcfR(paste0(survivor_dir,"/samples_merged_INS.vcf"), verbose = F)

ins_bed = as.data.frame(cbind(ins@fix[,c("CHROM","POS","ID")], extract.info(ins, "SVLEN")), stringsAsFactors = F)
colnames(ins_bed) = c("CHROM","POS","ID","SVLEN")
ins_bed$POS = as.numeric(ins_bed$POS)
ins_bed$END = ins_bed$POS + abs(as.numeric(ins_bed$SVLEN))
ins_bed = ins_bed[,c("CHROM","POS","END","ID","SVLEN")]
ins_bed = ins_bed[ order(ins_bed[,1],ins_bed[,2],ins_bed[,3]) , ]

## INS from MELT
melt_dir = paste0(work_dir, "/melt")
alu = read.vcfR(paste0(melt_dir,"/ALU/ALU.final_comp.vcf"), verbose = F)
line1 = read.vcfR(paste0(melt_dir,"/LINE1/LINE1.final_comp.vcf"), verbose = F)
sva = read.vcfR(paste0(melt_dir,"/SVA/SVA.final_comp.vcf"), verbose = F)
hervk = read.vcfR(paste0(melt_dir,"/HERVK/HERVK.final_comp.vcf"), verbose = F)

names_alu = colnames(alu@gt)
names_line1 = colnames(line1@gt)
names_sva = colnames(sva@gt)
names_hervk = colnames(hervk@gt)

colnames(alu@gt) = gsub("(.*?)\\.(.*)","\\1",colnames(alu@gt))
colnames(line1@gt) = gsub("(.*?)\\.(.*)","\\1",colnames(line1@gt))
colnames(sva@gt) = gsub("(.*?)\\.(.*)","\\1",colnames(sva@gt))
colnames(hervk@gt) = gsub("(.*?)\\.(.*)","\\1",colnames(hervk@gt))

alu.pass = alu[alu@fix[,7]=="PASS"]
line1.pass = line1[line1@fix[,7]=="PASS"]
sva.pass = sva[sva@fix[,7]=="PASS"]
hervk.pass = hervk[hervk@fix[,7]=="PASS"]

alu_bed = as.data.frame(cbind(alu.pass@fix[,c("CHROM","POS","ALT")], extract.info(alu.pass, "SVLEN")), stringsAsFactors = F)
colnames(alu_bed) = c("CHROM","POS","ID","SVLEN")
alu_bed$POS = as.numeric(alu_bed$POS)
alu_bed$END = alu_bed$POS + as.numeric(alu_bed$SVLEN)
alu_bed = alu_bed[,c("CHROM","POS","END","ID","SVLEN")]
alu_bed = alu_bed[ order(alu_bed[,1],alu_bed[,2],alu_bed[,3]) , ]

line1_bed = as.data.frame(cbind(line1.pass@fix[,c("CHROM","POS","ALT")], extract.info(line1.pass, "SVLEN")), stringsAsFactors = F)
colnames(line1_bed) = c("CHROM","POS","ID","SVLEN")
line1_bed$POS = as.numeric(line1_bed$POS)
line1_bed$END = line1_bed$POS + as.numeric(line1_bed$SVLEN)
line1_bed = line1_bed[,c("CHROM","POS","END","ID","SVLEN")]
line1_bed = line1_bed[ order(line1_bed[,1],line1_bed[,2],line1_bed[,3]) , ]

sva_bed = as.data.frame(cbind(sva.pass@fix[,c("CHROM","POS","ALT")], extract.info(sva.pass, "SVLEN")), stringsAsFactors = F)
colnames(sva_bed) = c("CHROM","POS","ID","SVLEN")
sva_bed$POS = as.numeric(sva_bed$POS)
sva_bed$END = sva_bed$POS + as.numeric(sva_bed$SVLEN)
sva_bed = sva_bed[,c("CHROM","POS","END","ID","SVLEN")]
sva_bed = sva_bed[ order(sva_bed[,1],sva_bed[,2],sva_bed[,3]) , ]

hervk_bed = as.data.frame(cbind(hervk.pass@fix[,c("CHROM","POS","ALT")], extract.info(hervk.pass, "SVLEN")), stringsAsFactors = F)
colnames(hervk_bed) = c("CHROM","POS","ID","SVLEN")
hervk_bed$POS = as.numeric(hervk_bed$POS)
hervk_bed$END = hervk_bed$POS + as.numeric(hervk_bed$SVLEN)
hervk_bed = hervk_bed[,c("CHROM","POS","END","ID","SVLEN")]
hervk_bed = hervk_bed[ order(hervk_bed[,1],hervk_bed[,2],hervk_bed[,3]) , ]

ins_alu = getoverlap(ins_bed,alu_bed)
ins_line1 = getoverlap(ins_bed,line1_bed)
ins_sva = getoverlap(ins_bed,sva_bed)
ins_hervk = getoverlap(ins_bed,hervk_bed)

# Flag IDs to remove
ins_to_remove = ins_bed[unique(c(queryHits(ins_alu),queryHits(ins_line1),queryHits(ins_sva)),queryHits(ins_hervk)),]
write.table(x = ins_to_remove$ID, file = paste0(survivor_dir,"/INS.toremove"), quote = F, row.names = F, col.names = F)

ins_clean = ins[!ins_bed$ID%in%ins_to_remove$ID]
samples = colnames(ins@gt)[-1]

alu.pass@gt = alu.pass@gt[,c("FORMAT",samples)]
line1.pass@gt = line1.pass@gt[,c("FORMAT",samples)]
sva.pass@gt = sva.pass@gt[,c("FORMAT",samples)]
hervk.pass@gt = hervk.pass@gt[,c("FORMAT",samples)]

merged = ins_clean

meta = merged@meta
meta2 = c(
  '##fileformat=VCFv4.1',
  paste0(meta[which(grepl("##source",meta))],";",gsub("##source=","",alu.pass@meta[which(grepl("##source",alu.pass@meta))])),
  paste0(meta[which(grepl("##fileDate",meta))],
         ";",gsub("##fileDate=","",alu.pass@meta[which(grepl("##fileDate",alu.pass@meta))]),
         ";",gsub("##fileDate=","",line1.pass@meta[which(grepl("##fileDate",line1.pass@meta))]),
         ";",gsub("##fileDate=","",sva.pass@meta[which(grepl("##fileDate",sva.pass@meta))]),
         ";",gsub("##fileDate=","",hervk.pass@meta[which(grepl("##fileDate",hervk.pass@meta))]))
  )
meta3 = c(meta2,
  c('##contig=<ID=chr1,length=248956422>',
  '##contig=<ID=chr10,length=133797422>',
  '##contig=<ID=chr11,length=135086622>',
  '##contig=<ID=chr11_KI270721v1_random,length=100316>',
  '##contig=<ID=chr12,length=133275309>',
  '##contig=<ID=chr13,length=114364328>',
  '##contig=<ID=chr14,length=107043718>',
  '##contig=<ID=chr14_GL000009v2_random,length=201709>',
  '##contig=<ID=chr14_GL000194v1_random,length=191469>',
  '##contig=<ID=chr14_GL000225v1_random,length=211173>',
  '##contig=<ID=chr14_KI270722v1_random,length=194050>',
  '##contig=<ID=chr14_KI270723v1_random,length=38115>',
  '##contig=<ID=chr14_KI270724v1_random,length=39555>',
  '##contig=<ID=chr14_KI270725v1_random,length=172810>',
  '##contig=<ID=chr14_KI270726v1_random,length=43739>',
  '##contig=<ID=chr15,length=101991189>',
  '##contig=<ID=chr15_KI270727v1_random,length=448248>',
  '##contig=<ID=chr16,length=90338345>',
  '##contig=<ID=chr16_KI270728v1_random,length=1872759>',
  '##contig=<ID=chr17,length=83257441>',
  '##contig=<ID=chr17_GL000205v2_random,length=185591>',
  '##contig=<ID=chr17_KI270729v1_random,length=280839>',
  '##contig=<ID=chr17_KI270730v1_random,length=112551>',
  '##contig=<ID=chr18,length=80373285>',
  '##contig=<ID=chr19,length=58617616>',
  '##contig=<ID=chr1_KI270706v1_random,length=175055>',
  '##contig=<ID=chr1_KI270707v1_random,length=32032>',
  '##contig=<ID=chr1_KI270708v1_random,length=127682>',
  '##contig=<ID=chr1_KI270709v1_random,length=66860>',
  '##contig=<ID=chr1_KI270710v1_random,length=40176>',
  '##contig=<ID=chr1_KI270711v1_random,length=42210>',
  '##contig=<ID=chr1_KI270712v1_random,length=176043>',
  '##contig=<ID=chr1_KI270713v1_random,length=40745>',
  '##contig=<ID=chr1_KI270714v1_random,length=41717>',
  '##contig=<ID=chr2,length=242193529>',
  '##contig=<ID=chr20,length=64444167>',
  '##contig=<ID=chr21,length=46709983>',
  '##contig=<ID=chr22,length=50818468>',
  '##contig=<ID=chr22_KI270731v1_random,length=150754>',
  '##contig=<ID=chr22_KI270732v1_random,length=41543>',
  '##contig=<ID=chr22_KI270733v1_random,length=179772>',
  '##contig=<ID=chr22_KI270734v1_random,length=165050>',
  '##contig=<ID=chr22_KI270735v1_random,length=42811>',
  '##contig=<ID=chr22_KI270736v1_random,length=181920>',
  '##contig=<ID=chr22_KI270737v1_random,length=103838>',
  '##contig=<ID=chr22_KI270738v1_random,length=99375>',
  '##contig=<ID=chr22_KI270739v1_random,length=73985>',
  '##contig=<ID=chr2_KI270715v1_random,length=161471>',
  '##contig=<ID=chr2_KI270716v1_random,length=153799>',
  '##contig=<ID=chr3,length=198295559>',
  '##contig=<ID=chr3_GL000221v1_random,length=155397>',
  '##contig=<ID=chr4,length=190214555>',
  '##contig=<ID=chr4_GL000008v2_random,length=209709>',
  '##contig=<ID=chr5,length=181538259>',
  '##contig=<ID=chr5_GL000208v1_random,length=92689>',
  '##contig=<ID=chr6,length=170805979>',
  '##contig=<ID=chr7,length=159345973>',
  '##contig=<ID=chr8,length=145138636>',
  '##contig=<ID=chr9,length=138394717>',
  '##contig=<ID=chr9_KI270717v1_random,length=40062>',
  '##contig=<ID=chr9_KI270718v1_random,length=38054>',
  '##contig=<ID=chr9_KI270719v1_random,length=176845>',
  '##contig=<ID=chr9_KI270720v1_random,length=39050>',
  '##contig=<ID=chrEBV,length=171823>',
  '##contig=<ID=chrM,length=16569>',
  '##contig=<ID=chrUn_GL000195v1,length=182896>',
  '##contig=<ID=chrUn_GL000213v1,length=164239>',
  '##contig=<ID=chrUn_GL000214v1,length=137718>',
  '##contig=<ID=chrUn_GL000216v2,length=176608>',
  '##contig=<ID=chrUn_GL000218v1,length=161147>',
  '##contig=<ID=chrUn_GL000219v1,length=179198>',
  '##contig=<ID=chrUn_GL000220v1,length=161802>',
  '##contig=<ID=chrUn_GL000224v1,length=179693>',
  '##contig=<ID=chrUn_GL000226v1,length=15008>',
  '##contig=<ID=chrUn_KI270302v1,length=2274>',
  '##contig=<ID=chrUn_KI270303v1,length=1942>',
  '##contig=<ID=chrUn_KI270304v1,length=2165>',
  '##contig=<ID=chrUn_KI270305v1,length=1472>',
  '##contig=<ID=chrUn_KI270310v1,length=1201>',
  '##contig=<ID=chrUn_KI270311v1,length=12399>',
  '##contig=<ID=chrUn_KI270312v1,length=998>',
  '##contig=<ID=chrUn_KI270315v1,length=2276>',
  '##contig=<ID=chrUn_KI270316v1,length=1444>',
  '##contig=<ID=chrUn_KI270317v1,length=37690>',
  '##contig=<ID=chrUn_KI270320v1,length=4416>',
  '##contig=<ID=chrUn_KI270322v1,length=21476>',
  '##contig=<ID=chrUn_KI270329v1,length=1040>',
  '##contig=<ID=chrUn_KI270330v1,length=1652>',
  '##contig=<ID=chrUn_KI270333v1,length=2699>',
  '##contig=<ID=chrUn_KI270334v1,length=1368>',
  '##contig=<ID=chrUn_KI270335v1,length=1048>',
  '##contig=<ID=chrUn_KI270336v1,length=1026>',
  '##contig=<ID=chrUn_KI270337v1,length=1121>',
  '##contig=<ID=chrUn_KI270338v1,length=1428>',
  '##contig=<ID=chrUn_KI270340v1,length=1428>',
  '##contig=<ID=chrUn_KI270362v1,length=3530>',
  '##contig=<ID=chrUn_KI270363v1,length=1803>',
  '##contig=<ID=chrUn_KI270364v1,length=2855>',
  '##contig=<ID=chrUn_KI270366v1,length=8320>',
  '##contig=<ID=chrUn_KI270371v1,length=2805>',
  '##contig=<ID=chrUn_KI270372v1,length=1650>',
  '##contig=<ID=chrUn_KI270373v1,length=1451>',
  '##contig=<ID=chrUn_KI270374v1,length=2656>',
  '##contig=<ID=chrUn_KI270375v1,length=2378>',
  '##contig=<ID=chrUn_KI270376v1,length=1136>',
  '##contig=<ID=chrUn_KI270378v1,length=1048>',
  '##contig=<ID=chrUn_KI270379v1,length=1045>',
  '##contig=<ID=chrUn_KI270381v1,length=1930>',
  '##contig=<ID=chrUn_KI270382v1,length=4215>',
  '##contig=<ID=chrUn_KI270383v1,length=1750>',
  '##contig=<ID=chrUn_KI270384v1,length=1658>',
  '##contig=<ID=chrUn_KI270385v1,length=990>',
  '##contig=<ID=chrUn_KI270386v1,length=1788>',
  '##contig=<ID=chrUn_KI270387v1,length=1537>',
  '##contig=<ID=chrUn_KI270388v1,length=1216>',
  '##contig=<ID=chrUn_KI270389v1,length=1298>',
  '##contig=<ID=chrUn_KI270390v1,length=2387>',
  '##contig=<ID=chrUn_KI270391v1,length=1484>',
  '##contig=<ID=chrUn_KI270392v1,length=971>',
  '##contig=<ID=chrUn_KI270393v1,length=1308>',
  '##contig=<ID=chrUn_KI270394v1,length=970>',
  '##contig=<ID=chrUn_KI270395v1,length=1143>',
  '##contig=<ID=chrUn_KI270396v1,length=1880>',
  '##contig=<ID=chrUn_KI270411v1,length=2646>',
  '##contig=<ID=chrUn_KI270412v1,length=1179>',
  '##contig=<ID=chrUn_KI270414v1,length=2489>',
  '##contig=<ID=chrUn_KI270417v1,length=2043>',
  '##contig=<ID=chrUn_KI270418v1,length=2145>',
  '##contig=<ID=chrUn_KI270419v1,length=1029>',
  '##contig=<ID=chrUn_KI270420v1,length=2321>',
  '##contig=<ID=chrUn_KI270422v1,length=1445>',
  '##contig=<ID=chrUn_KI270423v1,length=981>',
  '##contig=<ID=chrUn_KI270424v1,length=2140>',
  '##contig=<ID=chrUn_KI270425v1,length=1884>',
  '##contig=<ID=chrUn_KI270429v1,length=1361>',
  '##contig=<ID=chrUn_KI270435v1,length=92983>',
  '##contig=<ID=chrUn_KI270438v1,length=112505>',
  '##contig=<ID=chrUn_KI270442v1,length=392061>',
  '##contig=<ID=chrUn_KI270448v1,length=7992>',
  '##contig=<ID=chrUn_KI270465v1,length=1774>',
  '##contig=<ID=chrUn_KI270466v1,length=1233>',
  '##contig=<ID=chrUn_KI270467v1,length=3920>',
  '##contig=<ID=chrUn_KI270468v1,length=4055>',
  '##contig=<ID=chrUn_KI270507v1,length=5353>',
  '##contig=<ID=chrUn_KI270508v1,length=1951>',
  '##contig=<ID=chrUn_KI270509v1,length=2318>',
  '##contig=<ID=chrUn_KI270510v1,length=2415>',
  '##contig=<ID=chrUn_KI270511v1,length=8127>',
  '##contig=<ID=chrUn_KI270512v1,length=22689>',
  '##contig=<ID=chrUn_KI270515v1,length=6361>',
  '##contig=<ID=chrUn_KI270516v1,length=1300>',
  '##contig=<ID=chrUn_KI270517v1,length=3253>',
  '##contig=<ID=chrUn_KI270518v1,length=2186>',
  '##contig=<ID=chrUn_KI270519v1,length=138126>',
  '##contig=<ID=chrUn_KI270521v1,length=7642>',
  '##contig=<ID=chrUn_KI270522v1,length=5674>',
  '##contig=<ID=chrUn_KI270528v1,length=2983>',
  '##contig=<ID=chrUn_KI270529v1,length=1899>',
  '##contig=<ID=chrUn_KI270530v1,length=2168>',
  '##contig=<ID=chrUn_KI270538v1,length=91309>',
  '##contig=<ID=chrUn_KI270539v1,length=993>',
  '##contig=<ID=chrUn_KI270544v1,length=1202>',
  '##contig=<ID=chrUn_KI270548v1,length=1599>',
  '##contig=<ID=chrUn_KI270579v1,length=31033>',
  '##contig=<ID=chrUn_KI270580v1,length=1553>',
  '##contig=<ID=chrUn_KI270581v1,length=7046>',
  '##contig=<ID=chrUn_KI270582v1,length=6504>',
  '##contig=<ID=chrUn_KI270583v1,length=1400>',
  '##contig=<ID=chrUn_KI270584v1,length=4513>',
  '##contig=<ID=chrUn_KI270587v1,length=2969>',
  '##contig=<ID=chrUn_KI270588v1,length=6158>',
  '##contig=<ID=chrUn_KI270589v1,length=44474>',
  '##contig=<ID=chrUn_KI270590v1,length=4685>',
  '##contig=<ID=chrUn_KI270591v1,length=5796>',
  '##contig=<ID=chrUn_KI270593v1,length=3041>',
  '##contig=<ID=chrUn_KI270741v1,length=157432>',
  '##contig=<ID=chrUn_KI270742v1,length=186739>',
  '##contig=<ID=chrUn_KI270743v1,length=210658>',
  '##contig=<ID=chrUn_KI270744v1,length=168472>',
  '##contig=<ID=chrUn_KI270745v1,length=41891>',
  '##contig=<ID=chrUn_KI270746v1,length=66486>',
  '##contig=<ID=chrUn_KI270747v1,length=198735>',
  '##contig=<ID=chrUn_KI270748v1,length=93321>',
  '##contig=<ID=chrUn_KI270749v1,length=158759>',
  '##contig=<ID=chrUn_KI270750v1,length=148850>',
  '##contig=<ID=chrUn_KI270751v1,length=150742>',
  '##contig=<ID=chrUn_KI270752v1,length=27745>',
  '##contig=<ID=chrUn_KI270753v1,length=62944>',
  '##contig=<ID=chrUn_KI270754v1,length=40191>',
  '##contig=<ID=chrUn_KI270755v1,length=36723>',
  '##contig=<ID=chrUn_KI270756v1,length=79590>',
  '##contig=<ID=chrUn_KI270757v1,length=71251>',
  '##contig=<ID=chrX,length=156040895>',
  '##contig=<ID=chrY,length=57227415>',
  '##contig=<ID=chrY_KI270740v1_random,length=37240>',
  '##ALT=<ID=DEL,Description="Deletion">',
  '##ALT=<ID=DUP,Description="Duplication">',
  '##ALT=<ID=INV,Description="Inversion">',
  '##ALT=<ID=BND,Description="Translocation">',
  '##ALT=<ID=INS,Description="Insertion">',
  '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
  '##ALT=<ID=INS:ME:LINE1,Description="Insertion of LINE1 element">',
  '##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">',
  '##INFO=<ID=CIEND,Number=2,Type=String,Description="PE confidence interval around END">',
  '##INFO=<ID=CIPOS,Number=2,Type=String,Description="PE confidence interval around POS">',
  '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">',
  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
  '##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">',
  '##INFO=<ID=RE,Number=1,Type=Integer,Description="read support">',
  '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
  '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">',
  '##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">',
  '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Method for generating this merged VCF file.">',
  '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">',
  '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of supporting samples.">',
  '##INFO=<ID=SUPP,Number=1,Type=String,Description="Number of samples supporting the variant">',
  '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Indicating the direction of the reads with respect to the type and breakpoint.">',
  '##INFO=<ID=ASSESS,Number=1,Type=Integer,Description="Provides information on evidence availible to decide insertion site.0 = No overlapping reads at site;1 = Imprecise breakpoint due to greater than expected distance between evidence;2 = discordant pair evidence only -- No split read information;3 = left side TSD evidence only;4 = right side TSD evidence only;5 = TSD decided with split reads, highest possible quality">',
  '##INFO=<ID=TSD,Number=1,Type=String,Description="Precise Target Site Duplication for bases, if unknown, value will be NULL">',
  '##INFO=<ID=INTERNAL,Number=2,Type=String,Description="If insertion internal or close to a gene, listed here followed by a discriptor of the location in the gene (either INTRON, EXON_#, 5_UTR, 3_UTR, PROMOTER, or TERMINATOR). If multiple genes intersected, will be seperated by pipe">',
  '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY; If START or END is unknown, will be -1; If POLARITY is unknown, will be null">',
  '##INFO=<ID=DIFF,Number=.,Type=String,Description="Coverage and Differences in relation to the MEI reference. Form is %2XCoverage:Differences, with differences delimited by ,">',
  '##INFO=<ID=LP,Number=1,Type=Integer,Description="Total number of discordant pairs supporting the left side of the breakpont">',
  '##INFO=<ID=RP,Number=1,Type=Integer,Description="Total number of discordant pairs supporting the right side of the breakpont">',
  '##INFO=<ID=RA,Number=1,Type=Float,Description="Ratio between LP and RP, reported as log2(LP / RP)">',
  '##INFO=<ID=PRIOR,Number=1,Type=String,Description="True if this site was not discovered in this dataset, but was included on a provided priors list">',
  '##INFO=<ID=SR,Number=1,Type=Integer,Description="Total number of SRs at the estimated breakpoint for this site. Recomended to filter sites with <= 2 SRs">',
  '##FILTER=<ID=s25,Description="Greater than 25.0% of samples do not have data">',
  '##FILTER=<ID=rSD,Description="Ratio of LP to RP is greater than 2.0 standard deviations">',
  '##FILTER=<ID=hDP,Description="More than the expected number of discordant pairs at this site are also split">',
  '##FILTER=<ID=ac0,Description="No individuals in this VCF file were identified with this insertion">',
  '##FILTER=<ID=lc,Description="MEI is embeded in a low complexity region">',
  '##FILTER=<ID=PASS,Description="Passed all filters">',
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  '##FORMAT=<ID=PSV,Number=1,Type=String,Description="Previous support vector">',
  '##FORMAT=<ID=LN,Number=1,Type=Integer,Description="predicted length">',
  '##FORMAT=<ID=DR,Number=2,Type=Integer,Description="# supporting reference,variant reads in that order">',
  '##FORMAT=<ID=ST,Number=1,Type=String,Description="Strand of SVs">',
  '##FORMAT=<ID=QV,Number=1,Type=String,Description="Quality values: if not defined a . otherwise the reported value.">',
  '##FORMAT=<ID=TY,Number=1,Type=String,Description="Types">',
  '##FORMAT=<ID=ID,Number=1,Type=String,Description="Variant ID from input.">',
  '##FORMAT=<ID=RAL,Number=1,Type=String,Description="Reference allele sequence reported from input.">',
  '##FORMAT=<ID=AAL,Number=1,Type=String,Description="Alternative allele sequence reported from input.">',
  '##FORMAT=<ID=CO,Number=1,Type=String,Description="Coordinates">')
)
merged@meta = meta3
merged@gt = rbind(gsub("\\./\\.","0/0",ins_clean@gt), alu.pass@gt, line1.pass@gt, sva.pass@gt, hervk.pass@gt)
merged@fix = rbind(merged@fix, alu.pass@fix, line1.pass@fix, sva.pass@fix, hervk.pass@fix)
merged@fix[,8] = paste0("SVTYPE=",extract.info(merged,"SVTYPE"),";","SVLEN=",extract.info(merged,"SVLEN"))
merged@gt = gsub("(.*?):(.*)","\\1",merged@gt)

write.vcf(merged, file = output_vcf)