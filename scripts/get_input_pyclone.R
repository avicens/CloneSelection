#Set working directory
setwd("~/Dropbox/proyectos/dNdS_clones/")

#Load libraries
library("bedr")
library("dplyr")

#Load data
##Single nucleotide variations
snv.header<-read.table("TCGA-02-003/data/headers.maf", sep="\t", header = T)
snv<-read.table("data/maf/TCGA-OR-A5JI-01.maf",sep="\t",header=T, stringsAsFactors = F)
snv<-snv[,c(1,5:7,40:45)]

##Copy number variations
cnv.header<-read.table("TCGA-02-003/headers.seg.txt", sep="\t", header = T)
cnv<-read.table("data/cnv/TCGA-OR-A5JI-01.cnv",sep="\t", header = F, stringsAsFactors = F)

##SNV file
snv.bedr<-snv[,c(2:4,1,5:10)]
is.valid.region(snv.bedr, check.chr = FALSE, check.zero.based = FALSE)
snv.bedr.sorted<-bedr.sort.region(snv.bedr, check.zero.based = FALSE, check.chr = F)
snv.bedr.sorted.merged<-bedr.merge.region(snv.bedr.sorted, check.zero.based = FALSE, check.chr = FALSE)
snv.bedr.filtered<-snv.bedr.sorted[!duplicated(snv.bedr.sorted[c("Chromosome","Start_Position","End_Position")]),]

##CNV file
cnv.bedr<-cnv[,c(2:4,7,8)]
cnv.bedr$Chromosome<-as.character(cnv.bedr$Chromosome)
is.valid.region(cnv.bedr, check.chr = FALSE, check.zero.based = FALSE)
cnv.bedr.sorted<-bedr.sort.region(cnv.bedr, check.zero.based = FALSE, check.chr = F)

##Intersect
snv.int.cnv<-bedr.join.region(snv.bedr.filtered, cnv.bedr.sorted, check.chr = F, check.zero.based = F)

##Discard redundant chromosome coordinates columnd
snv.int.cnv<-snv.int.cnv[,-c(11:13)]

##Discard sex chromosomes and aberrant copy numbers
snv.int.cnv.filt<-snv.int.cnv %>% filter(Chromosome != "X" | Chromosome != "Y") %>% filter(Modal_HSCN_1 != "." & Modal_HSCN_2 != "-1")

##Get pyClone input
pyClone.input<-data.frame(mutation_id = paste(paste(snv.int.cnv.filt[,1], ":",snv.int.cnv.filt[,2],"-", snv.int.cnv.filt[,3],"_",snv.int.cnv.filt[,4], sep="")), 
                          ref_counts=snv.int.cnv.filt$t_ref_count, 
                          var_counts = snv.int.cnv.filt$t_alt_count, normal_cn = rep(2,nrow(snv.int.cnv.filt)), 
                          minor_cn = snv.int.cnv.filt$Modal_HSCN_1, major_cn = snv.int.cnv.filt$Modal_HSCN_2,
                          stringsAsFactors = F)

pyClone.input[,2:6]<-apply(pyClone.input[,2:6],2, function(x) as.integer(x))

#Save data frames
write.table(snv.int.cnv, "TCGA-02-003/TCGA-02-0003.tsv", row.names = F, col.names = T, sep="\t", quote = F)
write.table(pyClone.input, "TCGA-02-003/TCGA-02-0003.pyclone", row.names = F, col.names = T, sep="\t", quote = F)
