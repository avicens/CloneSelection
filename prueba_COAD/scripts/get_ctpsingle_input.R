library(stringr)
library(dplyr)

args<-commandArgs(TRUE)

intersFile<-args[1]

sample<-str_extract(intersFile,"TCGA-[0-9A-Z]{2}-[A-Z0-9]{4}")

snv.int.cnv<-read.table(intersFile, stringsAsFactors = F, header=F)
names(snv.int.cnv)<-c("Chr","Start","End",	"Reference",	"Alternate",	"Variant_freq",	"Ref_freq",
                 "Hugo_Symbol"	,"Variant_Classification","Chromosome",	"Start_position",	"End_position",
                 "nProbes"	,"Modal_HSCN_1",	"Modal_HSCN_2",
                 "Modal_Total_CN"	,"purity"	,"ploidy")


ctps.input<-snv.int.cnv[,c(1:2,5,4,6,7,16)]
names(ctps.input)<-c("Chromosome","Position","Mutant","Reference","Mcount","Rcount","Multiplier")
ctps.input$Gender<-rep("Unknown",nrow(ctps.input))

#Filter file removeing variants ins egments with CN != 2
ctps.input2 <-ctps.input %>% filter(Multiplier == 2)


write.table(ctps.input2,paste("data/ctpsingle/ctpsingle_input2/",sample,"_ctpsingle_input.txt",sep=""), sep="\t", row.names=F, col.names=T, quote=F)