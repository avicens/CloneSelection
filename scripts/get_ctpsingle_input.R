suppressMessages(library(stringr, quietly=TRUE))
suppressMessages(library(dplyr, quietly = TRUE))

args<-commandArgs(TRUE)

if (length(args) != 2) {
    cat("Input arguments must be inserted
    Usage: get_ctpsingle_input.R <intersect file> <output file name>\n")
    quit()
}
intersFile<-args[1]
outFile<-args[2]




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


write.table(ctps.input2,outFile, sep="\t", row.names=F, col.names=T, quote=F)