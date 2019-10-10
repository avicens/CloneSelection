#Load libraries
library(dndscv)
library(dplyr)

#Input arguments
args<-commandArgs(TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript get_clone_sequences.R  <in SNV file> <in Loci file> <dndscv input> <dNdScvout>\n")
  quit()
}

snv.file=args[1]
sample=args[2]
dndscv.infile=args[3]
dndscv.outfile=args[4]

#1)Variants file
snv<-read.table(gzfile(snv.file),header=T)

##2)Table with loci information from PyClone
#pyclone.loci<-read.table(in.loci.file, sep="\t", header=T, stringsAsFactors = F)

#Merge SNV and loci files (only matching mutations)
#snv$mutation_id<-with(snv,paste(Hugo_Symbol,":",Chr,"-",Start_Position,"_",End_Position,sep=""))
#muts<-merge(snv,pyclone.loci,by="mutation_id")

#Estimate global dN/dS for the bulk
dndscv.input <- snv %>% select("Chr","Start_Position","Reference","Alternate") %>% mutate(Sample=sample) %>% select(Sample, everything()) 
dndscv.output<-dndscv(dndscv.input, outp = 1)
global.dnds<-dndscv.output$globaldnds[c(1,5),]
#bulk.dnds$sample<-rep("Bulk",nrow(bulk.dnds))

#Save files
cat("Saving files\n")
write.table(dndscv.input, dndscv.infile, row.names=F, col.names=T, quote=F, sep="\t")
write.table(global.dnds, dndscv.outfile, row.names=F, col.names = TRUE, quote=F,sep="\t")

cat ("work finished\n")