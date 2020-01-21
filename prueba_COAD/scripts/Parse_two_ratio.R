library(stringr)
#Input arguments
args<-commandArgs(TRUE)

if(length(args) != 2) {
    cat("Usage: Rscript Parse_two_ratio.R <paml_dir> <output file>\nPlease provide the two required arguments\n")
    quit()
}
paml.dir<-args[1]
outfile<-args[2]

#Functions
##Function obtain_dnds_fgbg
obtain_dnds_fgbg <- function(sample) {
    paml.file.path<-paste(paml.dir,"/",sample,"/two-ratio/",sample,"_two-ratio.out", sep="")
    paml.file<-readLines(paml.file.path)
    dnds.line<-paml.file[str_detect(paml.file,"\\(dN/dS\\)")]
    dnds.values<-sapply(strsplit(dnds.line,split=":  "),"[",2)

    dnds.bg<-as.numeric(sapply(strsplit(dnds.values," "),"[",1))
    dnds.fg<-as.numeric(sapply(strsplit(dnds.values," "),"[",2))
    
    dnds.row<-c(sample,as.numeric(dnds.fg),as.numeric(dnds.bg))

    return(dnds.row)
}

samples <- list.dirs(paml.dir,recursive=FALSE, full.names=FALSE)
dnds_fgbg_df<-data.frame(sample=character(), fg.dnds=numeric(), bg.dnds=numeric(),stringsAsFactors=FALSE)

for (i in 1:length(samples)) {
    smp<-samples[i]
    dnds<-obtain_dnds_fgbg(smp)
    dnds_fgbg_df[i,]<-dnds
}

cat("Exporting cdata to file\n")
write.table(dnds_fgbg_df, outfile, sep="\t", quote=F, row.names=F, col.names=T)

cat("Work finished\n")

