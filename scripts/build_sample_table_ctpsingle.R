library(stringr)
library(dplyr)
library(tidyr)

setwd("/home/uvi/be/avs/store/dNdS_clones")

#Functional data
##Load purity and ploidy data
samples.info<-read.table("metadata/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                         sep="\t",header=T, stringsAsFactors = F)
names(samples.info)[2]<-"sample.code"

#Edit sample column to retrieve sample name like we call them
samples.info$sample<-str_extract(samples.info$sample.code,"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")

samples.info<-samples.info[,c(11,4:9)]

##Load clinical data
samples.clinical<-read.table("metadata/TCGA-CDR-Clinical_data.tsv",sep="\t", header=F, stringsAsFactors = T) %>% select(-c("X","race","initial_pathologic_dx_year",
                            "histological_grade","cause_of_death","menopause_status",
                            "margin_status","residual_tumor","Redaction")) %>% rename(sample = bcr_patient_barcode)
                

#Declare vector with sample names
samples<-samples.info$sample

#Add number of initial mutations per sample
snv.files<-list.files("data/snv/")
snv.samples<-gsub("_SNV.*","",snv.files)

mutations<-integer()

for (i in 1:length(snv.samples)) {
    sample<-snv.samples[i]
    #snv<-read.table(gzfile(paste("data/snv/",sample,"_SNV.tsv.gz",sep="")), header=T)
    mutations[i]<- read.table(gzfile(paste("data/snv/",sample,"_SNV.tsv.gz",sep="")), header=T) %>% nrow()
}

#Declare if the samples are hypermutators or not
quant.mut<-quantile(mutations, na.rm=TRUE)
q3<-quant.mut[4]
iqr<-(quant.mut["75%"] - quant.mut["25%"])
hm.treshold<- (q3+ 1.5*iqr)

is.hypermutator<-mutations > hm.treshold

samples.muts<-data.frame(sample=snv.samples,total_mutations=mutations, hypermutator = is.hypermutator, stringsAsFactors = F)

#Number of clusters inferred by CTPsingle
cluster.files<-list.files(path = "data/ctpsingle/ctpsingle_output",
                       pattern = "cluster_assignments.txt", 
                       full.names = T, recursive=TRUE)

cluster_samples<-str_extract(cluster.files,"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")

cluster_mutations<-integer()
cluster_number<-integer()

for (i in 1:length(cluster.files)) {
  
  cluster.file<-read.table(cluster.files[i],sep="\t",header=T,stringsAsFactors = F)
  cluster_mutations[i]<-nrow(subset(cluster.file))
  cluster_number[i]<-length(table(cluster.file$mostLikely))
}

cluster.df<-data.frame(sample=cluster_samples, clusters=cluster_number, clustered_mutations =  cluster_mutations, stringsAsFactors = F)

#Load topologies data
topo.data <-read.table("tables/TCGA_topologies_data_ctpsingle.txt",sep="\t",header=T, stringsAsFactors=F)
topo.data <- topo.data %>% separate(sample,c("sample","tree"),sep="_tree_") %>% group_by(sample) %>% summarize(seq.length= max(seqlength),global.dNdS = mean(global.dNdS),LRT_M0vsFB= mean(LRT_M0vsFB)) 

#Merge and save file
samples.data<-Reduce(function(x,y) merge(x,y, all=TR, by="sample"), list(samples.info, samples.clinical,samples.muts, cluster.df,topo.data))

write.table(samples.data,"tables/TCGA_samples_data_ctpsingle.txt",row.names = F, col.names = T,
            sep="\t", quote = FALSE)

#This code is for parsing samples with a single clonal tree topology
seq.files=list.files(path = "data/ctpsingle/seqs",pattern="_clone_seqs.fas", full.names=T, recursive = T)

seq_sample<-str_extract(seq.files, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
seq_number<-integer()

for (i in 1:length(seq.files)) {
  seq_file<-readLines(seq.files[i])
  seq_number[i]<-length(grep(">clone",seq_file))
}

seq.df<-data.frame(sample=seq_sample, sequences=seq_number, stringsAsFactors = F)

#Add alignment length
paml.folders<-list.dirs(path="data/ctpsingle/paml", recursive = F)

paml.out.files<-list.files("data/ctpsingle/paml", pattern = "out", 
                              full.names = T, recursive=T)

paml.M0.out.files<-paml.out.files[grepl("M0~.*",paml.out.files)]

paml.samples<-str_extract(paml.M0.out.files, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")

nseqs<-numeric()
seqlength<-numeric()

for (i in 1:length(paml.M0.out.files)) {
  m0.file<-readLines(paml.M0.out.files[i])
  firstLine<-m0.file[1]
  nseqs[i]<-as.integer(sapply(strsplit(sub("^\\s+","", firstLine),split="\\s+"),"[",1))
  seqlength[i]<-as.integer(sapply(strsplit(sub("^\\s+","", firstLine),split="\\s+"),"[",2))
  
}

seqlength.df<-data.frame(sample=paml.samples,nseqs, seqlength, stringsAsFactors = F)

#Add global dN/dS and LRT of free-branch vs M0 models of Codeml
paml.fb.files<-list.files(path="data/ctpsingle/paml",pattern="M0_vs_fb.txt",full.names = T, recursive = T)

paml.samples<-str_extract(paml.files, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
dNdS<-numeric()
lrt.fb<-numeric()

for (i in 1:length(paml.fb.files)) {
  
  paml.file<-readLines(paml.fb.files[i])
  
  #dN/dS
  lastLine<-paml.file[length(paml.file)]
  dNdS[i]<-gsub(".*: ","", lastLine)
  
  #LRT
  lrtLine<-paml.file[grep("p-value",paml.file)+2]
  lrtLineSplit<-strsplit(lrtLine,"\\| ")
  lrt.fb[i]<-lrtLineSplit[[1]][length(lrtLineSplit[[1]])]
}

paml.df<-data.frame(sample=paml.samples, global.dNdS =dNdS, LRT_M0vsFB=lrt.fb, signif.fb=lrt.fb < 0.05,
                    stringsAsFactors = F)
paml.df[,2:3]<-sapply(paml.df[,c(2,3)],function(x) as.numeric(x))

#Add global (bulk) dN/dS from from dndscv
dndscv.files=list.files(path = "data/dndscv/",pattern="bulk_dnds.txt",
                       full.names=T, recursive = T)

dndscv.samples<- str_extract(dndscv.files, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
dndscv_mis <-numeric()
dndscv_all <-numeric()

for (i in 1:length(dndscv.files)) {
dndscv.file<-read.table(dndscv.files[i],sep="\t",header=T)
dndscv_mis[i]<-dndscv.file[1,"mle"]
dndscv_all[i]<-dndscv.file[2,"mle"]
}

dndscv.df<-as.data.frame(cbind(dndscv.samples,dndscv_mis, dndscv_all))
names(dndscv.df)[2:3]<-c("bulk_dndscv_mis","bulk_dndscv_all")

##Molecular clock LRT
clock.dirs<-list.dirs("data/pyclone/molecular_clock",recursive=F)
samples.clocks<-str_extract(clock.dirs, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")

lnl.relaxed.clock<-numeric()
lnl.global.clock<-numeric()

for (i in 1:length(samples.clocks)) {
  
  #No-clock model
  noclock.file<-readLines(paste(clock.dirs[i],"/no_clock_pruned/",samples.clocks[i],"_noclock.out",sep = ""))
  noclock.lnl.line<-grep("lnL",noclock.file, value = T)
  lnl.relaxed.clock[i]<-as.numeric(sapply(strsplit(noclock.lnl.line,": +| +"),"[",5))
  
  #Global clock model
  globalclock.file<-readLines(paste(clock.dirs[i],"/clock_pruned/",samples.clocks[i],"_clock.out",sep = ""))
  globalclock.lnl.line<-grep("lnL",globalclock.file, value = T)
  lnl.global.clock[i]<-as.numeric(sapply(strsplit(globalclock.lnl.line,": +| +"),"[",5))
}

lrt.global.clock<-2*(lnl.relaxed.clock - lnl.global.clock)

clock.data<-data.frame(sample=samples.clocks,LRT_global_vs_relaxed_clock=lrt.global.clock)

np.clock<-(2*samples.data$nseqs)-1
np.noclock<-(3*samples.data$nseqs)-2
p.global.clock<-1-pchisq(as.numeric(samples.data$LRT_global_vs_relaxed_clock),df = (np.noclock -np.clock))
rejected.global.clock<-p.global.clock< 0.05

#Merge and save file
samples.data<-Reduce(function(x,y) merge(x,y, all=TRUE, by="sample"), list(samples.info, samples.clinical,samples.muts, cluster.df,seq.df, paml.df, clock.data))

write.table(samples.data,"tables/TCGA_samples_data_ctpsingle.txt",row.names = F, col.names = T,
            sep="\t", quote = FALSE)
