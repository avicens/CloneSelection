library(stringr)
library(plyr)

setwd("/home/uvi/be/avs/store/dNdS_clones/prueba_COAD")

#Load table with information of samples
samples.info<-read.table("metadata/COAD-TCGA_samples_info.txt",
                         sep="\t",header=T, stringsAsFactors = F)
names(samples.info)[2]<-"sample.code"

#Edit sample column to retrieve sample name like we call them
samples.info$sample<-gsub("-[A-Z,0-9]{3}-[A-Z,0-9]{3}-[A-Z,0-9]{4}-[A-Z,0-9]{2}$","",
                          samples.info$sample.code)

#Reorder and filter columns
samples.info<-samples.info[,c(1,2,11,4:9)]

#Declare vector with sample names
samples<-samples.info$sample

#Add clinical data
sample.clinical<-read.table("metadata/COAD-CDR-Clinical_data_compacted.tsv",sep="\t", header=T,
                            stringsAsFactors = F)

#Merge sample info and clinical data
samples.data<-merge(samples.info,sample.clinical[,c("bcr_patient_barcode","type","age_at_initial_pathologic_diagnosis", "ajcc_pathologic_tumor_stage",
                                                  "histological_type","birth_days_to","tumor_status","last_contact_days_to",
                                                  "treatment_outcome_first_course","OS.time","DFI.time","PFI.time")],
                    by.x="sample", by.y="bcr_patient_barcode",all.x=T)

samples.data[,c("OS.time","DFI.time","PFI.time")]<-sapply(samples.data[,c("OS.time","DFI.time","PFI.time")],function(x) as.numeric(x))

#Add number of mutations per sample
mutations<-integer()

for (i in 1:length(samples)) {
snv<-read.table(gzfile(paste("metadata/snv/",samples[i],"_SNV.tsv.gz",sep="")),
                header=T)
mutations[i]<-nrow(snv)
}

samples.muts<-data.frame(sample=samples,total_mutations=mutations, stringsAsFactors = F)
samples.data<-merge(samples.data, samples.muts, by="sample", all.x=TRUE)


#Number of clusters inferred by PyClone (Clusters with size > 1)
cluster.files<-list.files(path = "trace/pyclone_output",
                       pattern = "cluster.tsv", 
                       full.names = T, recursive=TRUE)

cluster_sample<-character()
cluster_number<-integer()

for (i in 1:length(cluster.files)) {
  cluster.file<-read.table(cluster.files[i],sep="\t",header=T,stringsAsFactors = F)
  cluster_sample[i]<-cluster.file$sample_id[1]
  cluster_number[i]<-nrow(subset(cluster.file, size >1))
}

cluster.df<-data.frame(sample=cluster_sample, clusters=cluster_number, stringsAsFactors = F)
samples.data<-merge(samples.data,cluster.df,by="sample",all.x=TRUE)

#Add number of clustered mutations (loci)

loci.files<-list.files("trace/loci/", pattern = "TCGA-.*", full.names = T)

loci_samples<-character()
cluster_mutations<-integer()

for (i in 1:length(loci.files)) {
  loci.file<-read.table(loci.files[i],sep="\t",header=T,stringsAsFactors = F)
  loci_samples[i]<-sapply(strsplit(loci.file$sample_id,"_"),"[",1)[1]
  cluster_mutations[i]<-nrow(loci.file)
}

loci.df<-data.frame(sample=loci_samples, clustered_mutations= cluster_mutations,
                    stringsAsFactors = F)
samples.data<-merge(samples.data, loci.df,by="sample")


#Number of clones deconvoluted (ClonEvol)
clone.tree.files=list.files("trace/clone_trees",pattern="*_tree.txt", full.names = T)
clone.tree.samples<-str_extract(clone.tree.files, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")

clones<-integer()
for (i in 1:length(clone.tree.files)) {
tree.file<-read.table(clone.tree.files[i], header=T)
clones[i]<-nrow(tree.file)
}

clones.df<-data.frame(sample=clone.tree.samples, clones=clones, stringsAsFactors = F)
samples.data<-merge(samples.data, clones.df,by="sample",all.x=T)


#Add number of clonal sequences
seq.files=list.files(path = "trace/seqs",pattern="_clone_ref_seqs.fas", full.names=T, recursive = T)

seq_sample<-str_extract(seq.files, "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
seq_number<-integer()

for (i in 1:length(seq.files)) {
  seq_file<-readLines(seq.files[i])
  seq_number[i]<-length(grep(">clone",seq_file))
}

seq.df<-data.frame(sample=seq_sample, sequences=seq_number, stringsAsFactors = F)
samples.d
#Add global (bulk) dN/dS from from dndscv
#dndscv.files=list.files(path = "data/dndscv/",pattern="bulk_dnds.txt",
 #                       full.names=T, recursive = T)

#dndscv.samples<- sapply(strsplit(dndscv.files,"//|_"),"[",2)
#dndscv_mis <-numeric()
#dndscv_all <-numeric()

#for (i in 1:length(dndscv.files)) {
#dndscv.file<-read.table(dndscv.files[i],sep="\t",header=T)
#dndscv_mis[i]<-dndscv.file[1,"mle"]
#dndscv_all[i]<-dndscv.file[2,"mle"]
#}

#dndscv.df<-as.data.frame(cbind(dndscv.samples,dndscv_mis, dndscv_all))
#names(dndscv.df)[2:3]<-c("bulk_dndscv_mis","bulk_dndscv_all")
#samples.data<-merge(samples.data, dndscv.df,by.x="sample",by.y="dndscv.samples",all.x=T)

#Add number of sequences and alignment length
paml.folders<-list.dirs(path="data/paml_2", recursive = F)

paml.M0.out.files<-list.files("data/paml_2", pattern = "out", 
                              full.names = T, recursive=T)

paml.M0.out.files<-paml.M0.out.files[grepl("M0",paml.M0.out.files)]

paml.samples<-sapply(strsplit(paml.M0.out.files,"/"),"[",3)
nseqs<-numeric()
seqlength<-numeric()

for (i in 1:length(paml.M0.out.files)) {
  m0.file<-readLines(paml.M0.out.files[i])
  firstLine<-m0.file[1]
  nseqs[i]<-as.integer(sapply(strsplit(sub("^\\s+","", firstLine),split="\\s+"),"[",1))
  seqlength[i]<-as.integer(sapply(strsplit(sub("^\\s+","", firstLine),split="\\s+"),"[",2))
  
}

seqlength.df<-data.frame(sample=paml.samples,nseqs, seqlength, stringsAsFactors = F)
samples.data<-merge(samples.data,seqlength.df, by ="sample", all.x=T)

#Add global dN/dS and LRT of free-branch vs M0 models of Codeml
paml.files<-list.files(path="data/paml_2",pattern="M0_vs_fb.txt",full.names = T, recursive = T)

paml.samples<-sapply(strsplit(paml.files,"/"),"[",3)
dNdS<-numeric()
lrt.fb<-numeric()

for (i in 1:length(paml.files)) {
  
  paml.file<-readLines(paml.files[i])
  
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

samples.data<-merge(samples.data,paml.df,by="sample",all.x=T)

##Molecular clock LRT
clock.dirs<-list.dirs("trace/molecular_clock",recursive=F)
samples.clocks<-sapply(strsplit(clock.dirs,"/"),"[",3)

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

samples.data<-merge(samples.data,clock.data,by="sample", all.x=TRUE)

np.clock<-(2*samples.data$nseqs)-1
np.noclock<-(3*samples.data$nseqs)-2
p.global.clock<-1-pchisq(as.numeric(samples.data$LRT_global_vs_relaxed_clock.y),df = (np.noclock -np.clock))
rejected.global.clock<-p.global.clock< 0.05

samples.data<-cbind(samples.data, rejected.global.clock)


#Save file
write.table(samples.data,"COAD_samples_data.txt",row.names = F, col.names = T,
            sep="\t", quote = FALSE)
