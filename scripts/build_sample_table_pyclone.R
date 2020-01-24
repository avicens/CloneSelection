library(stringr)
setwd("/home/uvi/be/avs/store/dNdS_clones/")

samples.info<-read.table(paste("samples/COAD/COAD_samples_info_278smp.txt",sep="")
                ,header=T, stringsAsFactors = F)

#Add clinical data
sample.clinical<-read.table("metadata/COAD-CDR-Clinical_data_compacted.tsv",sep="\t",
                     header=T)

#Merge sample info and clinical data
samples.data<-merge(samples.info,sample.clinical,by.x="sample", by.y="bcr_patient_barcode",all.x=T)

#Add number of mutations per sample
samples<-samples.data$sample
mutations<-integer()

for (i in 1:length(samples)) {
snv<-read.table(gzfile(paste("data/snv/",samples[i],"_SNV.tsv.gz",sep="")),
                header=T)
mutations[i]<-nrow(snv)
}

samples.data<-cbind(samples.data,mutations)

#Add number of mutations used for cluster infering
loci.files<-list.files(path = "data/pyclone_output",
                       pattern = "loci_processed.tsv", 
                       full.names = T, recursive=TRUE)

loci_samples<-character()
cluster_mutations<-integer()
for (i in 1:length(loci.files)) {
  loci.file<-read.table(loci.files[i],sep="\t",header=T,stringsAsFactors = F)
  loci_samples[i]<-sapply(strsplit(loci.file$sample_id,"_"),"[",1)[1]
  cluster_mutations[i]<-nrow(loci.file)
}

loci.df<-as.data.frame(cbind(loci_samples, cluster_mutations))

samples.data<-merge(samples.data,loci.df,by.x="sample",by.y="loci_samples",all.x=TRUE)

#Add number of clusters per samples
clone.tree.files=list.files("data/clone_trees",pattern="*_tree.txt", full.names = T)
clone.tree.samples<-sapply(strsplit(clone.tree.files,"/|_"),"[",4)

clusters<-integer()
for (i in 1:length(clone.tree.files)) {
tree.file<-read.table(clone.tree.files[i], header=T)
clusters[i]<-nrow(tree.file)
}

sample.clusters<-as.data.frame(cbind(clone.tree.samples, clusters))
samples.data<-merge(samples.data,sample.clusters,by.x="sample",by.y="clone.tree.samples",all.x=T)
samples.data$clusters<-as.integer(samples.data$clusters)

#Add global omega and LRT M0 vs Free-branch model
paml.dirs<-list.dirs(path="data/paml",recursive = F)
paml.sample<-sapply(strsplit(paml.dirs,"/"),"[",3)

dNdS<-numeric()
lrt<-numeric()

for (i in 1:length(paml.dirs)) {
  paml.file<-readLines(paste(paml.dirs[i],"/",paml.sample[i],"_M0_vs_fb.out",sep=""))
  lastLine<-paml.file[length(paml.file)]
  dNdS[i]<-gsub(".*: ","", lastLine)
  lrtLine<-paml.file[12]
  lrtLineSplit<-strsplit(lrtLine,"\\| ")
  lrt[i]<-lrtLineSplit[[1]][length(lrtLineSplit[[1]])]
}

paml.df<-data.frame(cbind(paml.sample,dNdS,lrt),stringsAsFactors = F)
paml.df[,2:3]<-sapply(paml.df[,2:3],function(x) as.numeric(x))

samples.data<-merge(samples.data,paml.df,by.x="sample", by.y="paml.sample",all.x=T)


write.table(samples.data,"metadata/COAD_samples_data.txt",row.names = F, col.names = T,
            sep="\t", quote = FALSE)
