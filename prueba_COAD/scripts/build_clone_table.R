library(reshape2)

args<-commandArgs(TRUE)

if (length(args) != 1) {
  cat("Usage: Rscript build_mutation_table.R  <sampleName>\n")
}

sample<-args[1]

workDir<-"/mnt/netapp2/Store_uni/home/uvi/be/avs/dNdS_clones/prueba_COAD"
setwd(workDir)
outFile<-paste("data/clones/",sample,"_clones_table.txt",sep="")

#Load files
cluster<-read.table(paste("trace/pyclone_output/",sample,"/tables/cluster.tsv",sep=""),
                    header=T, stringsAsFactors = F, sep="\t")

cluster<-cluster[cluster$size >=3,] #Filter out clones with less than 3 mutations

mutations<-read.table(paste("data/mutations/",sample,"_mutations_table.txt",sep=""),header=T, stringsAsFactors = T)

loci<-read.table(paste("trace/loci/",sample,"_loci.txt",sep=""),
                 header=T, stringsAsFactors = F)

tree<-read.table(paste("trace/clone_trees/",sample,"_clone_tree.txt",sep=""),
                 header=T, stringsAsFactors = F)

#Estimate normalized cluster size
nmut<-nrow(loci)
cluster$norm_size<-round(sapply(cluster$size, function(x) (x/nmut)*100),2)

#Get mutations
mf<-with(mutations, table(cluster_id, Variant_Classification))
cluster.mf<-merge(cluster, mf)
cluster.cast<- dcast(cluster.mf, cluster_id ~ Variant_Classification, value.var = "Freq")
cluster2<-merge(cluster, cluster.cast, by = "cluster_id")

#Get normalized clone IDs and ancestral state
nc<-as.data.frame(table(loci$cluster)) #Get data for clone size

cluster3<-merge(merge(cluster2,nc,by.x="size",by.y="Freq"), tree, by.x="Var1",by.y="lab")
cluster3$parental_clone<-cluster3$parent ==-1 #Add column with parental state

names(cluster3)[1]<-"clone_id" #Edit name of first column (Clone_id)

if (file.exists(paste("trace/dNdS/",sample,"_dnds_per_clone.tsv",sep=""))) {
  
  dnds<-read.table(paste("trace/dNdS/",sample,"_dnds_per_clone.tsv",sep=""),
                   header=T, stringsAsFactors = F)


  #Merge data frames (2)
  cluster4<-merge(cluster3, dnds, by.x="clone_id", by.y ="Clone")

  #Save file
  write.table(cluster4,outFile, col.names = T, row.names = F, quote =F, sep= "\t")
} else  {
  
  cluster3[,c("fg_dNdS", "bg_dNdS", "LRT_bfree_bneut")]<-"NA"
  write.table(cluster3,outFile, col.names = T, row.names = F, quote =F, sep= "\t")
  }