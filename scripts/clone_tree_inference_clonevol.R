#Load libraries
library(clonevol)
library(dplyr)

#Input arguments
if (length(args) != 1) {
  cat("Usage: Rscript clone_tree_inference_clonevol.R  <sample name (e.g: TCGA-2F-A9KT)>\n")
}

args<-commandArgs(TRUE)
sample<-args[1]
loci.output<-paste("data/pyclone/pyclone_output/",sample,"/",sample,"_loci.tsv",sep="")
tree.output<-paste("data/pyclone/clonevol/",sample,"_clone_tree.txt",sep="")

#Laod function for converting clone trees to branch
source("scripts/convert.clone.to.branch.R")

#Infer clone evolution from the loci table
pyclone.loci<-read.table(paste("data/pyclone/pyclone_output/",sample,"/tables/loci.tsv",sep=""), sep="\t", header=T, stringsAsFactors = F)                         

##Only for multiple samples
#Replicate variant frequency column with sample name
#vaf.col.names <- grep('.vaf', colnames(pyclone.loci), value=T)
#sample.names <- gsub('.vaf', '', vaf.col.names)
#pyclone.loci[,sample.names] <- pyclone.loci[, vaf.col.names]
#vaf.col.names <- sample.names

#pyclone.loci[,sample.names]<-pyclone.loci$variant_allele_frequency
#names(pyclone.loci)[6]<-paste(sample.names,"vaf",sep=".")

#sample.groups<-sample.names
#names(sample.groups)<-vaf.col.names

#Convert variant/cellular frequency to proportions
pyclone.loci$variant_allele_frequency <- pyclone.loci$variant_allele_frequency*100
pyclone.loci$cellular_prevalence <- pyclone.loci$cellular_prevalence*100

#Order data frame by cluster numbering
pyclone.loci<-pyclone.loci[order(pyclone.loci$cluster_id,decreasing = F),]

#Discard clusters with < 3 mutations
cluster.freq<-table(pyclone.loci$cluster_id)

  for (i in 1:length(cluster.freq)) {
    if (cluster.freq[i] < 3) {
      pyclone.loci<-subset(pyclone.loci, subset= cluster_id != i)
    }
  } 

#Add column with normalized cluster ID
  cluster.freq2<-table(pyclone.loci$cluster_id)
  cluster_id_norm<-integer()

  i<-1
  while(i <=  length(cluster.freq2)) {
    cluster_id_norm<-c(cluster_id_norm,rep(i, cluster.freq2[i]))
    i= i+1
    }

  pyclone.loci$cluster<-cluster_id_norm

  clone.tree = infer.clonal.models(variants = pyclone.loci,
                        cluster.col.name = 'cluster',
                        vaf.col.names = 'variant_allele_frequency',
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = order(table(pyclone.loci$cluster), decreasing = T)[1],
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        min.cluster.vaf = 0.01,
                        sum.p = 0.05,
                        alpha = 0.05)
  

clone.tree.branched<-convert.clone.to.branch(clone.tree$matched$merged.trees[[1]], 
                                               branch.lens=NULL, merged.tree.node.annotation = "")
  
clone.tree.top<- clone.tree.branched %>% mutate_all(funs(gsub("#","",.))) %>% select(lab,vaf,parent,ancestors)

#Export loci and  tree files
write.table(pyclone.loci,loci.output, sep="\t",row.names = F, col.names = T, quote = F) 
write.table(clone.tree.top,tree.output, sep="\t",row.names = F, col.names = T, quote = F)
  