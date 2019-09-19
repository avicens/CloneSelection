#Load libraries
library(dplyr)
library(stringr)

args<-commandArgs(TRUE)
pyclone.loci<-args[1]
clonevol.input<-args[2]

pyclone.loci<-read.table(pyclone.loci,
                         sep="\t", header=T, stringsAsFactors = F)


#Convert variant/cellular frequency to proportions
pyclone.loci$variant_allele_frequency <- pyclone.loci$variant_allele_frequency*100
pyclone.loci$cellular_prevalence <- pyclone.loci$cellular_prevalence*100

#Order data frame by cluster numbering
pyclone.loci<-pyclone.loci[order(pyclone.loci$cluster_id,decreasing = F),]

#Discard clusters with < 3 mutations
cluster.freq<-table(pyclone.loci$cluster_id)
idx<-names(which(cluster.freq >2))
pyclone.loci2<-pyclone.loci[pyclone.loci$cluster_id %in% idx,]

#Add column with normalized cluster ID
cluster.freq2<-table(pyclone.loci2$cluster_id)
cluster_id_norm<-integer()

i<-1
while(i <=  length(cluster.freq2)) {
  cluster_id_norm<-c(cluster_id_norm,rep(i, cluster.freq2[i]))
  i= i+1
}

pyclone.loci2$cluster<-cluster_id_norm

write.table(pyclone.loci2,clonevol.input, sep="\t",row.names = F, col.names = T, quote = F) 





