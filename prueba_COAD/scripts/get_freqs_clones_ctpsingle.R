#Libraries
library(stringr)

#Input arguments
args<-commandArgs(TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript get_freqs_clones_ctpsingle.R <ctpsingle_cluster_assignments_file>")
}

ctpsingle_input<-args[1]
ctpsingle_output<-args[2]

#ctpsingle_file<-"prueba_COAD/data/ctpsingle/ctpsingle_files/TCGA-3L-AA1B/TCGA-3L-AA1B_ctpsingle_cluster_assignments.txt"

#Functions
##Get_cluster_freq
get_cluster_freq<-function(ctpsingle_clusters, sample_name) {
  
  freq_clust<-table(ctpsingle_clusters$mostLikely)
  total<-nrow(ctpsingle_clusters)
  
  prop_clust<-round((freq_clust/total)*100,2)
  prop_clust<-as.data.frame(prop_clust)
  names(prop_clust)[1]<-"Cluster"
  prop_clust$sample<-rep(sample,nrow(prop_clust))
  
  return(prop_clust)
}


sample=str_match(ctpsingle_input,"TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}")

ctpsingle_clusters<-read.table(ctpsingle_input,
                               header=TRUE, stringsAsFactors = FALSE)

ctpsingle.clust.freq<-get_cluster_freq(ctpsingle_clusters,sample)

write.table(ctpsingle.clust.freq, ctpsingle_output, sep="\t",
            quote=F, row.names = F, col.names = T)

cat("Work finished!\n")