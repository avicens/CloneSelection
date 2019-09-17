library(dplyr)
library(ggplot2)

#Load data
##Clusters
clusters1<-read.table("prueba_COAD/tests/TEST1/cluster_frequencies.tsv",sep= "\t", header=T)
clusters4<-read.table("prueba_COAD/tests/TEST4/cluster_frequencies.tsv",sep= "\t", header=T)
#Discard clusters with a single mutation (size =1)
clusters1.filt <- clusters1 %>% filter(size > 1)
clusters4.filt <- clusters4 %>% filter(size > 1)

##Loci
loci1<-read.table("prueba_COAD/tests/TEST1/muts_clones_per_samples.tsv",sep="\t",header=T)
loci4<-read.table("prueba_COAD/tests/TEST4/muts_clones_per_samples.tsv",sep="\t",header=T)

#Bind tables
clusters<-rbind(clusters1.filt, clusters4.filt)
clusters$test<-as.factor(c(rep(1,nrow(clusters1.filt)),rep(4,nrow(clusters4.filt))))

loci<-rbind(loci1, loci4)
loci$test<-as.factor(c(rep(1,nrow(loci1)),rep(4,nrow(loci4))))

#Plots

##Mutations per sample
plot.mutations<-ggplot(loci,aes(x=mutations, fill=test)) + geom_histogram(alpha=0.3, color ="black", position="dodge") + theme_bw() + labs(x="Mutations per sample", y ="Count")

#Clusters per sample
plot.clusters<-ggplot(loci,aes(x=clusters, fill=test)) + geom_histogram(alpha=0.3, color ="black", position="dodge", binwidth = 1) + theme_bw() + labs(x="Clusters per sample", y ="Count")

##Cluster size
plot.cluster.size<-ggplot(subset(clusters,size < 500),aes(x=size, color=test)) + geom_density(alpha=0.5) + theme_bw() + labs(x="Cluster size", y = "Density")

##Cluster frequency (cellular prevalence)
plot.cluster.freq<-ggplot(subset(clusters),aes(x=mean, color=test)) + geom_density(alpha=0.5) + theme_bw() + labs(x="Cellular prevalence", y = "Density")

ggarrange(plot.mutations, plot.cluster.size, plot.clusters, plot.cluster.freq, ncol=2, nrow=2, common.legend = T)

