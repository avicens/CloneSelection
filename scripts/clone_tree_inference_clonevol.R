library(clonevol)

pyclone.output.dirs<-list.dirs(path = "data/pyclone_output",recursive = F)

for (i in 1:length(pyclone.output.dirs)) {

  sample<-sapply(strsplit(pyclone.output.dirs[i],"/"),"[",3)
  pyclone.output.dir<-pyclone.output.dirs[i]
  pyclone.loci<-read.table(paste(pyclone.output.dir,"tables/loci.tsv",sep="/"),
                         sep="\t", header=T, stringsAsFactors = F)

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

#  clone.tree.branch <- convert.consensus.tree.clone.to.branch(clone.tree, 
 #                                                             cluster.col = "cluster", branch.scale = "sqrt")

  
  if (!file.exists(paste(pyclone.output.dir,"tree",sep="/"))) {
    dir.create(paste(pyclone.output.dir,"tree",sep="/"))
  }

  write.table(clone.tree.branched,paste(pyclone.output.dir,"/tree/",sample,"_pyclone_tree.tsv",sep=""),
              sep="\t",row.names = F, col.names = T, quote = F)
}
