#Input arguments
args<-commandArgs(TRUE)
pyclone_out_dir<-args[1]

library(clonevol)
library(dplyr)

pyclone.output.dirs<-list.dirs(path = pyclone_out_dir,recursive = F)

#Function convert.clone.to.branch
convert.clone.to.branch <- function(t, branch.lens = NULL,
                                    merged.tree.node.annotation='sample.with.nonzero.cell.frac.ci'){
  # symbols used for each branch
  syms = c(seq(1,9), unlist(strsplit('abcdefghijklmnopqrstuvwxyz', '')))
  t = t[!is.na(t$parent) & !is.na(t$excluded) & !t$excluded,]
  t$branches = NA
  t$blengths = NA
  rownames(t) = t$lab
  grow.tree <- function(t, lab, parent.symbol=''){
    #print(paste0('---', lab))
    if (t[lab,'parent'] == '-1'){
      t[lab,'branches'] = 'Y'
    }
    children = !is.na(t$parent) & t$parent == lab
    if (any(children)){
      children.labs = t$lab[children]
      num.children = length(children.labs)
      if (num.children == 1){# && parent.symbol != ''){
        # one child, grow straight branch (no left right)
        children.syms = '0'
      }else{
        children.syms = syms[1:num.children]
      }
      children.syms = paste0(parent.symbol, children.syms)
      t$branches[children] = children.syms
      #print(children.labs)
      #print(children.syms)
      for (i in 1:length(children.labs)){
        t = grow.tree(t, children.labs[i], children.syms[i])
      }
    }
    return(t)
  }
  tg = grow.tree(t, t$lab[!is.na(t$parent) & t$parent == '-1'])
  if (merged.tree.node.annotation=='sample.with.nonzero.cell.frac.ci'){
    # remove ci info from sample annotation
    tg$samples.with.nonzero.cell.frac = gsub(',+$', '',
                                             gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.nonzero.cell.frac.ci))
    #tg$samples.with.nonzero.cell.frac = gsub(',+$', '', gsub('\u00B0[^,]+(,|$)', '',
    #     gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.cell.frac.ci)))
  }else{
    cat(paste0('WARN: merged.tree.node.annotation = ',
               merged.tree.node.annotation, ' not supported! No node annotation made.\n'))
  }
  if (is.null(branch.lens)){
    tg$blengths = 5
  }else{
    tg$blengths = branch.lens[tg$lab]
  }
  # color founding clone of met with diff. border
  #tg$node.border.color = ifelse(
  #    grepl('*', gsub('*P', '', tg[[merged.tree.node.annotation]], fixed=TRUE), fixed=T),
  #    'red', 'black')
  tg$node.border.color = 'black'
  tg$node.border.width = 1
  tg$branch.border.color = 'white'
  tg$branch.border.linetype = 'solid'
  tg$branch.border.width = 0.5
  
  return(tg)
  
}


for (i in 1:length(pyclone.output.dirs)) {

  sample<-sapply(strsplit(pyclone.output.dirs[i],"/"),"[",4)
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
  
  clone.tree.top<- clone.tree.branched %>% mutate_all(funs(gsub("#","",.))) %>% select(lab,vaf,parent,ancestors)

#  clone.tree.branch <- convert.consensus.tree.clone.to.branch(clone.tree, 
 #                                                             cluster.col = "cluster", branch.scale = "sqrt")

  
  if (!file.exists(paste(pyclone.output.dir,"tree",sep="/"))) {
    dir.create(paste(pyclone.output.dir,"tree",sep="/"))
  }

  write.table(clone.tree.branched,paste(pyclone.output.dir,"/tree/",sample,"_pyclone_tree.tsv",sep=""),
              sep="\t",row.names = F, col.names = T, quote = F)
  write.table(clone.tree.top,paste(pyclone.output.dir,"/tree/",sample,"_pyclone_tree_topology.tsv",sep=""),
              sep="\t",row.names = F, col.names = T, quote = F)
}
