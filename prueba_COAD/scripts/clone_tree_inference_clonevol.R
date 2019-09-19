#Load libraries
library(clonevol)
library(dplyr)
library(stringr)


#Input arguments
args<-commandArgs(TRUE)

if (length(args) != 3) {
  cat("Usage: Rscript clone_tree_inference_clonevol.R  [1]<clonevol.input> [2]<tree.output> 
      [3]<path to convert.clone.to.branch script>\n")
}

clonevol.input<-args[1]
tree.output<-args[2]
script.path<-args[3]

#Load function for converting clone trees to branch
source(paste(script.path,"/convert.clone.to.branch.R",sep=""))

loci<-read.table(clonevol.input, header=T, stringsAsFactors = F)
  
  clone.tree = infer.clonal.models(variants = loci,
                        cluster.col.name = 'cluster',
                        vaf.col.names = 'variant_allele_frequency',
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = order(table(loci$cluster), decreasing = T)[1],
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        min.cluster.vaf = 0.01,
                        sum.p = 0.05,
                        alpha = 0.05)
  

  clone.tree.branched<-convert.clone.to.branch(clone.tree$matched$merged.trees[[1]], 
                                               branch.lens=NULL, merged.tree.node.annotation = "")
  
  clone.tree.top<- clone.tree.branched %>% mutate_all(funs(gsub("#","",.))) %>% select(lab,vaf,parent,ancestors)

  #Save files
write.table(clone.tree.top,tree.output, sep="\t",row.names = F, col.names = T, quote = F)
