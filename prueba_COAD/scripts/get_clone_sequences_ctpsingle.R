#Input arguments
args<-commandArgs(TRUE)

if (length(args) != 2) {
  cat("Usage Rscript get_clone_sequences_ctpsingle.R <sample_ctpsingle_directory> <sequence_directory>\n")
  cat("Exiting\n")
  quit()
}
#Load libraries
library(dplyr)
library(tidyr)
library(seqinr)
library(stringr)

#Assigning inputs to objects
sample.dir=args[1]
seq.dir=args[2]

sample=str_match(sample.dir,"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")

#Load functions
source(paste("scripts/annotate_mutations.R",sep=""))

inferseq<-function(ref.pos, new.pos) {
  seqs<-new.pos %>% select(chr,pos,codon_alt) %>%
    mutate(coord=paste(chr,pos,sep=":")) %>% select(codon_alt,coord) %>% 
    right_join(ref.pos, by = "coord") 
  
  gaps<-which(is.na(seqs$codon_alt))
  seqs$codon_alt[gaps]<-seqs$codon[gaps]
  seqs<-seqs[seqs$codon_alt !=".",]
  inf.seq<-seqs[,c("coord","codon_alt")]
  names(inf.seq)[2] <- "codon"
  
  return(inf.seq)
} 

#Read cluster assignment file
cluster.assign.file<-read.table(paste(sample.dir,"/",sample,"_ctpsingle_cluster_assignments.txt",sep=""),header=T)

cluster.assign.file2<-cluster.assign.file %>% 
  select(c("Chromosome","Position","Mutant","Reference","mostLikely")) %>% 
  arrange(mostLikely) %>% 
  mutate(sample_id=rep(sample,nrow(cluster.assign.file))) %>%
  select(c("sample_id","Chromosome","Position","Reference","Mutant","mostLikely")) 


cluster.assign.list<-split(cluster.assign.file2, cluster.assign.file2$mostLikely)

cluster.annots.list<-list()

for (i in 1:length(cluster.assign.list)) {
  cluster.annotmuts<-annotate.mutations(cluster.assign.list[[i]])
  cluster.annotmuts<-separate(cluster.annotmuts,"codonsub",into = c("codon_ref","codon_alt"),sep=">")
  cluster.annots.list<-c(cluster.annots.list,list(cluster.annotmuts))
}

names(cluster.annots.list)<-names(cluster.assign.list)
clusters.annotmuts.df<-bind_rows(cluster.annots.list, .id="column_label")

#Reconstruct sequences

##Reference sequence
ref.seq<-clusters.annotmuts.df %>% 
select(chr,pos,codon_ref) %>% 
mutate(chr=as.integer(chr), pos=as.integer(pos)) %>%
arrange(chr,pos) %>% mutate(coord=paste(chr,pos,sep=":")) %>% filter(codon_ref !=".")
ref.seq<-ref.seq[!duplicated(ref.seq$coord),] 

names(ref.seq)[3]<-"codon"

##Crear loop para analizar todos los árboles
topologies<-list.files(sample.dir, pattern = "*num*")

#Ejecutar loop solo para archivos con score = 0
for (top in topologies) {
  
  clone.tree<-read.table(paste(sample.dir,top,sep="/"),header=F) 
  clone.tree= clone.tree %>% select((ncol(clone.tree)-3):ncol(clone.tree))
  names(clone.tree)=c("parent_node","child_node","cell_fraction_child_node","score")

  tip.nodes<-clone.tree %>% filter(!(child_node %in% parent_node)) %>% select(child_node)
  
  if ( clone.tree[1,"score"] > 0) {
    cat("Topology ", top, "discarded for having score > 0\n")
    next
  } else {
    cat("Inferring sequences for topology",top,"\n")
    tree<-str_match(top,"tree_[0-99]")
    
    #Parental clone
    parental.clone.lab<-clone.tree[clone.tree$parent_node == 0, "child_node"]
    parental.clone.muts<-cluster.annots.list[[parental.clone.lab]]
    parental.clone.muts<-parental.clone.muts[!duplicated(parental.clone.muts$pos) & !duplicated(parental.clone.muts$chr),]
    
    assign(paste("clone",parental.clone.lab,"seq",sep=""),
           inferseq(ref.seq,parental.clone.muts))
    
    parent.node=parental.clone.lab
    child.nodes=clone.tree[clone.tree$parent_node == parent.node, "child_node"]
    
    while (length(child.nodes) > 0) {
      
      for (node in child.nodes) {
        child.clone.muts<-cluster.annots.list[[node]]
        child.clone.muts<-child.clone.muts[!duplicated(child.clone.muts$pos) & !duplicated(child.clone.muts$chr),]
        assign(paste("clone",node,"seq",sep=""),
               inferseq(get(paste("clone",parent.node,"seq",sep="")),
                        child.clone.muts))
      }
      
      parent.nodes=child.nodes
      child.nodes=clone.tree[clone.tree$parent_node %in% parent.nodes, "child_node"]
      
    }
    
    #Bind the sequences
    cat("Binding sequences...\n")
    clone.seqs<-apropos("clone[1-99]seq")
    all.seqs<-ref.seq %>% select(1:4)
    names(all.seqs)[3]<-"ref"
    
    cl<-1
    while (cl <= length(clone.seqs)) {
      
      cl.seq<-get(clone.seqs[cl])
      all.seqs<-merge(all.seqs, cl.seq, by="coord", all.x=TRUE)
      names(all.seqs)[ncol(all.seqs)]<-gsub("seq","",clone.seqs[cl])
      cl=cl+1
    }
    
    #Convert to fasta
    cat("Saving all sequences in FASTA format...\n")
    
    all.seq.file=paste(seq.dir,"/",sample,"_",tree,"_clone_seqs.fas",sep="")
    
    all.seqs.cod<-all.seqs %>% arrange(chr,pos) %>% select(ref:ncol(all.seqs))
    
    write.fasta(all.seqs.cod, names=names(all.seqs.cod),
                file.out=all.seq.file, open="w")
    
    #Sequences for tip clones
    cat("Saving tip sequences in FASTA format...\n")
    tip.nodes.labels<-sapply(tip.nodes$child_node, function(x) paste("clone",x,sep=""))
    tip.seqs.cod<- all.seqs.cod[,c("ref",tip.nodes.labels)]
    
    tip.seq.file=paste(seq.dir,"/",sample,"_",tree,"_clone_tip_seqs.fas",sep="")
    
    write.fasta(tip.seqs.cod, names=names(tip.seqs.cod),
                file.out=tip.seq.file, open="w")
    
    rm(list=clone.seqs)
  }
}

cat("Work finished successfully!\n")
