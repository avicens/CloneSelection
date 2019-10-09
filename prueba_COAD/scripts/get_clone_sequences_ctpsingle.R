library(dplyr)
library(seqinr)
#Load functions
source(paste(script.path,"/annotate_mutations.R",sep=""))

sample="TCGA-3L-AA1B"

inferseq<-function(ref.pos, new.pos) {
  seqs<-new.pos %>% select(chr,pos,codon_alt) %>%
    mutate(coord=paste(chr,pos,sep=":")) %>% select(codon_alt,coord) %>% 
    right_join(ref.pos, by = "coord") 
  
  gaps<-which(is.na(seqs$codon_alt))
  seqs$codon_alt[gaps]<-seqs$codon[gaps]
  
  inf.seq<-seqs[,c("coord","codon_alt")]
  names(inf.seq)[2] <- "codon"
  
  return(inf.seq)
} 


cluster.assign.file<-read.table("prueba_COAD/data/ctpsingle/ctpsingle_files/TCGA-3L-AA1B/TCGA-3L-AA1B_ctpsingle_cluster_assignments.txt",
                                header=T)

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
#Reference sequence
ref.seq<-clusters.annotmuts.df %>% 
  select(chr,pos,codon_ref) %>% 
  mutate(chr=as.integer(chr), pos=as.integer(pos)) %>%
  arrange(chr,pos) %>% 
  mutate(coord=paste(chr,pos,sep=":")) %>%  
  rename(codon = codon_ref)

##Crear loop para analizar todos los árboles
topologies<-list.files(path=paste("prueba_COAD/data/ctpsingle/ctpsingle_files/",sample, sep=""),
                       pattern = "*num*")

#Ejecutar loop solo para archivos con score = 0
for (top in topologies) {
  
  clone.tree<-read.table(paste("prueba_COAD/data/ctpsingle/ctpsingle_files",sample,top,sep="/"),header=F) %>% select(6:9)
  names(clone.tree)=c("parent_node","child_node","cell_fraction_child_node","score")
  if ( clone.tree[1,"score"] > 0) {
    cat("Topology ", top, "discarded for having score > 0\n")
    next
  } else {
    cat("Inferring sequences for topology",top,"\n")
    tree<-str_match(top,"tree_[0-99]")
    
    #Parental clone
    parental.clone.lab<-clone.tree[clone.tree$parent_node == 0, "child_node"]
    parental.clone.muts<-cluster.annots.list[[parental.clone.lab]]
    
    assign(paste("clone",parental.clone.lab,"seq",sep=""),
           inferseq(par.seq,clone.annotmuts.df))
    
    parent.node=parental.clone.lab
    child.nodes=clone.tree[clone.tree$parent_node == parent.node, "child_node"]
    
    while (length(child.nodes) > 0) {
      
      for (node in child.nodes) {
        child.clone.muts<-cluster.annots.list[[node]]
        assign(paste("clone",node,"seq",sep=""),
               inferseq(get(paste("clone",parent.node,"seq",sep="")),
                        child.clone.muts))
      }
      
      parent.node=child.nodes
      child.nodes=clone.tree[clone.tree$parent_node == parent.node, "child_node"]
      
    }
    
    #Bind the sequences
    cat("Binding sequences...\n")
    clone.seqs<-apropos("clone[1-99]seq")
    all.seqs<-ref.seq %>% select(1:3) %>% rename(ref=codon)
    
    cl<-1
    while (cl <= length(clone.seqs)) {
      
      cl.seq<-as.character(get(clone.seqs[cl])[,"codon"])
      all.seqs<-cbind(all.seqs, cl.seq)
      names(all.seqs)[ncol(all.seqs)]<-gsub("seq","",clone.seqs[cl])
      cl=cl+1
    }
    
    #Convert to fasta
    cat("Saving sequences in FASTA format...\n")
    
    if (!file.exists(paste("prueba_COAD/data/ctpsingle/seqs/",sample,sep=""))) {
      dir.create(paste("prueba_COAD/data/ctpsingle/seqs/",sample,sep=""))
    }
    
    seq.file=paste("prueba_COAD/data/ctpsingle/seqs/",sample,"/",sample,"_",tree,"_clone_seqs.fas",sep="")
    
    all.seqs.cod<-all.seqs %>% arrange(chr,pos) %>% select(ref:ncol(all.seqs))
    write.fasta(all.seqs.cod, 
                names=names(all.seqs.cod),
                file.out=seq.file,open="w")
    
    rm(list=clone.seqs)
  }
}
