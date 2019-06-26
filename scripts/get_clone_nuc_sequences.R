#Load libraries
library(dplyr)
library(tidyr)
library(dndscv)
library(seqinr)

#Load functions
source("scripts/annotate_mutations.R")

#Input arguments
args<-commandArgs(TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript get_clone_sequences.R <in SNV file> <in Loci file>
      <in Tree file> <out Fasta file>\n")
}


in.snv.file=args[1]
in.loci.file=args[2]
in.tree.file=args[3]
out.seq.file=args[4]

#Input files
##Variants file
snv<-read.table(gzfile(in.snv.file),header=T)

##Table with loci information from PyClone
#pyclone.loci<-read.table(paste("data/pyclone_output",sample,"tables/loci_processed.tsv",sep="/"), sep="\t", header=T, stringsAsFactors = F)
pyclone.loci<-read.table(in.loci.file, sep="\t", header=T, stringsAsFactors = F)

##Table with clone tree inferred by ClonEvol using PyClone output
#pyclone.tree<-read.table(paste("data/clone_trees/",sample,"_tree.txt",sep=""), sep="\t", header=T, stringsAsFactors = F)
pyclone.tree<-read.table(in.tree.file, sep="\t", header=T, stringsAsFactors = F)

#Merge snv and loci tables by mutation code
snv$mutation_id<-with(snv,paste(Hugo_Symbol,":",Chr,"-",Start_Position,"_",End_Position,sep=""))
muts<-merge(snv,pyclone.loci,by="mutation_id")
cluster.muts.list<-split(muts, muts$cluster)
clones.muts.df<-bind_rows(cluster.muts.list, .id = "column_label")

#Reference sequence
#ref.cod<-clones.annotmuts.df[,c("chr","pos","codon_ref")];
ref.nuc<-clones.muts.df[,c("Chr","Start_Position","Reference")]
ref.nuc<-ref.nuc[with(ref.nuc,order(Chr,Start_Position)),]
ref.nuc$coord<-with(ref.nuc,paste(Chr,Start_Position,sep=":"))
names(ref.nuc)[3]<-"ref"


#Clone sequences
##Parental clone
parental.clones<-pyclone.tree[pyclone.tree$parent == -1, "lab"]

for (p in parental.clones) {
clone.muts.df<-cluster.muts.list[[p]]

seq.nuc<-clone.muts.df[,c("Chr","Start_Position","Alternate")]
seq.nuc$coord<-with(seq.nuc,paste(Chr,Start_Position,sep=":"))
seqs.nuc<-merge(ref.nuc,seq.nuc[,c("coord","Alternate")],by="coord", all.x=T)
gaps<-which(is.na(seqs.nuc$Alternate))
seqs.nuc$Alternate[gaps]<-seqs.nuc$ref[gaps]
names(seqs.nuc)[ncol(seqs.nuc)]<-paste("clone",p,sep="")
}
##Derived clones
###First order clones (Clones directly derived from parental)
#derived.clones.lab<-pyclone.tree[pyclone.tree$parent!=-1,"lab"]
first.order.clones<-pyclone.tree[pyclone.tree$parent == p,"lab"]

for (fo in first.order.clones) { 
  clone.muts.df<-cluster.muts.list[[fo]]
  
  seq.nuc<-clone.muts.df[,c("Chr","Start_Position","Alternate")]
  seq.nuc$coord<-with(seq.nuc,paste(Chr,Start_Position,sep=":"))
  seqs.nuc<-merge(seqs.nuc,seq.nuc[,c("coord","Alternate")],by="coord", all.x=T)
  gaps<-which(is.na(seqs.nuc$Alternate))
  seqs.nuc$Alternate[gaps]<-seqs.nuc[gaps,paste("clone",p,sep="")]
  names(seqs.nuc)[ncol(seqs.nuc)]<-paste("clone",fo,sep="")
  }

###Second order clones
second.order.clones<-pyclone.tree[pyclone.tree$parent != -1 & pyclone.tree$parent != p,"lab"]

for(so in second.order.clones) { 
  ancestral.clone<-pyclone.tree[pyclone.tree$lab == so,"parent"]
  clone.muts.df<-cluster.muts.list[[so]]
  
  seq.nuc<-clone.muts.df[,c("Chr","Start_Position","Alternate")]
  seq.nuc$coord<-with(seq.nuc,paste(Chr,Start_Position,sep=":"))
  seqs.nuc<-merge(seqs.nuc,seq.nuc[,c("coord","Alternate")],by="coord", all.x=T)
  gaps<-which(is.na(seqs.nuc$Alternate))
  seqs.nuc$Alternate[gaps]<-seqs.nuc[gaps,paste("clone",ancestral.clone,sep="")]
  names(seqs.nuc)[ncol(seqs.nuc)]<-paste("clone",so,sep="")
}


if (!dir.exists("data/nuc_seqs")) {
  dir.create("data/nuc_seqs")
}

seqs.nuc<-seqs.nuc[order(seqs.nuc$Chr,seqs.nuc$Start_Position,decreasing = F),]
write.fasta(seqs.nuc[,4:ncol(seqs.nuc)], names = names(seqs.nuc)[4:ncol(seqs.nuc)], 
            file.out = paste(out.seq.file,sep=""),open= "w")
