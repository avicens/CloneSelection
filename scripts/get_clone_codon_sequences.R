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

#parent.clone.lab<-pyclone.tree[pyclone.tree$parent==-1,"lab"]
cluster.dndscv.list<-list()

for (i in 1:length(cluster.muts.list)) {
clone.variants<-cluster.muts.list[[i]]
clone.dndscv<-clone.variants[,c("sample_id","Chr","Start_Position","Reference","Alternate")]
clone.annotmuts<-annotate.mutations(clone.dndscv)
clone.annotmuts<-separate(clone.annotmuts,"codonsub",into = c("codon_ref","codon_alt"),sep=">")
#clone.annotmuts$clone<-clone.lab
cluster.dndscv.list<-c(cluster.dndscv.list,list(clone.annotmuts))
}

names(cluster.dndscv.list)<-names(cluster.muts.list)

clones.annotmuts.df<-bind_rows(cluster.dndscv.list, .id = "column_label")
#clone.annotmuts.df<-separate(clones.annotmust.df,"codonsub",into = c("codon_ref","codon_alt"),sep=">")

#Reference sequence
ref.cod<-clones.annotmuts.df[,c("chr","pos","codon_ref")];
ref.cod[,c(1,2)]<-apply(ref.cod[,c(1,2)],2,function(x) as.integer(x))
ref.cod<-ref.cod[with(ref.cod,order(chr,pos)),]
ref.cod$coord<-with(ref.cod,paste(chr,pos,sep=":"))
names(ref.cod)[3]<-"ref"


#Clone sequences
##Parental clone
parental.clones<-pyclone.tree[pyclone.tree$parent == -1, "lab"]

for (p in parental.clones) {
clone.dndscv.df<-cluster.dndscv.list[[p]]

seq.cod<-clone.dndscv.df[,c("chr","pos","codon_alt")]
seq.cod$coord<-with(seq.cod,paste(chr,pos,sep=":"))
seqs.cod<-merge(ref.cod,seq.cod[,c("coord","codon_alt")],by="coord", all.x=T)
gaps<-which(is.na(seqs.cod$codon_alt))
seqs.cod$codon_alt[gaps]<-seqs.cod$ref[gaps]
names(seqs.cod)[ncol(seqs.cod)]<-paste("clone",p,sep="")
}
##Derived clones
###First order clones (Clones directly derived from parental)
#derived.clones.lab<-pyclone.tree[pyclone.tree$parent!=-1,"lab"]
first.order.clones<-pyclone.tree[pyclone.tree$parent == p,"lab"]

for (fo in first.order.clones) { 
  clone.dndscv.df<-cluster.dndscv.list[[fo]]
  
  seq.cod<-clone.dndscv.df[,c("chr","pos","codon_alt")]
  seq.cod$coord<-with(seq.cod,paste(chr,pos,sep=":"))
  seqs.cod<-merge(seqs.cod,seq.cod[,c("coord","codon_alt")],by="coord", all.x=T)
  gaps<-which(is.na(seqs.cod$codon_alt))
  seqs.cod$codon_alt[gaps]<-seqs.cod[gaps,paste("clone",p,sep="")]
  names(seqs.cod)[ncol(seqs.cod)]<-paste("clone",fo,sep="")
  }

###Second order clones
second.order.clones<-pyclone.tree[pyclone.tree$parent != -1 & pyclone.tree$parent != p,"lab"]

for(so in second.order.clones) { 
  ancestral.clone<-pyclone.tree[pyclone.tree$lab == so,"parent"]
  clone.dndscv.df<-cluster.dndscv.list[[so]]
  
  seq.cod<-clone.dndscv.df[,c("chr","pos","codon_alt")]
  seq.cod$coord<-with(seq.cod,paste(chr,pos,sep=":"))
  seqs.cod<-merge(seqs.cod,seq.cod[,c("coord","codon_alt")],by="coord", all.x=T)
  gaps<-which(is.na(seqs.cod$codon_alt))
  seqs.cod$codon_alt[gaps]<-seqs.cod[gaps,paste("clone",ancestral.clone,sep="")]
  names(seqs.cod)[ncol(seqs.cod)]<-paste("clone",so,sep="")
}



if (!dir.exists("data/codon_seqs")) {
  dir.create("data/codon_seqs")
}

seqs.cod<-seqs.cod[order(seqs.cod$chr,seqs.cod$pos,decreasing = F),]
write.fasta(seqs.cod[,4:ncol(seqs.cod)], names = names(seqs.cod)[4:ncol(seqs.cod)], 
            file.out = paste(out.seq.file,sep=""),open= "w")
