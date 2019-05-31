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
  cat("Usage: Rscript get_clone_sequences.R  <in SNV file> <in Loci file>
      <in Tree file> <out Fasta file>\n")
}


in.snv.file=args[1]
in.loci.file=args[2]
in.tree.file=args[3]
out.seq.file=args[4]

#Input files
##1)Variants file
#snv<-read.table(gzfile(paste("data/snv/",sample,"_SNV.tsv.gz",sep="")),header=T)
snv<-read.table(gzfile(in.snv.file),header=T)
##Table with loci information from PyClone
#pyclone_out_dir<-"data/pyclone_output"
#pyclone.loci<-read.table(paste("data/pyclone_output",sample,"tables/loci_processed.tsv",sep="/"), sep="\t", header=T, stringsAsFactors = F)
pyclone.loci<-read.table(in.loci.file, sep="\t", header=T, stringsAsFactors = F)

##Table with clone tree inferred by ClonEvol using PyClone output
#pyclone.tree<-read.table(paste("data/trees/",sample,"_tree.txt",sep=""), sep="\t", header=T, stringsAsFactors = F)
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
#clone.dndsout<-dndscv(clone.dndscv, outp = 1)
#clone.annotmuts<-clone.dndsout$annotmuts
clone.annotmuts<-annotate.mutations(clone.dndscv)
clone.annotmuts<-separate(clone.annotmuts,"codonsub",into = c("codon_ref","codon_alt"),sep=">")
#clone.annotmuts$clone<-clone.lab
cluster.dndscv.list<-c(cluster.dndscv.list,list(clone.annotmuts))
}

names(cluster.dndscv.list)<-names(cluster.muts.list)

clones.annotmuts.df<-bind_rows(cluster.dndscv.list, .id = "column_label")
#clone.annotmuts.df<-separate(clones.annotmust.df,"codonsub",into = c("codon_ref","codon_alt"),sep=">")

#Reference sequence
ref<-clones.annotmuts.df[,c("chr","pos","codon_ref")];
ref[,c(1,2)]<-apply(ref[,c(1,2)],2,function(x) as.integer(x))
ref<-ref[with(ref,order(chr,pos)),]
ref$coord<-with(ref,paste(chr,pos,sep=":"))
names(ref)[3]<-"ref"


#Clone sequences
##Parental clone
parental.clone.lab<-pyclone.tree[pyclone.tree$parent == -1, "lab"]

for (i in parental.clone.lab) {

clone.dndscv.df<-cluster.dndscv.list[[i]]
seq<-clone.dndscv.df[,c("chr","pos","codon_alt")]
seq$coord<-with(seq,paste(chr,pos,sep=":"))
seqs<-merge(ref,seq[,c("coord","codon_alt")],by="coord", all.x=T)
gaps<-which(is.na(seqs$codon_alt))
seqs$codon_alt[gaps]<-seqs$ref[gaps]
names(seqs)[ncol(seqs)]<-paste("clone",i,sep="")

}

##Derived clones
###First order clones (Clones directly derived from parental)
#derived.clones.lab<-pyclone.tree[pyclone.tree$parent!=-1,"lab"]
first.order.clones<-pyclone.tree[pyclone.tree$parent == parental.clone.lab,"lab"]

for(i in first.order.clones) { 
  clone.dndscv.df<-cluster.dndscv.list[[i]]
  seq<-clone.dndscv.df[,c("chr","pos","codon_alt")]
  seq$coord<-with(seq,paste(chr,pos,sep=":"))
  seqs<-merge(seqs,seq[,c("coord","codon_alt")],by="coord", all.x=T)
  gaps<-which(is.na(seqs$codon_alt))
  seqs$codon_alt[gaps]<-seqs[gaps,paste("clone",parental.clone.lab,sep="")]
  names(seqs)[ncol(seqs)]<-paste("clone",i,sep="")
  }

###Second order clones
second.order.clones<-pyclone.tree[pyclone.tree$parent != -1 & pyclone.tree$parent != parental.clone.lab,"lab"]

for(i in second.order.clones) { 
  ancestral.clone<-pyclone.tree[pyclone.tree$lab == i,"parent"]
  clone.dndscv.df<-cluster.dndscv.list[[i]]
  seq<-clone.dndscv.df[,c("chr","pos","codon_alt")]
  seq$coord<-with(seq,paste(chr,pos,sep=":"))
  seqs<-merge(seqs,seq[,c("coord","codon_alt")],by="coord", all.x=T)
  gaps<-which(is.na(seqs$codon_alt))
  seqs$codon_alt[gaps]<-seqs[gaps,paste("clone",ancestral.clone,sep="")]
  names(seqs)[ncol(seqs)]<-paste("clone",i,sep="")
}

if (!file.exists("data/seqs")) {
  dir.create("data/seqs")
}

seqs<-seqs[order(seqs$chr,seqs$pos,decreasing = F),]
write.fasta(seqs[,4:ncol(seqs)], names = names(seqs)[4:ncol(seqs)], 
            file.out = paste(out.seq.file,sep=""),open= "w")
