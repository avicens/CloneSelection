#Load libraries
library(dplyr)
library(tidyr)
library(dndscv)
library(seqinr)


#Input arguments
args<-commandArgs(TRUE)
if (length(args) != 5) {
  cat("Usage: Rscript get_clone_sequences.R <in SNV file> <in Loci file>
      <in Tree file> <out Fasta file> <script.path>\n")
}


in.snv.file=args[1]
in.loci.file=args[2]
in.tree.file=args[3]
out.seq.file=args[4]
script.path=args[5]

#Load functions
source(paste(script.path,"/annotate_mutations.R",sep=""))

#Input files
##Variants file
snv<-read.table(gzfile(in.snv.file),header=T)

##Table with variant information (ClonEvol input)
#pyclone.loci<-read.table(paste("data/pyclone_output",sample,"tables/loci_processed.tsv",sep="/"), sep="\t", header=T, stringsAsFactors = F)
loci<-read.table(in.loci.file, sep="\t", header=T, stringsAsFactors = F)

##Table with clone tree inferred by ClonEvol using PyClone output
#pyclone.tree<-read.table(paste("data/clone_trees/",sample,"_tree.txt",sep=""), sep="\t", header=T, stringsAsFactors = F)
clone.tree<-read.table(in.tree.file, sep="\t", header=T, stringsAsFactors = F)

#Merge snv and loci tables by mutation code
snv$mutation_id<-with(snv,paste(Hugo_Symbol,":",Chr,"-",Start_Position,"_",End_Position,sep=""))
muts<-merge(snv,pyclone.loci,by="mutation_id")
cluster.muts.list<-split(muts, muts$cluster)

#parent.clone.lab<-pyclone.tree[pyclone.tree$parent==-1,"lab"]
cluster.dndscv.list<-list()

for (i in 1:length(cluster.muts.list)) {
clone.variants<-cluster.muts.list[[i]]
clone.variants2<-clone.variants[,c("sample_id","Chr","Start_Position","Reference","Alternate")]
clone.annotmuts<-annotate.mutations(clone.dndscv)
clone.annotmuts2<-separate(clone.annotmuts,"codonsub",into = c("codon_ref","codon_alt"),sep=">")
cluster.dndscv.list<-c(cluster.dndscv.list,list(clone.annotmuts))
}

names(cluster.dndscv.list)<-names(cluster.muts.list)

clones.annotmuts.df<-bind_rows(cluster.dndscv.list, .id = "column_label")

#Reference sequence
ref.cod<-clones.annotmuts.df[,c("chr","pos","codon_ref")];
ref.cod[,c(1,2)]<-apply(ref.cod[,c(1,2)],2,function(x) as.integer(x))
ref.cod<-ref.cod[with(ref.cod,order(chr,pos)),]
ref.cod$coord<-with(ref.cod,paste(chr,pos,sep=":"))
names(ref.cod)[3]<-"ref"


#Clone sequences
##Parental clone
parental.clone<-clone.tree[clone.tree$parent == -1, "lab"]

for (p in parental.clone) {
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
first.order.clones<-clone.tree[clone.tree$parent == parental.clone,"lab"]

for (fo in first.order.clones) { 
  clone.dndscv.df<-cluster.dndscv.list[[fo]]
  
  seq.cod<-clone.dndscv.df[,c("chr","pos","codon_alt")]
  seq.cod$coord<-with(seq.cod,paste(chr,pos,sep=":"))
  seqs.cod<-merge(seqs.cod,seq.cod[,c("coord","codon_alt")],by="coord", all.x=T)
  gaps<-which(is.na(seqs.cod$codon_alt))
  seqs.cod$codon_alt[gaps]<-seqs.cod[gaps,paste("clone",parental.clone,sep="")]
  names(seqs.cod)[ncol(seqs.cod)]<-paste("clone",fo,sep="")
  }

###Second order clones
second.order.clones<-clone.tree[clone.tree$parent != -1 & clone.tree$parent != parental.clone,"lab"]

for(so in second.order.clones) { 
  ancestral.clone<-clone.tree[clone.tree$lab == so,"parent"]
  clone.dndscv.df<-cluster.dndscv.list[[so]]
  
  seq.cod<-clone.dndscv.df[,c("chr","pos","codon_alt")]
  seq.cod$coord<-with(seq.cod,paste(chr,pos,sep=":"))
  seqs.cod<-merge(seqs.cod,seq.cod[,c("coord","codon_alt")],by="coord", all.x=T)
  gaps<-which(is.na(seqs.cod$codon_alt))
  seqs.cod$codon_alt[gaps]<-seqs.cod[gaps,paste("clone",ancestral.clone,sep="")]
  names(seqs.cod)[ncol(seqs.cod)]<-paste("clone",so,sep="")
}


if (!file.exists("data/seqs")) {
  dir.create("data/seqs")
}

seqs<-seqs[order(seqs$chr,seqs$pos,decreasing = F),]
write.fasta(seqs[,4:ncol(seqs)], names = names(seqs)[4:ncol(seqs)], 
            file.out = paste(out.seq.file,sep=""),open= "w")
