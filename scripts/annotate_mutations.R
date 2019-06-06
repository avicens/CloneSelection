annotate.mutations = function(mutations,  refdb = "hg19", gene_list = NULL,numcode = 1) {   
 ## 1. Environment
  message("[1] Loading the environment...")
  
 #Loading genome reference
   #[Input] Reference database
  if (refdb == "hg19") {
    data("refcds_hg19", package="dndscv")
    if (any(gene_list=="CDKN2A")) { # Replace CDKN2A in the input gene list with two isoforms
      gene_list = unique(c(setdiff(gene_list,"CDKN2A"),"CDKN2A.p14arf","CDKN2A.p16INK4a"))
    }
  } else {
    load(refdb)
  }
  
  # Expanding the reference sequences [for faster access]
  for (j in 1:length(RefCDS)) {
    RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
    RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
    RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
    if (!is.null(RefCDS[[j]]$seq_splice)) {
      RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
      RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
      RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
    }
  }
  
   ## 2. Load mutations
  mutations = mutations[,1:5] # Restricting input matrix to first 5 columns
  mutations[,c(1,2,3,4,5)] = lapply(mutations[,c(1,2,3,4,5)], as.character) # Factors to character
  mutations[[3]] = as.numeric(mutations[[3]]) # Chromosome position as numeric
  mutations = mutations[mutations[,4]!=mutations[,5],] # Removing mutations with identical reference and mutant base
  colnames(mutations) = c("sampleID","chr","pos","ref","mut")
  
  # Removing NA entries from the input mutation table
  indna = which(is.na(mutations),arr.ind=T)
  if (nrow(indna)>0) {
    mutations = mutations[-unique(indna[,1]),] # Removing entries with an NA in any row
    warning(sprintf("%0.0f rows in the input table contained NA entries and have been removed. Please investigate.",length(unique(indna[,1]))))
  }
  
      ## 3. Mutation annotation
  message("[2] Annotating the mutations...")
  
  nt = c("A","C","G","T")
  trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
  trinucinds = setNames(1:64, trinucs)
  trinucsubs = NULL
  for (j in 1:length(trinucs)) {
    trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
  }
  trinucsubsind = setNames(1:192, trinucsubs)
  
  ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
  gr_genes_ind = ind[gr_genes$names]
  
  # Warning about possible unannotated dinucleotide substitutions
  if (any(diff(mutations$pos)==1)) {
    warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
  }
  
  # Warning about multiple instances of the same mutation in different sampleIDs
  if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
    warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
  }
  
  # Start and end position of each mutation
  mutations$end = mutations$start = mutations$pos
  l = nchar(mutations$ref)-1 # Deletions of multiple bases
  mutations$end = mutations$end + l
  ind = substr(mutations$ref,1,1)==substr(mutations$mut,1,1) & nchar(mutations$ref)>nchar(mutations$mut) # Position correction for deletions annotated in the previous base (e.g. CA>C)
  mutations$start = mutations$start + ind
  
  # Mapping mutations to genes
  gr_muts = GenomicRanges::GRanges(mutations$chr, IRanges::IRanges(mutations$start,mutations$end))
  ol = as.data.frame(GenomicRanges::findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene
  mutations$geneind = gr_genes_ind[ol[,2]]
  mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]
  mutations = unique(mutations)
  
  # Additional annotation of substitutions
  
  mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
  snv = (mutations$ref %in% nt & mutations$mut %in% nt)
  if (!any(snv)) { stop("Zero coding substitutions found in this dataset. Unable to run dndscv.") }
  indels = mutations[!snv,]
  mutations = mutations[snv,]
  mutations$ref_cod = mutations$ref
  mutations$mut_cod = mutations$mut
  compnt = setNames(rev(nt), nt)
  isminus = (mutations$strand==-1)
  mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
  mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]
  
 # for (j in 1:length(RefCDS)) {
#    RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
#  }
  
  # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
  
  chr2cds = function(pos,cds_int,strand) {
    if (strand==1) {
      return(which(unlist(apply(cds_int, 1, function(x) x[1]:x[2])) %in% pos))
    } else if (strand==-1) {
      return(which(rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2]))) %in% pos))
    }
  }
  
  # Annotating the functional impact of each substitution and populating the N matrices
  
  ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = codonsub = array(NA, nrow(mutations))
  
  for (j in 1:nrow(mutations)) {
    
    geneind = mutations$geneind[j]
    pos = mutations$pos[j]
    
    if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution
      
      impact[j] = "Essential_Splice"; impind = 4
      pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
      cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
      aachange[j] = ntchange[j] = codonsub[j] = "."
      
    } else { # Coding substitution
      
      pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
      cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
      ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
      codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
      old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
      pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
      new_codon = old_codon; new_codon[pos_in_codon] = mutations$mut_cod[j]
      old_aa = seqinr::translate(old_codon, numcode = numcode)
      new_aa = seqinr::translate(new_codon, numcode = numcode)
      aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
      ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
      codonsub[j] = sprintf('%s>%s',paste(old_codon,collapse=""),paste(new_codon,collapse=""))
      
      # Annotating the impact of the mutation
      if (new_aa == old_aa){ 
        impact[j] = "Synonymous"; impind = 1
      } else if (new_aa == "*"){
        impact[j] = "Nonsense"; impind = 3
      } else if (old_aa != "*"){
        impact[j] = "Missense"; impind = 2
      } else if (old_aa=="*") {
        impact[j] = "Stop_loss"; impind = NA
      }
    }
    
    if (mutations$ref_cod[j] != as.character(cdsnt)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
      wrong_ref[j] = 1
    } else if (!is.na(impind)) { # Correct base annotation in the input mutation file
      trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
      RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices
    }
    
    if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g%% ...', round(j/nrow(mutations),2)*100)) }
  }
  
  mutations$ref3_cod = ref3_cod
  mutations$mut3_cod = mut3_cod
  mutations$aachange = aachange
  mutations$ntchange = ntchange
  mutations$codonsub = codonsub
  mutations$impact = impact
  mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]
  return(mutations)
  }
