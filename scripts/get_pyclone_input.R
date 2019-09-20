library(dplyr, warn.conflicts = FALSE)
#Input arguments
args<-commandArgs(TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript get_pyclone_input_Script.R  <intersFile> <outFile>\n")
  
}
intersFile<-args[1]
outFile<-args[2]

snv.int.cnv<-read.table(intersFile, stringsAsFactors = F, header=F)
names(snv.int.cnv)<-c("Chr","Start","End",	"Reference",	"Alternate",	"Variant_freq",	"Ref_freq",
                 "Hugo_Symbol"	,"Variant_Classification","Chromosome",	"Start_position",	"End_position",
                 "nProbes"	,"Modal_HSCN_1",	"Modal_HSCN_2",
                 "Modal_Total_CN"	,"purity"	,"ploidy")

pyclone.input<-data.frame(mutation_id = paste(paste(snv.int.cnv[,8], ":",snv.int.cnv[,1],"-", 
                                                    snv.int.cnv[,2],"_",snv.int.cnv[,3], sep="")), 
                          ref_counts=snv.int.cnv$Ref_freq, var_counts = snv.int.cnv$Variant_freq, 
                          normal_cn = rep(2,nrow(snv.int.cnv)),
                          minor_cn = snv.int.cnv$Modal_HSCN_1, major_cn = snv.int.cnv$Modal_HSCN_2,
                          total_cn = (snv.int.cnv$Modal_HSCN_1 + snv.int.cnv$Modal_HSCN_2),
                          purity = snv.int.cnv$purity, ploidy = snv.int.cnv$ploidy,
                          var_freq = (snv.int.cnv$Variant_freq / (snv.int.cnv$Variant_freq + snv.int.cnv$Ref_freq)),                        
                          depth = (snv.int.cnv$Variant_freq + snv.int.cnv$Ref_freq),
                          norm_depth = (snv.int.cnv$Variant_freq + snv.int.cnv$Ref_freq)*snv.int.cnv$purity,
                          stringsAsFactors = F,
                          Chr=snv.int.cnv$Chr, Start=snv.int.cnv$Start, End= snv.int.cnv$End)


#(Optional) Apply filters
pyclone.input.filtered<-pyclone.input %>% filter(minor_cn == 1) %>% filter(major_cn ==1) %>% filter(norm_depth > 20) %>% filter(var_freq > 0.1)

if (nrow(pyclone.input.filtered) > 0) {

  write.table(pyclone.input.filtered, outFile,sep="\t", row.names = F,
              quote = F)
} else {
  cat ("File without mutations after filtering\n")

}

