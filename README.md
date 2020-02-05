Clonal selection
=======
# Pan-cancer analysis of selection
In this project, I conducted a Pan-Cancer analysis to examine selection during clonal evolution. 

## Repository content
The repository contains the scripts written to perform the computational workflow of this study, as well as some example files required as input in some steps.

The scripts were written to upload jobs to a SLURM-based scheduling system (CESGA FinisTerrae II). The scripts whose filename end with *array.sh* call to launch multiple jobs with the same parameters. The scripts whose filename end with *parallel.sh* run  parrallel jobs as a single task using GNU parallel.

## Cancer genomic data
I used meta-cancer data of [PanCanAtlas](https://gdc.cancer.gov/about-data/publications/pancanatlas). This project assemble information from 10,824 tumor samples across 33 cancer types from the Cancer Genome Atlas (TCGA) (Study Accession: phs000178.v10.p8).

Somatic mutations for each tumor sample were obtained from the metafile **mc3.v0.2.8.PUBLIC.maf.gz** available in PanCanAtlas. This metafile was generated by the Multi-Center Mutation Calling in Multiple Cancers project [MC3](https://www.sciencedirect.com/science/article/pii/S2405471218300966?via%3Dihub). For extracting mutations per tumor sample, I used the wrapper script **Wrapper_SNV_TCGA.py**, which calls the script **Obtain_SNV_TCGA.py**. 

I obtained copy number state, ploidy and purity information for 9786 tumor samples from a metafile (**TCGA_mastercalls.abs_segtabs.fixed.txt**) computed by [ABSOLUTE](https://software.broadinstitute.org/cancer/cga/absolute) tool. For extracting copy number state per tumor sample, I used the script **Wrapper_CNV_TCGA.py**, which calls the script **Obtain_CNV_TCGA.py**.

## Clonal deconvolution

I used two softwares to infer clones based on allele frequency, copy number state and ploidy. *PyClone* and *CTPsingle*.

To ensure we use high-quality variant data for clone inference, I retained the variants with variant frequency > 0.1, normalized depth > 20 and Copy Number = 2 (i.e. minor CN = 1 and, major CN =1) (this is because PyClone is not able to incorporate subclonal SNVs).

The input for running these two programs was obtained with the R scripts **get_pyclone_input.R** and **get_ctpsingle_input.R**, and were called in loop with **get_pyclone_input_array.R** and **get_ctpsingle_input_array.R**

## Get clonal sequences
I retrieved sequences containing the mutated codons for each clone using the script **get_clone_sequences.R**. This script load a function based on the [dNdScv](https://github.com/im3sanger/dndscv/tree/master/R) library to extract the codons containing each mutation, taking the human genome hg19 as reference.
In CESGA, the script was executed for every samples in an array job (**get_clone_seqs_array.sh**)

## Reconstruct phylogenetic trees of clones
I reconstruct phylogentic trees using the clonal sequences (plus reference) with [RAxML-ng](https://github.com/amkozlov/raxml-ng). Both rooted as unrooted trees were buil. Rooted trees were used for the global/ (strict) molecular clock model, and unrooted for the relaxed clock model and dN/dS estimation.

In CESGA, the script was executed for every samples in an array job (**raxml_array.sh**).

## Estimation of selective pressures (dN/dS)
I estimated dN/dS in tumor samples using **codeml**. I estimated global dN/dS for tumor sample with the M0 model, and I estimated a dN/dS per branch using the *fb* model. PAML analysis were run with the [ETE toolkit](http://etetoolkit.org/). 

To assess the hypothesis of dN/dS being constant or variable across branches, I applied a Likelihood Ratio TEST (LRT) comparing the
In CESGA, the script was executed for every samples in an array job (**ete_codmel_array.sh**).




>>>>>>> 342596516e4213931c7a790152b8c37e09008ed1
