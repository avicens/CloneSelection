import pandas as pd
from Obtain_CNV_TCGA import ObtainCNV
from Obtain_SNV_TCGA import ObtainSNV

#CNV file
input_path = '/home/uvi/be/avs/store/dNdS_clones/metadata/mc3.v0.2.8.PUBLIC.maf.gz'
CNA_TCGA = pd.read_csv(input_path,sep='\t')
out_cnv1 = '/home/uvi/be/avs/store/dNdS_clones/data/cnv/'
out_cnv2 = '_CNV.tsv.gz'

#Purity files
input_path = '/home/uvi/be/avs/store/dNdS_clones/metadata/TCGA-mastercalls.abs_tables_JSedit.fixed.txt'
purity = pd.read_csv(input_path,sep='\t')

#SNV files
input_path = '/home/uvi/be/avs/store/dNdS_clones/metadata/TCGA-mastercalls.abs_segtabs.fixed'
SNV_TCGA = pd.read_csv(input_path,compression='gzip',sep='\t',low_memory=False)
out_snv1 = '~/Dropbox/proyectos/dNdS_clones/data/snv/'
out_snv2 = '_SNV.tsv.gz'

ObtainCNV(CNA_TCGA,purity,out_cnv1,out_cnv2)
ObtainSNV(SNV_TCGA,out_snv1,out_snv2)
