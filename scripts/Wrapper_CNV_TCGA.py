import pandas as pd
from Obtain_CNV_TCGA import ObtainCNV

#CNV file
input_path = '/home/uvi/be/avs/store/dNdS_clones/metadata/TCGA_mastercalls.abs_segtabs.fixed.txt'
CNA_TCGA = pd.read_csv(input_path,sep='\t')
out_cnv1 = '/home/uvi/be/avs/store/dNdS_clones/metadata/cnv/'
out_cnv2 = '_CNV.tsv.gz'

#Purity files
input_path = '/home/uvi/be/avs/store/dNdS_clones/metadata/TCGA_mastercalls.abs_tables_JSedit.fixed.txt'
purity = pd.read_csv(input_path,sep='\t')

ObtainCNV(CNA_TCGA,purity,out_cnv1,out_cnv2)
