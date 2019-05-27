import pandas as pd
from Obtain_SNV_TCGA import ObtainSNV

input_path = '~/store/dNdS_clones/metadata/mc3.v0.2.8.PUBLIC.maf.gz'
SNV_TCGA = pd.read_csv(input_path,compression='gzip',sep='\t',low_memory=False)
out_snv1 = '~/store/dNdS_clones/data/snv/'
out_snv2 = '_SNV.tsv.gz'

ObtainSNV(SNV_TCGA,out_snv1,out_snv2)
