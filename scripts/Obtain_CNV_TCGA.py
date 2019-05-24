import pandas as pd

def ObtainCNV(CNA_TCGA,purity,outpath_CNV1,outpath_CNV2):

    CNA_TCGA[['Chromosome', 'Start','End','Modal_Total_CN']] = CNA_TCGA[['Chromosome','Start','End','Modal_Total_CN']].fillna(value=0)
    CNA_TCGA['Chromosome'] = CNA_TCGA['Chromosome'].astype(int)
    CNA_TCGA['Start'] = CNA_TCGA['Start'].astype(int)
    CNA_TCGA['End'] = CNA_TCGA['End'].astype(int)

    CNA_TCGA2 = CNA_TCGA[(CNA_TCGA['Chromosome'] != 0) | (CNA_TCGA['End'] != 0)
                        | (CNA_TCGA['Start'] != 0)]

    CNA_TCGA3 = CNA_TCGA2[CNA_TCGA2.columns[0]]
    CNAV = CNA_TCGA3.str.rsplit('-', expand=True,n=1)
    CNAV.columns = ['Sample','None']
    SampleID = CNAV[['Sample']]
    CNA_TCGA2 = CNA_TCGA2.drop(['Sample'],axis=1)

    SegMatCopy = pd.concat([SampleID, CNA_TCGA2],axis=1)
    SegMatCopy = SegMatCopy.rename(columns = {'Sample':'SampleID','Chromosome':'Chr','Num_Probes':'nProbes'})

    Purity2 = purity[purity.columns[1]]
    PPurity = Purity2.str.rsplit('-', expand=True,n=4)
    PPurity.columns = ['SampleID','None','None','None','None']
    SampleIDPP = PPurity[['SampleID']]
    PurityCon = purity[['purity','ploidy']]
    PurityValues = pd.concat([SampleIDPP, PurityCon],axis=1)

    TotalSegMatCopyTable = pd.merge(SegMatCopy, PurityValues, on='SampleID')

    TotalSegMatCopyTable = TotalSegMatCopyTable[TotalSegMatCopyTable['Chr'] != 23]
    TotalSegMatCopyTable = TotalSegMatCopyTable.drop(['Length','Subclonal_HSCN_a1','Subclonal_HSCN_a2','Ccf_ci95_low_a1','Ccf_ci95_high_a1',
                                                     'Cancer_cell_frac_a2','Ccf_ci95_low_a2','Ccf_ci95_high_a2',
                                                     'LOH','Homozygous_deletion','solution','Cancer_cell_frac_a1'],axis=1)

    for my_index, my_df in TotalSegMatCopyTable.groupby('SampleID'):
        my_df = my_df.drop(['SampleID'],axis=1)
        my_df.to_csv(outpath_CNV1 + my_index + outpath_CNV2,sep='\t', index=None,compression='gzip')