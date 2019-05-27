import pandas as pd

def countStr(s):
    return len(s)

def ObtainSNV(SNV_TCGA_NoF,outpath_snv1,outpath_snv2):

    #Take filtered ones
    SNV_TCGA = SNV_TCGA_NoF[SNV_TCGA_NoF['FILTER'] == 'PASS']

    #reduce table and modify column names
    FinalTableSNV = SNV_TCGA[['Matched_Norm_Sample_Barcode','Chromosome','Start_Position','End_Position','Reference_Allele','Allele',
                           't_alt_count', 't_ref_count','Hugo_Symbol','Variant_Classification']]
    FinalTableSNV = FinalTableSNV.rename(columns = {'Matched_Norm_Sample_Barcode':'Sample','Chromosome':'Chr',
                                              'Reference_Allele':'Reference','Allele':'Alternate','t_alt_count':'Variant_freq',
                                             't_ref_count':'Ref_freq'})

    #Arrange Sample Names
    ArrangePatient = FinalTableSNV[FinalTableSNV.columns[0]]
    ArrangePatient2= ArrangePatient.str.rsplit('-', expand=True,n=4)
    ArrangePatient2.columns = ['Sample','None','None','None','None']
    SampleID = ArrangePatient2[['Sample']]
    FinalTableSNV = FinalTableSNV.drop(['Sample'],axis=1)
    MutationTableSNV = pd.concat([SampleID, FinalTableSNV],axis=1)

    #Delete Indels and select only SNV
    MutationTableSNV['lenREF'] = MutationTableSNV['Reference'].astype(str).apply(countStr)
    MutationTableSNV['lenALT'] = MutationTableSNV['Alternate'].astype(str).apply(countStr)

    MutationTableSNV = MutationTableSNV[(MutationTableSNV['lenREF'] == 1)&(MutationTableSNV['lenALT'] == 1)]
    MutationTableSNV = MutationTableSNV[(MutationTableSNV['Alternate'] != '-')&(MutationTableSNV['Reference'] != '-')]
    MutationTableSNV = MutationTableSNV[(MutationTableSNV['Chr'] != 'X')&(MutationTableSNV['Chr'] != 'Y')] #No sexual chromosomes
    MutationTableSNV = MutationTableSNV.drop(['lenREF','lenALT'],axis=1)

    for my_index, my_df in MutationTableSNV.groupby('Sample'):
        my_df = my_df.drop(['Sample'],axis=1)
        my_df.to_csv(outpath_snv1 + my_index + outpath_snv2 ,sep='\t', index=None,compression='gzip')