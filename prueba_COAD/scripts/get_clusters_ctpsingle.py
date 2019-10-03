#!/usr/bin/env python

import os
import pandas as pd

cwd=os.getcwd()
dirs = os.listdir(os.path.join(cwd,'data','ctpsingle'))
samples = [s for s in dirs if "TCGA" in s]

#Open files
num_clusters_file=os.path.join(cwd,'data','ctpsingle','num_clusters_per_sample.tsv')
fh=open(num_clusters_file,'w')
fh.write("sample\tclusters\n")

cluster_size_file = os.path.join(cwd,'data','ctpsingle','clusters_size_per_sample.tsv')
fi=open(cluster_size_file,'w')
fi.write("sample\tcluster\tsize\n")

for sample in samples[0:2]:
    
    input_path=os.path.join(cwd,'data','ctpsingle',sample,(sample + "_ctpsingle_cluster_assignments.txt"))
    clusters_file=pd.read_csv(input_path, sep='\t')
    
    #Add lines with number of clusters per sample
    num_clusters=clusters_file['mostLikely'].nunique()
    print(sample+"\t"+str(num_clusters), file=fh)

    #Add lines with sizes of clusters from each sample
    cluster_size=clusters_file['mostLikely'].value_counts()
    cluster_size_df=cluster_size.to_frame()
    cluster_size_df.columns = ['size']
    cluster_size_df.insert(0, "cluster", cluster_size.index) #Add column with cluster id (index)
    cluster_size_df.insert(0, "sample", sample) #Add colum with sample code
    cluster_size_df.to_csv(fi,sep='\t',header=False, index =False)
fh.close()
fi.close()

