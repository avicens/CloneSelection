#!/usr/bin/env python

import os
import sys
import pandas as pd

ctpdir = sys.argv[1]
num_clusters_file=sys.argv[2]
cluster_size_file=sys.argv[3]

#ctpdir = os.path.join(cwd,'data','ctpsingle','ctpsingle_output2')
dirs = os.listdir(ctpdir)
samples = [s for s in dirs if "TCGA" in s]

#Open files
#num_clusters_file=os.path.join(ctpdir,'num_clusters_per_sample.tsv')
fh=open(num_clusters_file,'w')
fh.write("sample\tclusters\n")

#cluster_size_file = os.path.join(ctpdir,'clusters_size_per_sample.tsv')
fi=open(cluster_size_file,'w')
fi.write("sample\tcluster\tsize\n")

for sample in samples:
    input_path=os.path.join(ctpdir,sample,(sample + "_ctpsingle_cluster_assignments.txt"))
    if os.path.exists(input_path):
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
    else:
        continue

fh.close()
fi.close()

print("Work finished, good bye!")

