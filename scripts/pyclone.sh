#!/bin/bash

WORKDIR="/home/uvi/be/avs/store/dNdS_clones"
SCRIPTDIR="$WORKDIR/scripts"
TRACEDIR="$WORKDIR/trace"
PYCLONEDIR="$TRACEDIR/pyclone"

SAMPLE=$1
PURITY=$2

#Load miniconda module
module load miniconda3

#Activate environment with python 2.7
source activate avicens2.7

echo "Running pyclone for sample ${SAMPLE}"
PyClone run_analysis_pipeline --in_files $PYCLONEDIR/pyclone_input/${SAMPLE}_pyclone.txt  --working_dir $PYCLONEDIR/pyclone_output/${SAMPLE} --tumour_contents ${PURITY}  --samples ${SAMPLE} --density pyclone_binomial --prior major_copy_number --min_cluster_size 2
