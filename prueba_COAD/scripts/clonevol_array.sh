#!/bin/bash
##SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT_80,TIME_LIMIT,END
#SBATCH --mail-user=avicens@uvigo.es
#SBATCH -J get_clone_seqs 
#SBATCH -o %x-%A_%a.out
#SBATCH --time=00:30:00
#SBATCH -n 4
#SBATCH -p shared --qos=shared
#SBATCH -c 2

WORKDIR="/home/uvi/be/avs/store/dNdS_clones/prueba_COAD"
DATADIR="$WORKDIR/data"
SCRIPTDIR="$WORKDIR/scripts"
LOCIDIR="$DATADIR/loci"
CLONEVOLDIR="$DATADIR/clonevol"
SEQDIR="$DATADIR/seqs"

SAMPLES="$DATADIR/pyclone_output/finished_samples_pyclone.sh"

#Create folders where saving the files for the specific test
mkdir -p $CLONEVOLDIR

#Load modules
module load gcc/6.4.0 R/3.5.3

#SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLES}`
SAMPLE=`sed "1q;d" ${SAMPLES}`

echo "Analyzing ${SAMPLE}"
echo "Getting clone evolution tree and processing loci" 
LOCI="$LOCIDIR/${SAMPLE}_loci.txt"
CLONETREE="$CLONEVOLDIR/${SAMPLE}_clone_tree.txt"

Rscript $SCRIPTDIR/clone_tree_inference_clonevol.R $LOCI $CLONETREE $SCRIPTDIR

echo "Work finished"

