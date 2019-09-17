#!/bin/bash
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT_80,TIME_LIMIT,END
#SBATCH --mail-user=avicens@uvigo.es
#SBATCH -J raxml_arraye 
#SBATCH -o %x-%A_%a.out
#SBATCH --time=01:00:00
#SBATCH -n 4
#SBATCH -p shared --qos=shared
#SBATCH -c 2


WORKDIR="/home/uvi/be/avs/store/dNdS_clones/prueba_COAD"
DATADIR="$WORKDIR/data"
SCRIPTDIR="$WORKDIR/scripts"
SAMPLES="$DATADIR/seqs_2/clone_seqs_samples.txt"

module load gcc/6.4.0 openmpi/2.1.1 raxml-ng/0.8.1

SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLES}`
#SAMPLE=`sed "1q;d" ${SAMPLES}`

echo "Analyzing ${SAMPLE}"

CLONESEQS="${DATADIR}/seqs_2/${SAMPLE}_clone_seqs.fas"
echo "Building phylogenetic tree (without outgroup)"
mkdir -p $DATADIR/raxml_trees/unrooted/
UNROOTTREEDIR="$DATADIR/raxml_trees/unrooted/${SAMPLE}_unrooted"
raxml-ng-mpi --msa $CLONESEQS --model GTR+G --precision 12 --blmin 0.000000001 --redo --prefix $UNROOTTREEDIR

echo "Building phylogenetic tree (assigning healthy as outgroup)"
mkdir -p $DATADIR/raxml_trees/rooted/
ROOTTREEDIR="$DATADIR/raxml_trees/rooted/${SAMPLE}_rooted"
raxml-ng-mpi --msa $CLONESEQS --model GTR+G --precision 12 --blmin 0.000000001 --outgroup ref --redo --prefix  $ROOTTREEDIR

echo "Work finished"
