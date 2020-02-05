#!/bin/bash
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT_80,TIME_LIMIT,END
#SBATCH --mail-user=avicens@uvigo.es
#SBATCH -J dndscv_bulk 
#SBATCH -o %x-%A_%a.out
#SBATCH --time=00:10:00
#SBATCH -n 1
#SBATCH -p shared --qos=shared
#SBATCH -c 2
#SBATCH --mem-per-cpu=10GB

#source readConfig.sh $1

WORKDIR="/home/uvi/be/avs/store/dNdS_clones/prueba_COAD"
DATADIR="$WORKDIR/data"
SCRIPTDIR="$WORKDIR/scripts"
SAMPLEDIR="$WORKDIR/samples"
DNDSCVDIR="$DATADIR/dndscv"

SAMPLES="$SAMPLEDIR/COAD_samples.txt"

module load gcc/6.4.0 R/3.5.3

SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLES}`

mkdir -p $DNDSCVDIR/${SAMPLE}
echo "Running dNdScv script for  ${SAMPLE}"
SNV="$DATADIR/snv/${SAMPLE}_SNV.tsv.gz"
LOCI="$DATADIR/pyclone/loci/${SAMPLE}_loci.txt"
INFILE="$DNDSCVDIR/${SAMPLE}/${SAMPLE}_dndscv_input.txt"
OUTFILE="$DNDSCVDIR/${SAMPLE}/${SAMPLE}_bulk_dnds.txt"

Rscript $SCRIPTDIR/dndscv_samples.R $SNV $SAMPLE $INFILE $OUTFILE


