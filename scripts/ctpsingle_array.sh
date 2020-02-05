#!/bin/bash
# SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=avicens@uvigo.es
#SBATCH -J run_ctpsingle_array
#SBATCH -o %x-%A_%a.out
#SBATCH --time=100:00:00
#SBATCH -n 1
#SBATCH -p shared --qos=shared
#SBATCH -c 4

CTPSDIR="/home/uvi/be/avs/store/MyRlibs/ctpsingle"
WORKDIR="/home/uvi/be/avs/store/dNdS_clones"
DATADIR="$WORKDIR/data"
SCRIPTDIR="$WORKDIR/scripts"
CTPSINPUTDIR="$DATADIR/ctpsingle/ctpsingle_input"
CTPSOUTPUTDIR="$DATADIR/ctpsingle/ctpsingle_output"

module load gcc/6.4.0 R/3.5.3

SAMPLES=$(ls ${CTPSINPUTDIR} | egrep -o "TCGA-[0-9A-Z]{2}-[A-Z0-9]{4}")

#SAMPLE=`echo "${SAMPLES}" | sed -n ${SLURM_ARRAY_TASK_ID}p`
SAMPLE=`echo "${SAMPLES}" | sed -n 401p`
echo "RUNNING CTPSINGLE FOR ${SAMPLE}"
mkdir -p $CTPSOUTPUTDIR/${SAMPLE}
Rscript $CTPSDIR/CTPsingle2.R -f $CTPSINPUTDIR/${SAMPLE}_ctpsingle_input.txt -o $CTPSOUTPUTDIR/${SAMPLE}/${SAMPLE}_ctpsingle -m $CTPSDIR/GammaAdjMatrices

#done <<< "${SAMPLES}"
echo "Work finished"
