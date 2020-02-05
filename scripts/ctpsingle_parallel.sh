#!/bin/bash

#SBATCH -J ctpsingle_parallel
#SBATCH -o %x-%J.out
#SBATCH --time=100:00:00
#SBATCH -n 6
#SBATCH -p shared --qos=shared
#SBATCH -c 8
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT_80,END
#SBATCH --mail-user=avicens@uvigo.es


CTPSDIR="/home/uvi/be/avs/store/MyRlibs/ctpsingle"
WORKDIR="/home/uvi/be/avs/store/dNdS_clones"
DATADIR="$WORKDIR/data"
SCRIPTDIR="$WORKDIR/scripts"
CTPSINPUTDIR="$DATADIR/ctpsingle/ctpsingle_input"
CTPSOUTPUTDIR="$DATADIR/ctpsingle/ctpsingle_output"

module load parallel
LOGDIR=$WORKDIR/jobs/logs
mkdir -p $LOGDIR

SAMPLESFILE="${WORKDIR}/ctpsingle_input_files.txt"
mkdir -p ${CTPSOUTPUTDIR}

#shared partitions
MEMPERCORE=$(eval $(scontrol show partition $SLURM_JOB_PARTITION -o);echo $DefMemPerCPU)
if [ -z "$MEMPERCORE" ]
  then
  #exclusive partitions
  MEMPERCORE=$(( $(sinfo -e -p $SLURM_JOB_PARTITION -o "%m/%c" -h) ))
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK=1}

SRUN="srun -o $LOGDIR/%x-%J-%t.out -N1 -n1 --mem=$(( $MEMPERCORE*$OMP_NUM_THREADS )) -c $OMP_NUM_THREADS --cpu_bind=none"

parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog $LOGDIR/runtask.log --resume-failed"

$parallel --colsep '\t' "$SRUN bash ${SCRIPTDIR}/ctpsingle.sh {1} {#}" :::: ${SAMPLESFILE}

