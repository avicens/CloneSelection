#!/bin/bash
#SBATCH -J pyclone_parallel -o %x-%J.out
#SBATCH -t 100:00:00
#SBATCH -n 6
#SBATCH -p shared --qos=shared
#SBATCH -c 8
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT_50,TIME_LIMIT_80,END 
#SBATCH --mail-user=avicens@uvigo.es

WORKDIR="/home/uvi/be/avs/store/dNdS_clones"
SCRIPTDIR="$WORKDIR/scripts"
TRACEDIR="$WORKDIR/trace"
PYCLONEDIR="$TRACEDIR/pyclone"
CONFIGFILE="$WORKDIR/pending_samples_pyclone_21-oct-2019.txt"

module load parallel
LOGDIR=$WORKDIR/jobs/logs
mkdir -p $LOGDIR

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

$parallel --header : --colsep '\t' "$SRUN bash $SCRIPTDIR/pyclone.sh  {7} {1} {#}" :::: ${CONFIGFILE}
