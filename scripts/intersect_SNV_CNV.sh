#!/bin/bash
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=avicens@uvigo.es
#SBATCH -J intersect_snv_cnv -o %x-%J.out
#SBATCH --time=01:30:00
#SBATCH -n 1
#SBATCH -p shared --qos=shared
#SBATCH -c 2

WORKDIR="/home/uvi/be/avs/store/dNdS_clones"
DATADIR="$WORKDIR/data"
SAMPLEDIR="$WORKDIR/samples"
SAMPLES="$SAMPLEDIR/intersect_samples.txt"

mkdir -p $DATADIR/intersect

module load gcccore/6.4.0 bedtools/2.27.1
while read -r SAMPLE
do
SNV="$DATADIR/snv/${SAMPLE}_SNV.tsv.gz"
CNV="$DATADIR/cnv/${SAMPLE}_CNV.tsv.gz"
INTERSFILE="$DATADIR/intersect/${SAMPLE}_intersect_SNV_CNV.txt"

echo "creating intersect for ${SAMPLE}"

bedtools intersect -a $SNV -b $CNV -wa -wb > $INTERSFILE
done <  ${SAMPLES}
