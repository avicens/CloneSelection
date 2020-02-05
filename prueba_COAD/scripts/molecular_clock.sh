#!/bin/bash
#SBATCH --mail-type=BEGIN,FAIL,TIME_LIMIT_80,TIME_LIMIT,END
#SBATCH --mail-user=avicens@uvigo.es
#SBATCH -J molecular_clock
#SBATCH -o %x-%A_%a.out
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -p shared
#SBATCH --qos shared
#SBATCH --mem 4G

#Load modules
module load gcccore/6.4.0 paml/4.9i

#Set directories
WORKDIR="/home/uvi/be/avs/store/dNdS_clones/prueba_COAD"
DATADIR="$WORKDIR/data"
SEQSDIR="${TRACEDIR}/seqs"
RAXMLDIR="${TRACEDIR}/raxml_trees"
CLOCKDIR="${TRACEDIR}/molecular_clock"

#Set file of samples for molecular clock test (those for which an alignment and a tree were built)
SAMPLES="$RAXMLDIR/raxml_trees_samples.txt"

#Set template Codeml file
CODEMLFILE="$CLOCKDIR/codeml.ctl"

#SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLES}`
SAMPLE=`sed  "1q;d" ${SAMPLES}`

echo "Analyzing sample ${SAMPLE}"
echo "Creating Codeml control file"
SEQFILE="$SEQSDIR/${SAMPLE}_clone_seqs.fas"
ROOTREE="$RAXMLDIR/rooted/${SAMPLE}_rooted.raxml.bestTree"
UNROOTREE="$RAXMLDIR/unrooted/${SAMPLE}_unrooted.raxml.bestTree"
CODEMLFILE="$CLOCKDIR/codeml.ctl"

#Create directory for molecular clock analysis
mkdir -p ${CLOCKDIR}/${SAMPLE}

#1.Setting models
#1.1. Model that does not assume molecular clock
mkdir -p ${CLOCKDIR}/${SAMPLE}/no_clock
 
##Copy input files to Molecular Clock directory
###Alignment file
cp ${SEQFILE} ${CLOCKDIR}/${SAMPLE}/no_clock

###Unrooted tree topology
sed 's/:0.[0-9]*//g' ${UNROOTREE} > ${CLOCKDIR}/${SAMPLE}/no_clock/${SAMPLE}_unrooted.tre

###Control file
echo "Editing control file"

sed -e "s/seqfile = .*$/seqfile = ${SAMPLE}_clone_seqs.fas/g" ${CODEMLFILE} | 
sed -e "s/treefile = .*$/treefile = ${SAMPLE}_unrooted.tre/"  |
sed -e "s/outfile = .*$/outfile = ${SAMPLE}_noclock.out/" |
sed -e "s/clock = [0-3]/clock = 0/" > ${CLOCKDIR}/${SAMPLE}/no_clock/${SAMPLE}_noclock.ctl

##1.2.Setting model that assumes constant molecular clock
mkdir -p ${CLOCKDIR}/${SAMPLE}/global_clock
 
##Copy input files to MC directory
###Alignment file
cp ${SEQFILE} ${CLOCKDIR}/${SAMPLE}/global_clock

###Rooted tree topology
sed 's/:0.[0-9]*//g' ${ROOTREE} > ${CLOCKDIR}/${SAMPLE}/global_clock/${SAMPLE}_rooted.tre

###Control file
echo "Editing control file"

sed -e "s/seqfile = .*$/seqfile = ${SAMPLE}_clone_seqs.fas/g" ${CODEMLFILE} | 
sed -e "s/treefile = .*$/treefile = ${SAMPLE}_rooted.tre/"  |
sed -e "s/outfile = .*$/outfile = ${SAMPLE}_global_clock.out/" |
sed -e "s/clock = [0-3]/clock = 1/" > ${CLOCKDIR}/${SAMPLE}/global_clock/${SAMPLE}_global_clock.ctl


#2.Running clock models
echo "Running Codeml"

##2.1. No-clock model
echo "Running no-clock model"
cd ${CLOCKDIR}/${SAMPLE}/no_clock
codeml ${SAMPLE}_noclock.ctl


##2.2 Global-clock model
echo "Running Global-clock model"
cd ${CLOCKDIR}/${SAMPLE}/global_clock
codeml ${SAMPLE}_global_clock.ctl

echo "Work finished"
