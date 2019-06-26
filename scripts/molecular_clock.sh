#!/usr/bin/bash

WORKDIR="/Users/avicens/Dropbox/proyectos/dNdS_clones"
DATADIR="$WORKDIR/data"
CLOCKDIR="$DATADIR/molecular_clock"
BASEMLFILE="$WORKDIR/baseml.ctl"

SAMPLE="TCGA-02-0003"

SEQFILE="$DATADIR/seqs/${SAMPLE}_clone_seqs.fas"
ROOTREE="$DATADIR/raxml_trees/rooted/$SAMPLE_rooted.raxml.bestTree"
UNROOTREE="$DATADIR/raxml_trees/unrooted/$SAMPLE.raxml.bestTree"

#Create directory for molecular clock analysis
mkdir -p ${CLOCKDIR}/${SAMPLE}/

#1.Setting models
#1.1. Model that does not assume molecular clock
mkdir -p ${CLOCKDIR}/${SAMPLE}/NOCLOCK
 
##Copy input files to MC directory
###Alignment file
cp ${SEQFILE} ${CLOCKDIR}/${SAMPLE}/NOCLOCK

###Unrooted tree topology
sed 's/:0.[0-9]*//g' ${UNROOTREE} > ${CLOCKDIR}/${SAMPLE}/NOCLOCK/${SAMPLE}_unrooted.tre

###Control file
echo "Editing control file"

sed -e "s/seqfile = .*$/seqfile = ${SAMPLE}_clone_seqs.fas/g" ${BASEMLFILE} | 
sed -e "s/treefile = .*$/treefile = ${SAMPLE}_unrooted.tre/"  |
sed -e "s/outfile = .*$/outfile = ${SAMPLE}_noclock.out/" |
sed -e "s/clock = [0-3]/clock = 0/" > ${CLOCKDIR}/${SAMPLE}/NOCLOCK/${SAMPLE}_noclock.ctl

##1.2.Setting model that assumes constant molecular clock
mkdir -p ${CLOCKDIR}/${SAMPLE}/GLOBALCLOCK
 
##Copy input files to MC directory
###Alignment file
cp ${SEQFILE} ${CLOCKDIR}/${SAMPLE}/GLOBALCLOCK

###Rooted tree topology
sed 's/:0.[0-9]*//g' ${ROOTREE} > ${CLOCKDIR}/${SAMPLE}/GLOBALCLOCK/${SAMPLE}_rooted.tre


###Control file
echo "Editing control file"

sed -e "s/seqfile = .*$/seqfile = ${SAMPLE}_clone_seqs.fas/g" ${BASEMLFILE} | 
sed -e "s/treefile = .*$/treefile = ${SAMPLE}_rooted.tre/"  |
sed -e "s/outfile = .*$/outfile = ${SAMPLE}_global_clock.out/" |
sed -e "s/clock = [0-3]/clock = 1/" > ${CLOCKDIR}/${SAMPLE}/GLOBALCLOCK/${SAMPLE}_global_clock.ctl

#2.Run baseml
##2.1. No-clock model
cd ${CLOCKDIR}/${SAMPLE}/NOCLOCK
baseml ${SAMPLE}_noclock.ctl
echo "Work finished"

#2.2 Global-clock model
cd ${CLOCKDIR}/${SAMPLE}/GLOBALCLOCK
baseml ${SAMPLE}_global_clock.ctl
