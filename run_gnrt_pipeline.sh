#!/bin/bash

SAMPLE_ID=$1
BAMS=$2
X10_BAM=$3
WFOLDER=$4
OUT_SCRTP=run_jobs.sh
REP_LIB=$5
REF=$6
XTEA=$7

python gnrt_pipeline_cloud.py -i ${SAMPLE_ID} -b ${BAMS} -x ${X10_BAM} -p ${WFOLDER} -o ${OUT_SCRTP} -n 8 \
-l ${REP_LIB} -r ${REF} -x ${XTEA} --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 5 --flklen 3000 -f 19

#1. tar rep lib folder to one file
#2. tar reference genome to one file
#3. change to bam file
#python ./../xTEA/gnrt_pipeline_cloud.py -b input.bam -p /home/ec2-user/results2/ -o run_jobs.sh -n 16 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -x /home/ec2-user/xTEA/ --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 5 --flklen 3000 -f 19

####