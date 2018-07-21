#!/bin/bash

SAMPLE_ID=$1
BAMS=$2
X10_BAM=$3
WFOLDER=$4
OUT_SCRTP=run_jobs.sh
REP_LIB=$5
REF=$6
XTEA=$7

python gnrt_pipeline_cloud.py -i ${SAMPLE_ID} -b ${BAMS} -x ${X10_BAM} -p ${WFOLDER} -o ${OUT_SCRTP} -n 16 \
-l ${REP_LIB} -r ${REF} -x ${XTEA} --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 4 --flklen 3000 -f 19
