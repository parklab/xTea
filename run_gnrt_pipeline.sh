#!/bin/bash

 SAMPLE_ID=$1
 BAMS=$2
 X10_BAM=$3
 WFOLDER=$4
 OUT_SCRTP=run_jobs.sh
 TIME=2-5:00
 REP_LIB=$5
 python gnrt_pipeline.py -i ${SAMPLE_ID} -b ${BAMS} -x ${X10_BAM} -p ${WFOLDER} -o ${OUT_SCRTP} -q park -n 16 -m 70 -t ${TIME} \
 -l ${REP_LIB} --nclip 5 --cr 3 --nd 8 --nfclip 3 --nfdisc 6 --flklen 3000 -f 19
