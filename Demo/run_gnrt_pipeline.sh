#!/bin/bash

SAMPLE_ID=sample_id.txt
BAMS=sample_bam.txt
X10_BAM=null
WFOLDER=[abosolute-path-of-working-folder]/
OUT_SCRTP=submit_jobs.sh
TIME=0-05:00
REP_LIB=[abosolute-path-of-repeat-library-folder]/rep_lib_annotation/
REF=[abosolute-path-of-reference-path]/refdata-b37-2.1.0/fasta/genome.fa
GENE=[abosolute-path-of-gene-annotation-file-path]/gene_annotation/gencode.v28lift37.annotation.gff3
XTEA=[xTea-folder-path]/xTea/
BLK_LIST=[abosolute-path-of-repeat-library-folder]/rep_lib_annotation/rep_lib_annotation/blacklist/hg19/centromere.bed

#
python ${XTEA}"gnrt_pipeline_local.py" -i ${SAMPLE_ID} -b ${BAMS} -x ${X10_BAM} -p ${WFOLDER} -o ${OUT_SCRTP} -q short -n 8 -m 16 -t ${TIME} \
-l ${REP_LIB} -r ${REF} -g ${GENE} --xtea ${XTEA}  --nclip 4 --cr 2 --nd 5 --nfclip 4 --nfdisc 5 --flklen 3000 -f 5907  -y 7  --blacklist ${BLK_LIST}
