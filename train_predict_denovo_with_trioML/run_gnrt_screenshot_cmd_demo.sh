xTEA=/n/projects/XTEA/code/xTEA/
WFOLDER=/n/denovo_MEI/gnrt_training_data/Ewing_Sarcoma/wfolder/
REF=/n/projects/data/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
TYPE=L1
TP=/n/denovo_MEI/gnrt_training_data/TP_manual_inspection.txt

python ${xTEA}"prep_denovo_training_data_from_trio.py" --pa pa_folder_list.txt --ma ma_folder_list.txt --proband proband_folder_list.txt --bam Ewing_Sarcoma_sample_bams_all.txt -p ${WFOLDER} -r ${REF} -e 450 --type ${TYPE} -o "generate_positive_training_images_"${TYPE}".sh" --tp ${TP}
