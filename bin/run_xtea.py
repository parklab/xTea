#!/usr/bin/env python

##05/08/2024
##@@author: Corinne Sexton, DBMI, Harvard Medical School
##@@contact: corinne_sexton@hms.harvard.edu

from pathlib import Path
import configargparse
from xtea.locate_insertions import get_clip_sites,get_disc_sites,filter_csn,get_transduction,get_sibling,filter_sites_post,annotate_genes,call_genotypes, generate_VCF


## run_xtea -c config.toml -i bam_list (or file) -o output_dir
##          --repeat_type <list> --run_type germline (default)

# CURRENT SETUP:
# time python ${XTEA_PATH}"x_TEA_main.py" -C -i ${BAM_LIST} --lc 3 --rc 3 --cr 1  
# -r ${L1_COPY_WITH_FLANK}  -a ${ANNOTATION} --cns ${L1_CNS} --ref ${REF} -p ${TMP} 
# -o ${PREFIX}"candidate_list_from_clip.txt"  -n 8 
# --cp /n/data1/hms/dbmi/park/DATA/Kamihara_DFCI_WGS/thyroid_cancer/analysis_terra_2023/xtea/results/TP003/pub_clip/

# time python ${XTEA_PATH}"x_TEA_main.py"  -D -i ${PREFIX}"candidate_list_from_clip.txt" 
# --nd 5 --ref ${REF} -a ${ANNOTATION} -b ${BAM_LIST} -p ${TMP} -o ${PREFIX}"candidate_list_from_disc.txt" -n 8

def parse_toml_args():
    p = configargparse.ArgParser(default_config_files=['params.toml'],
                                 config_file_parser_class=configargparse.IniConfigParser(
                                     ['required_arguments','insertion_detection_params',
                                      'annotation_directories','output_directories'],split_ml_text_to_list=False),)


    p.add('-c', '--config', required=True, is_config_file=True, help='TOML config file path')
    p.add('--input_bams','-i',nargs = '+',required = True, help = "input bam(s)")
    p.add('--sample_name','-s',required = True, help = "Sample identifier")
    p.add('--repeat_type',required = False, default = ['ALU','L1','SVA'], nargs='+',
           help = 'Type of repeats to detect. Options include: L1, ALU, SVA, HERV, Mitochondrial (Default = ALU L1 SVA)')
    p.add('-m','--mode',required = False, default = 'germline',
           help = 'Which mode to run. Options: germline, mosaic, case-control, denovo (Default: germline)')
    p.add('-g','--genome',required = False, default = 'hg38',
          help = "Genome version used. Options: hg38 (default), hg19, chm13")

    # insertion detection paramaters:
    p.add('--cr',
          help="When specified, override default automatic calculation: cutoff of minimum # of clipped parts fall in repeats")
    p.add('--clip_cutoff',
          help="When specified, override default automatic calculation: cutoff of minimum # of clipped reads (used for both left & right sides)")
    p.add('--nd',
          help="When specified, override default automatic calculation: cutoff of minimum # of discordant pair")
    p.add('--tumor',default=False,
          help="Working on tumor samples")
    p.add('--purity', default=0.45,
          help="Tumor purity")
    p.add('--single', default=False,
          help="Call clip positions from single-end reads")
    p.add('--mosaic', default=False,
          help="Call mosaic events")

    # annotation directories:
    p.add("--rep_lib_annot_dir",required = True, help = 'Path to rep_lib_annotation/ directory')
    p.add("--genome_reference",required = True, help = 'Path to genome fasta file')
    p.add("--genome_gff3",required = True, help = 'Path to genome annotation file (.gff3 format)')

    # output directories:
    p.add("-o","--output_dir",required = True, help = 'Output directory',default = '.')
    p.add("--tmp_dir",required = False, help = 'Temporary directory',default = '.')


    p.add('-n',dest = 'cores',default = 1, help = "number of cores")
    p.add('--resume',default = True, help = "resume previous run if available")
    p.add('--save-intermediate-files',dest = "int_files",
           action='store_true', help = "Save all intermediate files")
    p.add('-v', help='version', action='store_true')

    return p.parse_args()

def setup_annotation_paths(rep,rep_lib_annot_dir,genome_reference,genome):

    annotation_paths = dict()

    if rep == 'ALU':
        r = 'Alu'
        annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/{r}/{genome}/{genome}_AluJabc_copies_with_flank.fa"
        annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus/{rep}.fa" #options.cns
        annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/{r}/{genome}/{genome}_{r}.out" #options.annotation
        annotation_paths["sf_anno1"] = f"{rep_lib_annot_dir}/{r}/{genome}/{genome}_Alu.out"
        annotation_paths["sf_flank"] = "null"
    elif rep == 'L1':
        annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/LINE/{genome}/{genome}_L1HS_copies_larger_5K_with_flank.fa"
        annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus/LINE1.fa" #options.cns
        annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/LINE/{genome}/{genome}_L1_larger_500_with_all_L1HS.out" #options.annotation
        annotation_paths["sf_anno1"] = f"{rep_lib_annot_dir}/LINE/{genome}/{genome}_L1.fa.out"
        annotation_paths["sf_flank"] = f"{rep_lib_annot_dir}/LINE/{genome}/{genome}_FL_L1_flanks_3k.fa"
    elif rep == 'SVA':
        annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_SVA_copies_with_flank.fa"
        annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus/{rep}.fa" #options.cns
        annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_{rep}.out" #options.annotation
        annotation_paths["sf_anno1"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_SVA.out"
        annotation_paths["sf_flank"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_FL_SVA_flanks_3k.fa"
    elif rep == 'HERV':
        annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_HERV_copies_with_flank.fa"
        annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus/{rep}.fa" #options.cns
        annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_{rep}.out" #options.annotation
        annotation_paths["sf_anno1"] = f"{rep_lib_annot_dir}/{r}/{genome}/{genome}_HERV.out"
        annotation_paths["sf_flank"] = "null"
    # elif rep == "MSTA":  ## TODO NOT TESTED
    #     annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_MSTA_copies_with_flank.fa"
    #     annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus/{rep}.fa" #options.cns
    #     annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_{rep}.out" #options.annotation
    #     #     sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "MSTA/hg38/hg38_MSTA.out"
    #     annotation_paths["sf_flank"] = "SF_FLANK null\n"
    # elif rep == "Pseudogene":
    #     annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/{rep}/{genome}/gencode_v28_GRCh38_transcript_masked.fa"
    #     annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus/gencode_v28_GRCh38_transcript_masked.fa" #options.cns
    #     annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/Pseudogene/{genome}/gencode_v28_GRCh38_exon_annotation.out" #options.annotation
    #     #     sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "Pseudogene/hg38/gencode_v28_GRCh38_exon_annotation.out"
    #     annotation_paths["sf_flank"] = "SF_FLANK null\n"
    # elif rep == "Mitochondrion":
    #     annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/Mitochondrion/hg38/hg38_mitochondrion_copies_with_flank.fa"
    #     annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus//mitochondrion.fa" #options.cns
    #     annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/Mitochondrion/hg38/hg38_numtS.out" #options.annotation
        
    annotation_paths["sf_ref"] = genome_reference

    return annotation_paths

def setup_output_dir(out_dir,tmp_dir,sample_name,repeat):
    o_dir = f'{out_dir}/{sample_name}/{repeat}/'
    t_dir = f'{out_dir}/{sample_name}/tmp/'

    Path(o_dir).mkdir(parents=True, exist_ok=True)
    Path(t_dir).mkdir(parents=True, exist_ok=True)

    o_abs_dir = str(Path(o_dir).resolve()) + '/'
    t_abs_dir = str(Path(t_dir).resolve()) + '/'

    return o_abs_dir, t_abs_dir



# BIG TODOS
#   - NO 10x support
#   - fflank file not used at all

if __name__ == '__main__':
    options = parse_toml_args()

    germline = True

    repeats = options.repeat_type

    if germline:
        for r in repeats:
            output_dir, sample_public_dir = setup_output_dir(options.output_dir,options.tmp_dir,options.sample_name,r)

            annot_path_dict = setup_annotation_paths(r,options.rep_lib_annot_dir,
                                                     options.genome_reference,
                                                     options.genome)
            
            #perform clipped step:
            print("Clipped reads step...")
            rcd,basic_rcd = get_clip_sites(options,annot_path_dict,output_dir,sample_public_dir)

            # perform discordant step:
            print("Discordant reads step...")
            get_disc_sites(options,annot_path_dict,output_dir,rcd,basic_rcd)

            # perform filter based on consensus seq:
            print("Consensus filter step...")
            filter_csn(options,annot_path_dict,output_dir,rcd,basic_rcd)

            #perform transduction step:
            print("Transduction step...")
            get_transduction(r,options,annot_path_dict,output_dir,rcd,basic_rcd)

            #sibling
            print("Sibling step...")
            get_sibling(r,options,annot_path_dict,output_dir,rcd,basic_rcd)
            
            #filter ( WHY DOES THIS HAPPEN 2x??? )
            print("Filter step...")
            filter_sites_post(r,options,annot_path_dict,output_dir,basic_rcd)
            # filter_sites_post(r,options,annot_path_dict,output_dir,basic_rcd)


            #annotate
            # time python ${XTEA_PATH}"x_TEA_main.py" --gene -a ${GENE} 
            # -i ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt"  -n 8 
            # -o ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"

            # time python ${XTEA_PATH}"x_TEA_main.py" --gntp_classify 
            # -i ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"  
            # -n 1 --model ${XTEA_PATH}"genotyping/DF21_model_1_2"  
            # -o ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"
            
            
            #output
            # time python ${XTEA_PATH}"x_TEA_main.py" --gVCF 
            # -i ${PREFIX}"candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"  
            # -o ${PREFIX} -b ${BAM_LIST} --ref ${REF} --rtype 2

    # elif 'mosaic':

    # elif 'case-control':

    # elif 'de-novo':