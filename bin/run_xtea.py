##05/08/2024
##@@author: Corinne Sexton, DBMI, Harvard Medical School
##@@contact: corinne_sexton@hms.harvard.edu

from pathlib import Path
import configargparse
from xtea.locate_insertions import get_clipped_reads


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
           help = 'Type of repeats to detect. Options include: L1, ALU, SVA, HERV, Mitochondrial (Default = Alu L1 SVA)')
    p.add('-m','--mode',required = False, default = 'germline',
           help = 'Which mode to run. Options: germline, mosaic, case-control, denovo (Default: germline)')
    p.add('-g','--genome',required = False, default = 'hg38',
          help = "Genome version used. Options: hg38 (default), hg19, chm13")

    # insertion detection paramaters:
    p.add('--cr',
          help="When specified, override default automatic calculation: cutoff of minimum # of clipped parts fall in repeats")
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

    # output directories:
    p.add("-o","--output_dir",required = True, help = 'Output directory',default = '.')
    p.add("--tmp_dir",required = False, help = 'Temporary directory',default = '.')


    p.add('-n',dest = 'cores',default = 1, help = "number of cores")
    p.add('--resume',default = True, help = "resume previous run if available")
    p.add('-v', help='version', action='store_true')

    return p.parse_args()

def setup_annotation_paths(rep,rep_lib_annot_dir,genome_reference,genome):
    annotation_paths = dict()
    annotation_paths["sf_rep"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_AluJabc_copies_with_flank.fa" #options.reference TODO!!!!
    annotation_paths["sf_rep_cns"] = f"{rep_lib_annot_dir}/consensus/{rep}.fa" #options.cns
    annotation_paths["sf_annotation"] = f"{rep_lib_annot_dir}/{rep}/{genome}/{genome}_{rep}.out" #options.annotation
    annotation_paths["sf_ref"] = genome_reference

    return annotation_paths

def setup_output_dir(out_dir,tmp_dir,sample_name,repeat):
    if out_dir != '.':
        o_dir = f'{out_dir}/{sample_name}/{repeat}'
    else:
        o_dir = f'{sample_name}/{repeat}'

    if tmp_dir != '.':
        t_dir = f'{out_dir}/{sample_name}/tmp'
    else:
        t_dir = f'{sample_name}/tmp'

    Path(o_dir).mkdir(parents=True, exist_ok=True)
    Path(t_dir).mkdir(parents=True, exist_ok=True)

    return o_dir, t_dir


if __name__ == '__main__':
    options = parse_toml_args()

    germline = True

    repeats = options.repeat_type

    if germline:
        for r in repeats:
            output_dir, tmp_dir = setup_output_dir(options.output_dir,options.tmp_dir,options.sample_name,r)

            annot_path_dict = setup_annotation_paths(r,options.rep_lib_annot_dir,
                                                     options.genome_reference,
                                                     options.genome)
            
            #perform clip step:
            get_clipped_reads(options,r,annot_path_dict,output_dir,tmp_dir)

            #clip
            #disc
            #transduction
            #sibling
            #filter
            #annotate
            #output

    # elif 'mosaic':

    # elif 'case-control':

    # elif 'de-novo':