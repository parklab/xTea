##11/04/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
from subprocess import *
from optparse import OptionParser
import ntpath
import global_values

####
PUB_CLIP="pub_clip"
REP_TYPE_L1="L1"
REP_TYPE_ALU="Alu"
REP_TYPE_SVA="SVA"
REP_TYPE_HERV="HERV"
REP_TYPE_MIT="Mitochondrion"
REP_TYPE_MSTA="MSTA"

def gnrt_script_head(spartition, ncores, stime, smemory, s_id):
    s_head = "#!/bin/bash\n\n"
    s_head += "#SBATCH -n {0}\n".format(ncores)
    s_head += "#SBATCH -t {0}\n".format(stime)
    s_head += "#SBATCH --mem={0}G\n".format(smemory)
    s_head += "#SBATCH -p {0}\n".format(spartition)
    s_head += "#SBATCH -o {0}_%j.out\n".format(s_id)
    s_head += "#SBATCH --mail-type=NONE\n"
    s_head += "#SBATCH --mail-user=chong.simonchu@gmail.com\n"
    if spartition == "park" or spartition == "priopark":
        s_head += "#SBATCH --account=park_contrib\n\n"
    return s_head


# load in the parameter file or the configuration file
def load_par_config(sf_par_config):
    # by default, SF_FLANK is set to null, as Alu no need for SF_FLANK, as we don't check transduction for Alu
    l_pars = []
    with open(sf_par_config) as fin_par_config:
        for line in fin_par_config:
            if len(line) > 0 and line[0] == "#":
                continue
            fields = line.split()
            l_pars.append((fields[0], fields[1]))
    return l_pars

# gnrt parameters
def gnrt_parameters(l_pars):
    s_pars = ""
    for rcd in l_pars:
        sid = rcd[0]
        svalue = str(rcd[1])
        sline = sid + "=" + svalue + "\n"
        s_pars += sline
    return s_pars

####
# grnt calling steps
def gnrt_calling_command(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len, min_tei_len, iflag,
                         b_mosaic, b_user_par, b_force, i_rep_type, s_cfolder, b_SVA=False):
    s_user=""
    if b_user_par==True:
        s_user="--user"
    s_clean=""
    if b_force==True:
        s_clean="--force"
    sclip_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -C -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                 "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n".format(iclip_c, iclip_c, iclip_rp,
                                                                                            ncores, s_cfolder, s_user, s_clean)
    if b_SVA is True:
        sclip_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -C --sva -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                     "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n".format(iclip_c, iclip_c,
                                                                                                    iclip_rp, ncores,
                                                                                                    s_cfolder, s_user, s_clean)
    sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                 "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    if b_SVA == True:
        sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --sva -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                     "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    if b_mosaic == True:
        sclip_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --mosaic -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                     "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n".format(iclip_c, iclip_c,
                                                                                            iclip_rp, ncores, s_cfolder,
                                                                                            s_user, s_clean)
        sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --somatic -i ${{PREFIX}}\"candidate_list_from_clip.txt\" " \
                     "--nd {0} --ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    if b_mosaic==True and b_SVA ==True:
        sclip_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --mosaic --sva -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                     "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5} {6}\n".format(iclip_c, iclip_c,
                                                                                                    iclip_rp, ncores,
                                                                                                    s_cfolder, s_user, s_clean)
        sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --sva --somatic -i ${{PREFIX}}\"candidate_list_from_clip.txt\" " \
                     "--nd {0} --ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    sbarcode_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -B -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 " \
                    "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} " \
                    "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {0} {1}\n".format(ncores, s_user)
    sfilter_10x = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
                  "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
                  "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
                  "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len,
                                                                                   ncores, s_user)
    if b_SVA==True:
        sfilter_10x = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --sva --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
                      "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
                      "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
                      "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len,
                                                                                       ncores, s_user)
    s_filter = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
               "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
               "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
               "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    if b_SVA==True:
        s_filter = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --sva --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
                   "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
                   "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
                   "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len,
                                                                                    ncores, s_user)
    sf_collect = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -E --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} " \
                 "-i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} " \
                 "--flklen {1}\n".format(ncores, iflk_len)
    sf_asm = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -A -L -p ${{TMP}} --ref ${{REF}} -n {0} " \
             "-i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n".format(ncores)
    sf_alg_ctg = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -M -i ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\" " \
                 "--ref ${{REF}} -n {0} -p ${{TMP}} -r ${{L1_CNS}} " \
                 "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n".format(ncores)
    sf_mutation = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -I -p ${{TMP}} -n {0} " \
                  "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -r ${{L1_CNS}} " \
                  "--teilen {1} -o ${{PREFIX}}\"internal_snp.vcf.gz\"\n".format(ncores, min_tei_len)
    sf_gene="python ${{XTEA_PATH}}\"x_TEA_main.py\" --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns_with_gene.txt\"\n".format(ncores)

    ####

    sf_post_filter = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --rtype {0} -p ${{TMP_CNS}} " \
                     "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
                     "-a ${{ANNOTATION1}} " \
                     "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n".format(i_rep_type, ncores)

    if b_mosaic is True:
        sf_post_filter = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --postFmosaic --blacklist ${{BLACK_LIST}} " \
                         "--rtype {0} -p ${{TMP_CNS}} " \
                         "-n {1} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
                         "-a ${{ANNOTATION1}} " \
                         "-o ${{PREFIX}}\"candidate_disc_filtered_cns_post_filtering.txt\"\n".format(i_rep_type, ncores)

    sf_post_filter_hc = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --rtype {0} -p ${{TMP_CNS}} -n {1} " \
                        "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident\" -a ${{ANNOTATION1}} " \
                        "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n" \
        .format(i_rep_type, ncores)
    if b_mosaic is True:
        sf_post_filter_hc = "python ${{XTEA_PATH}}\"x_TEA_main.py\" --postF --postFmosaic --blacklist ${{BLACK_LIST}} " \
                            "--rtype {0} -p ${{TMP_CNS}} -n {1} " \
                            "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident\" -a ${{ANNOTATION1}} " \
                            "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt\"\n" \
            .format(i_rep_type, ncores)

    #sf_clean_tmp = "find ${TMP} -type f -name \'*tmp*\' -delete\n"
    # sf_clean_sam="find ${TMP} -type f -name \'*.sam\' -delete\n"
    # sf_clean_fa="find ${TMP} -type f -name \'*.fa\' -delete\n"
    # sf_clean_fq = "find ${TMP} -type f -name \'*.fq\' -delete\n"

    ####
    ####
    s_cmd = ""
    if iflag & 1 == 1:
        s_cmd += sclip_step
    if iflag & 2 == 2:
        s_cmd += sdisc_step
    if iflag & 4 == 4:
        s_cmd += sbarcode_step
    if iflag & 8 == 8:
        s_cmd += sfilter_10x
    if iflag & 16 == 16:
        s_cmd += s_filter
    if iflag & 32 == 32:
        s_cmd += sf_collect
    if iflag & 64 == 64:
        s_cmd += sf_asm
    if iflag & 128 == 128:
        s_cmd += sf_alg_ctg
    if iflag & 256 == 256:
        s_cmd += sf_mutation
    if iflag & 512 == 512:
        s_cmd += sf_gene
    if iflag & 1024 == 1024:
        s_cmd += sf_post_filter
        s_cmd += sf_post_filter_hc
    # s_cmd+=sf_clean_tmp
    # s_cmd +=sf_clean_sam
    # s_cmd +=sf_clean_fa
    # s_cmd +=sf_clean_fq
    return s_cmd
####
# grnt calling steps
def gnrt_calling_command_non_RNA(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len, min_tei_len, iflag,
                                 b_mosaic, b_user_par, s_cfolder):
    s_user = ""
    if b_user_par == True:
        s_user = "--user"
    sclip_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -C --dna -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                 "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5}\n".format(iclip_c, iclip_c, iclip_rp,
                                                                                            ncores, s_cfolder, s_user)
    sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --dna -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                 "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    if b_mosaic == True:
        sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --somatic --dna -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                     "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    sbarcode_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -B --dna -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 " \
                    "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} " \
                    "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {0} {1}\n".format(ncores, s_user)
    sfilter_10x = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
                  "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
                  "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
                  "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    s_filter = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
               "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
               "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
               "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    sf_collect = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -E --dna --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} " \
                 "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} " \
                 "--flklen {1}\n".format(ncores, iflk_len)
    sf_asm = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -A -L --dna -p ${{TMP}} --ref ${{REF}} -n {0} " \
             "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\"\n".format(ncores)
    sf_alg_ctg = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -M --dna -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
                 "--ref ${{REF}} -n {0} -p ${{TMP}} -r ${{L1_CNS}} " \
                 "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n".format(ncores)
    sf_mutation = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -I --dna -p ${{TMP}} -n {0} " \
                  "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -r ${{L1_CNS}} " \
                  "--teilen {1} -o ${{PREFIX}}\"internal_snp.vcf.gz\"\n".format(ncores, min_tei_len)
    sf_gene="python ${{XTEA_PATH}}\"x_TEA_main.py\" --gene -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns_with_gene.txt\"\n".format(ncores)
    #sf_clean_tmp = "find ${TMP} -type f -name \'*tmp*\' -delete\n"
    # sf_clean_sam="find ${TMP} -type f -name \'*.sam\' -delete\n"
    # sf_clean_fa="find ${TMP} -type f -name \'*.fa\' -delete\n"
    # sf_clean_fq = "find ${TMP} -type f -name \'*.fq\' -delete\n"

    ####
    ####
    s_cmd = ""
    if iflag & 1 == 1:
        s_cmd += sclip_step
    if iflag & 2 == 2:
        s_cmd += sdisc_step
    if iflag & 4 == 4:
        s_cmd += sbarcode_step
    if iflag & 8 == 8:
        s_cmd += sfilter_10x
    if iflag & 16 == 16:
        s_cmd += s_filter
    if iflag & 32 == 32:
        s_cmd += sf_collect
    if iflag & 64 == 64:
        s_cmd += sf_asm
    if iflag & 128 == 128:
        s_cmd += sf_alg_ctg
    if iflag & 256 == 256:
        s_cmd += sf_mutation
    if iflag & 512 == 512:
        s_cmd += sf_gene

    # s_cmd+=sf_clean_tmp
    # s_cmd +=sf_clean_sam
    # s_cmd +=sf_clean_fa
    # s_cmd +=sf_clean_fq
    return s_cmd

# grnt calling steps
def gnrt_calling_command_MT(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len, min_tei_len, iflag,
                            b_mosaic, b_user_par, s_cfolder):
    s_user = ""
    if b_user_par == True:
        s_user = "--user"
    sclip_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -C --dna --mit -i ${{BAM_LIST}} --lc {0} --rc {1} --cr {2}  " \
                 "-r ${{L1_COPY_WITH_FLANK}}  -a ${{ANNOTATION}} --cns ${{L1_CNS}} --ref ${{REF}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {3} --cp {4} {5}\n".format(iclip_c, iclip_c, iclip_rp,
                                                                                            ncores, s_cfolder, s_user)
    sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --dna --mit -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                 "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    if b_mosaic == True:
        sdisc_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\"  -D --somatic --dna --mit -i ${{PREFIX}}\"candidate_list_from_clip.txt\" --nd {0} " \
                     "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM_LIST}} -p ${{TMP}} " \
                     "-o ${{PREFIX}}\"candidate_list_from_disc.txt\" -n {1} {2}\n".format(idisc_c, ncores, s_user)
    sbarcode_step = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -B --dna --mit -i ${{PREFIX}}\"candidate_list_from_disc.txt\" --nb 400 " \
                    "--ref ${{REF}} -a ${{ANNOTATION}} -b ${{BAM1}} -d ${{BARCODE_BAM}} -p ${{TMP}} " \
                    "-o ${{PREFIX}}\"candidate_list_barcode.txt\" -n {0} {1}\n".format(ncores, s_user)
    sfilter_10x = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --mit --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
                  "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_barcode.txt\" " \
                  "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
                  "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    s_filter = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -N --dna --mit --cr {0} --nd {1} -b ${{BAM_LIST}} -p ${{TMP_CNS}} " \
               "--fflank ${{SF_FLANK}} --flklen {2} -n {3} -i ${{PREFIX}}\"candidate_list_from_disc.txt\" " \
               "-r ${{L1_CNS}} --ref ${{REF}} -a ${{ANNOTATION}} " \
               "-o ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" {4}\n".format(iflt_clip, iflt_disc, iflk_len, ncores, s_user)
    sf_collect = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -E --dna --mit --nb 500 -b ${{BAM1}} -d ${{BARCODE_BAM}} --ref ${{REF}} " \
                 "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -p ${{TMP}} -a ${{ANNOTATION}} -n {0} " \
                 "--flklen {1}\n".format(ncores, iflk_len)
    sf_asm = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -A -L --dna --mit -p ${{TMP}} --ref ${{REF}} -n {0} " \
             "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\"\n".format(ncores)
    sf_alg_ctg = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -M --dna --mit -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
                 "--ref ${{REF}} -n {0} -p ${{TMP}} -r ${{L1_CNS}} " \
                 "-o ${{PREFIX}}\"candidate_list_asm.txt\"\n".format(ncores)
    sf_mutation = "python ${{XTEA_PATH}}\"x_TEA_main.py\" -I --dna --mit -p ${{TMP}} -n {0} " \
                  "-i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" -r ${{L1_CNS}} " \
                  "--teilen {1} -o ${{PREFIX}}\"internal_snp.vcf.gz\"\n".format(ncores, min_tei_len)
    sf_gene="python ${{XTEA_PATH}}\"x_TEA_main.py\" --gene --dna --mit -a ${{GENE}} -i ${{PREFIX}}\"candidate_disc_filtered_cns.txt\" " \
            " -n {0} -o ${{PREFIX}}\"candidate_disc_filtered_cns_with_gene.txt\"\n".format(ncores)
    #sf_clean_tmp = "find ${TMP} -type f -name \'*tmp*\' -delete\n"
    # sf_clean_sam="find ${TMP} -type f -name \'*.sam\' -delete\n"
    # sf_clean_fa="find ${TMP} -type f -name \'*.fa\' -delete\n"
    # sf_clean_fq = "find ${TMP} -type f -name \'*.fq\' -delete\n"

    ####
    ####
    s_cmd = ""
    if iflag & 1 == 1:
        s_cmd += sclip_step
    if iflag & 2 == 2:
        s_cmd += sdisc_step
    if iflag & 4 == 4:
        s_cmd += sbarcode_step
    if iflag & 8 == 8:
        s_cmd += sfilter_10x
    if iflag & 16 == 16:
        s_cmd += s_filter
    if iflag & 32 == 32:
        s_cmd += sf_collect
    if iflag & 64 == 64:
        s_cmd += sf_asm
    if iflag & 128 == 128:
        s_cmd += sf_alg_ctg
    if iflag & 256 == 256:
        s_cmd += sf_mutation
    if iflag & 512 == 512:
        s_cmd += sf_gene

    # s_cmd+=sf_clean_tmp
    # s_cmd +=sf_clean_sam
    # s_cmd +=sf_clean_fa
    # s_cmd +=sf_clean_fq
    return s_cmd

####

####gnrt the whole pipeline
def gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_id, sf_bams, sf_bams_10X, sf_root_folder, rep_type):
    l_sbatch_files=[]
    sf_working_folder=sf_root_folder+rep_type
    if sf_working_folder[-1] != "/":
        sf_working_folder += "/"

    m_id = {}
    with open(sf_id) as fin_id:
        for line in fin_id:
            sid = line.rstrip()
            m_id[sid] = 1

            sf_folder = sf_working_folder
            if os.path.exists(sf_folder)==False:
                cmd = "mkdir {0}".format(sf_folder)
                run_cmd(cmd)
            # create the temporary folders
            sf_tmp=sf_folder + "/tmp"
            if os.path.exists(sf_tmp)==False:
                cmd = "mkdir {0}".format(sf_tmp)
                run_cmd(cmd)
            sf_tmp_clip=sf_folder + "/tmp/clip"
            if os.path.exists(sf_tmp_clip)==False:
                cmd = "mkdir {0}".format(sf_tmp_clip)
                run_cmd(cmd)
            sf_tmp_cns=sf_folder + "/tmp/cns"
            if os.path.exists(sf_tmp_cns)==False:
                cmd = "mkdir {0}".format(sf_tmp_cns)
                run_cmd(cmd)
    m_bams = {}
    if sf_bams != "null":
        with open(sf_bams) as fin_bams:
            for line in fin_bams:
                fields = line.split()
                sid = fields[0]
                s_bam = fields[1]
                if sid not in m_bams:
                    m_bams[sid]=[]
                m_bams[sid].append(s_bam)

    m_bams_10X = {}
    if sf_bams_10X != "null":
        with open(sf_bams_10X) as fin_bams_10X:
            for line in fin_bams_10X:
                fields = line.split()
                sid = fields[0]
                if sid not in m_id:
                    continue
                s_bam = fields[1]
                s_barcode_bam = fields[2]
                m_bams_10X[sid] = s_bam

                # soft-link the bams
                sf_10X_bam = sf_working_folder + "10X_phased_possorted_bam.bam"
                if os.path.isfile(sf_10X_bam) == False:
                    cmd = "ln -s {0} {1}".format(s_bam, sf_10X_bam)
                    run_cmd(cmd)

                sf_10X_barcode_bam = sf_working_folder + "10X_barcode_indexed.sorted.bam"
                if os.path.isfile(sf_10X_barcode_bam) == False:
                    cmd = "ln -s {0} {1}".format(s_barcode_bam, sf_10X_barcode_bam)  #
                    run_cmd(cmd)
                # soft-link the bai
                sf_10X_bai = sf_working_folder + "10X_phased_possorted_bam.bam.bai"
                if os.path.isfile(sf_10X_bai) == False:
                    cmd = "ln -s {0} {1}".format(s_bam + ".bai", sf_10X_bai)
                    run_cmd(cmd)

                sf_10X_barcode_bai = sf_working_folder + "10X_barcode_indexed.sorted.bam.bai"
                if os.path.isfile(sf_10X_barcode_bai) == False:
                    cmd = "ln -s {0} {1}".format(s_barcode_bam + ".bai", sf_10X_barcode_bai)
                    run_cmd(cmd)
                    ####
####
    for sid in m_id:
        sf_folder = sf_working_folder
        if os.path.exists(sf_folder) == False:
            continue

        ####gnrt the bam list file
        sf_bam_list = sf_folder + "bam_list.txt"
        with open(sf_bam_list, "w") as fout_bam_list:
            if sid in m_bams:
                for sf_tmp_bam in m_bams[sid]:
                    fout_bam_list.write(sf_tmp_bam + "\t" + global_values.ILLUMINA + "\n")
            if sid in m_bams_10X:
                fout_bam_list.write(m_bams_10X[sid] + "\t" + global_values.X10 + "\n")

        ####gnrt the pipeline file
        sf_out_sh = sf_folder + "run_xTEA_pipeline.sh"
        #sf_out_sh = sf_folder + "run_xTEA_pipeline1.sh"
        with open(sf_out_sh, "w") as fout_sh:  ###gnrt the pipeline file
            fout_sh.write(s_head)
            s_prefix = "PREFIX={0}\n".format(sf_folder)
            fout_sh.write(s_prefix)
            fout_sh.write("############\n")
            fout_sh.write("############\n")
            fout_sh.write(s_libs)
            fout_sh.write("############\n")
            fout_sh.write("############\n")
            fout_sh.write(s_calling_cmd)
        l_sbatch_files.append(sf_out_sh)
    return l_sbatch_files

#
def write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene, sf_black_list, sf_copy_with_flank, sf_flank, sf_cns,
                    sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config):
    with open(sf_config, "w") as fout_L1:
        fout_L1.write(sf_anno)
        fout_L1.write(sf_anno1)
        fout_L1.write(sf_ref)
        fout_L1.write(sf_gene)
        fout_L1.write(sf_black_list)
        fout_L1.write(sf_copy_with_flank)
        fout_L1.write(sf_flank)
        fout_L1.write(sf_cns)
        fout_L1.write(sf_xtea)
        fout_L1.write(s_bl)
        fout_L1.write(s_bam1)
        fout_L1.write(s_bc_bam)
        fout_L1.write(s_tmp)
        fout_L1.write(s_tmp_clip)
        fout_L1.write(s_tmp_cns)
####
#generate library configuration files
def gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, sf_config_prefix):
    if sf_folder_rep[-1] != "/":
        sf_folder_rep += "/"
    if sf_folder_xtea[-1] != "/":
        sf_folder_xtea += "/"
    if sf_config_prefix[-1] != "/":
        sf_config_prefix += "/"

    s_bl = "BAM_LIST ${PREFIX}\"bam_list.txt\"\n"
    s_bam1 = "BAM1 ${PREFIX}\"10X_phased_possorted_bam.bam\"\n"
    s_bc_bam = "BARCODE_BAM ${PREFIX}\"10X_barcode_indexed.sorted.bam\"\n"
    s_tmp = "TMP ${PREFIX}\"tmp/\"\n"
    s_tmp_clip = "TMP_CLIP ${PREFIX}\"tmp/clip/\"\n"
    s_tmp_cns = "TMP_CNS ${PREFIX}\"tmp/cns/\"\n"
    sf_ref = "REF " + sf_ref + "\n"
    sf_xtea = "XTEA_PATH " + sf_folder_xtea + "\n"
    sf_gene_anno="GENE " + sf_gene + "\n"
    sf_black_list="BLACK_LIST "+sf_black_list+"\n"
####
    for s_rep_type in l_rep_type:
        if s_rep_type=="L1":# for L1
            sf_config_L1 = sf_config_prefix + "L1.config"
            sf_anno = "ANNOTATION " + sf_folder_rep + "LINE/hg19/hg19_L1_larger500_with_all_L1HS.out\n"
            sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "LINE/hg19/hg19_L1.fa.out\n"
            sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep + "LINE/hg19/hg19_L1HS_copies_larger_5K_with_flank.fa\n"
            sf_flank = "SF_FLANK " + sf_folder_rep + "LINE/hg19/hg19_FL_L1_flanks_3k.fa\n"
            sf_cns = "L1_CNS " + sf_folder_rep + "consensus/LINE1.fa\n"
            write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_L1)
        elif s_rep_type=="Alu":#### for Alu
            sf_config_Alu = sf_config_prefix + "Alu.config"
            sf_anno = "ANNOTATION " + sf_folder_rep + "Alu/hg19/hg19_Alu.out\n"
            sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "Alu/hg19/hg19_Alu.out\n"
            sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep + "Alu/hg19/hg19_AluJabc_copies_with_flank.fa\n"
            sf_flank = "SF_FLANK null\n"
            sf_cns = "L1_CNS " + sf_folder_rep + "consensus/ALU.fa\n"
            write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_Alu)
        elif s_rep_type=="SVA":####for SVA
            sf_config_SVA = sf_config_prefix + "SVA.config"
            sf_anno = "ANNOTATION " + sf_folder_rep + "SVA/hg19/hg19_SVA.out\n"
            sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "SVA/hg19/hg19_SVA.out\n"
            sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep + "SVA/hg19/hg19_SVA_copies_with_flank.fa\n"
            sf_flank = "SF_FLANK " + sf_folder_rep + "SVA/hg19/hg19_FL_SVA_flanks_3k.fa\n"
            sf_cns = "L1_CNS " + sf_folder_rep + "consensus/SVA.fa\n"
            write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_SVA)
        elif s_rep_type=="HERV":####HERV
            sf_config_HERV = sf_config_prefix + "HERV.config"
            sf_anno = "ANNOTATION " + sf_folder_rep + "HERV/hg19/hg19_HERV.out\n"
            sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "HERV/hg19/hg19_HERV.out\n"
            sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep + "HERV/hg19/hg19_HERV_copies_with_flank.fa\n"
            sf_flank = "SF_FLANK null\n"
            sf_cns = "L1_CNS " + sf_folder_rep + "consensus/HERV.fa\n"
            write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_HERV)
        elif s_rep_type=="MSTA":####MSTA
            sf_config_MSTA = sf_config_prefix + "MSTA.config"
            sf_anno = "ANNOTATION " + sf_folder_rep + "MSTA/hg19/hg19_MSTA.out\n"
            sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "MSTA/hg19/hg19_MSTA.out\n"
            sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep + "MSTA/hg19/hg19_MSTA_copies_with_flank.fa\n"
            sf_flank = "SF_FLANK null\n"
            sf_cns = "L1_CNS " + sf_folder_rep + "consensus/MSTA.fa\n"
            write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_MSTA)

        elif s_rep_type=="Mitochondrion":####for Mitochondrion
            sf_config_Mitochondrion = sf_config_prefix + "Mitochondrion.config"
            sf_anno = "ANNOTATION " + sf_folder_rep + "Mitochondrion/hg19/hg19_numtS.out\n"
            sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "Mitochondrion/hg19/hg19_numtS.out\n"
            sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep +\
                                 "Mitochondrion/hg19/hg19_mitochondrion_copies_with_flank.fa\n"
            sf_flank = "SF_FLANK null\n"
            sf_cns = "L1_CNS " + sf_folder_rep + "consensus/mitochondrion.fa\n"
            write_to_config(sf_anno, sf_anno1, sf_ref, sf_gene_anno, sf_black_list, sf_copy_with_flank, sf_flank,
                            sf_cns, sf_xtea, s_bl, s_bam1, s_bc_bam, s_tmp, s_tmp_clip, s_tmp_cns, sf_config_Mitochondrion)

def cp_file(sf_from, sf_to):
    cmd = "cp {0} {1}".format(sf_from, sf_to)
    if os.path.isfile(sf_from)==False:
        return
    run_cmd(cmd)

def run_cmd(cmd):
    print cmd
    Popen(cmd, shell=True, stdout=PIPE).communicate()

def get_sample_id(sf_bam):
    fname = ntpath.basename(sf_bam)
    fname_fields = fname.split(".")
    if fname_fields[-1] != "bam" and fname_fields[-1] != "cram":
        print "Alignment is not end with .bam"
        return None
    sample_id = ".".join(fname_fields[:-1])
    return sample_id


####gnrt the running shell
def gnrt_running_shell(sf_ids, sf_bams, sf_10X_bams, l_rep_type, b_mosaic, b_user_par, b_force, s_wfolder, sf_folder_rep,
                       sf_ref, sf_gene, sf_black_list, sf_folder_xtea, spartition, stime, smemory, ncores, sf_submit_sh):
    if s_wfolder[-1] != "/":
        s_wfolder += "/"
    if os.path.exists(s_wfolder) == False:
        scmd = "mkdir {0}".format(s_wfolder)
        Popen(scmd, shell=True, stdout=PIPE).communicate()

    m_id = {}
    with open(sf_ids) as fin_id:
        for line in fin_id:
            sid = line.rstrip()
            m_id[sid] = 1
            sf_folder = s_wfolder + sid  # first creat folder
            if os.path.exists(sf_folder) == True:
                continue
            cmd = "mkdir {0}".format(sf_folder)
            Popen(cmd, shell=True, stdout=PIPE).communicate()
            sf_pub_clip = sf_folder+ "/" + PUB_CLIP + "/"
            cmd = "mkdir {0}".format(sf_pub_clip)
            Popen(cmd, shell=True, stdout=PIPE).communicate()

    ####gnrt the sample, bam, x10 bam files
    split_bam_list(m_id, sf_bams, sf_10X_bams, l_rep_type, s_wfolder)

    l_sh=[]
    for sid_tmp in m_id:
        sf_sample_folder=s_wfolder + sid_tmp + "/"
        sf_pub_clip = sf_sample_folder + PUB_CLIP + "/"
        gnrt_lib_config(l_rep_type, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, sf_sample_folder)
        for rep_type in l_rep_type:
            i_rep_type = get_flag_by_rep_type(rep_type)
            sf_config = sf_sample_folder+rep_type+".config"

            #s_wfolder_rep=s_wfolder+rep_type
            s_wfolder_rep = sf_sample_folder + rep_type
            if os.path.exists(s_wfolder_rep)==False:
                cmd="mkdir {0}".format(s_wfolder_rep)
                run_cmd(cmd)

            sf_rep_sample_id = s_wfolder_rep + "/sample_id.txt"
            sf_rep_bam = s_wfolder_rep + "/bam_list1.txt"
            if os.path.isfile(sf_rep_bam)==False:
                sf_rep_bam="null"
            sf_rep_x10_bam = s_wfolder_rep + "/x10_bam_list1.txt"
            if os.path.exists(sf_rep_x10_bam)==False:
                sf_rep_x10_bam="null"


            s_head = gnrt_script_head(spartition, ncores, stime, smemory, sid_tmp)
            l_libs = load_par_config(sf_config)
            s_libs = gnrt_parameters(l_libs)
            ##
            iclip_c = options.nclip
            iclip_rp = options.cliprep
            idisc_c = options.ndisc
            iflt_clip = options.nfilterclip
            iflt_disc = options.nfilterdisc
            iflk_len = options.flklen
            itei_len = options.teilen
            iflag = options.flag

            s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
                                                 itei_len, iflag, b_mosaic, b_user_par, b_force, i_rep_type, sf_pub_clip)
            if rep_type==REP_TYPE_SVA:
                s_calling_cmd = gnrt_calling_command(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
                                                     itei_len, iflag, b_mosaic, b_user_par, b_force, i_rep_type, sf_pub_clip, True)
            if rep_type==REP_TYPE_MSTA or rep_type==REP_TYPE_HERV:
                s_calling_cmd = gnrt_calling_command_non_RNA(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores,
                                                             iflk_len, itei_len, iflag, b_mosaic, b_user_par, sf_pub_clip)
            elif rep_type==REP_TYPE_MIT:
                s_calling_cmd = gnrt_calling_command_MT(iclip_c, iclip_rp, idisc_c, iflt_clip, iflt_disc, ncores, iflk_len,
                                                     itei_len, iflag, b_mosaic, b_user_par, sf_pub_clip)

            l_tmp_sh=gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_rep_sample_id, sf_rep_bam, sf_rep_x10_bam,
                           sf_sample_folder, rep_type)
            for tmp_sh in l_tmp_sh:
                l_sh.append(tmp_sh)
    with open(sf_submit_sh, "w") as fout_submit:
        fout_submit.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            fout_submit.write("sbatch < "+s_sh+"\n")
####

####Input:
# m_ids: sample id dictionary
def split_bam_list(m_ids, sf_bams, sf_x10_bams, l_rep_type, s_wfolder):
    #load in sf_bams
    m_bams = {}
    if sf_bams != "null":
        with open(sf_bams) as fin_bams:
            for line in fin_bams:
                fields = line.split()
                sid = fields[0]
                if sid not in m_ids:
                    continue
                s_bam = fields[1]
                if sid not in m_bams:
                    m_bams[sid]=[]
                m_bams[sid].append(s_bam)

    #load in sf_x10_bams
    m_bams_10X = {}
    if sf_bams_10X != "null":
        with open(sf_x10_bams) as fin_bams_10X:
            for line in fin_bams_10X:
                fields = line.split()
                sid = fields[0]
                if sid not in m_ids:
                    continue
                s_bam = fields[1]
                s_barcode_bam = fields[2]
                m_bams_10X[sid] = (s_bam, s_barcode_bam)

    for sid_tmp in m_ids:
        sf_sample_folder = s_wfolder + sid_tmp + "/"
        if os.path.exists(sf_sample_folder) == False:
            cmd = "mkdir {0}".format(sf_sample_folder)
            run_cmd(cmd)
        for rep_type in l_rep_type:
            s_wfolder_rep = sf_sample_folder + rep_type
            if os.path.exists(s_wfolder_rep)==False:
                cmd="mkdir {0}".format(s_wfolder_rep)
                run_cmd(cmd)
            # need to generate two files for each sample
            sf_rep_sample_id = s_wfolder_rep + "/sample_id.txt"
            with open(sf_rep_sample_id, "w") as fout_rep_sample_id:
                fout_rep_sample_id.write(sid_tmp)
            if sid_tmp in m_bams:
                sf_rep_bam = s_wfolder_rep + "/bam_list1.txt"
                with open(sf_rep_bam, "w") as fout_rep_bams:
                    for sf_tmp_bam in m_bams[sid_tmp]:
                        fout_rep_bams.write(sid_tmp+"\t"+sf_tmp_bam+"\n")
            if sid_tmp in m_bams_10X:
                sf_rep_x10_bam = s_wfolder_rep + "/x10_bam_list1.txt"
                with open(sf_rep_x10_bam, "w") as fout_rep_x10_bams:
                    fout_rep_x10_bams.write(sid_tmp+"\t"+m_bams_10X[sid_tmp][0]+"\t"+m_bams_10X[sid_tmp][1]+"\n")


####run the pipelines
def run_pipeline(l_rep_type, sample_id, s_wfolder):
    for rep_type in l_rep_type:
        sf_sbatch_sh_rep = s_wfolder + sample_id + "/run_{0}.sh".format(rep_type)
        cmd="sh {0}".format(sf_sbatch_sh_rep)
        run_cmd(cmd)

def decompress(sf_in_tar, sf_out):
    cmd="tar -zxvf {0} -C {1}".format(sf_in_tar, sf_out)
    run_cmd(cmd)

def cp_compress_results(s_wfolder, l_rep_type, sample_id):
    # create a "results" folder
    sf_rslts = s_wfolder + "results/"
    if os.path.exists(sf_rslts)==False:
        cmd = "mkdir {0}".format(sf_rslts)
        run_cmd(cmd)

    for rep_type in l_rep_type:
        sf_rslts_rep_folder=sf_rslts+rep_type+"/"
        if os.path.exists(sf_rslts_rep_folder)==False:
            cmd = "mkdir {0}".format(sf_rslts_rep_folder)
            run_cmd(cmd)
        sf_samp_folder = sf_rslts_rep_folder + sample_id + "/"
        if os.path.exists(sf_samp_folder)==False:
            cmd="mkdir {0}".format(sf_samp_folder)
            run_cmd(cmd)

        sf_source_folder=s_wfolder+sample_id+"/"+rep_type+"/"
        sf_rslt1=sf_source_folder+"candidate_disc_filtered_cns.txt"
        cp_file(sf_rslt1, sf_samp_folder)
        sf_rslt2 = sf_source_folder + "candidate_list_from_clip.txt"
        cp_file(sf_rslt2, sf_samp_folder)
        sf_rslt2_1 = sf_source_folder + "candidate_list_from_clip.txt_tmp"
        cp_file(sf_rslt2_1, sf_samp_folder)
        sf_rslt3 = sf_source_folder + "candidate_list_from_disc.txt"
        cp_file(sf_rslt3, sf_samp_folder)

        # s_tmp1=sf_source_folder+"tmp/cns/candidate_sites_all_disc.fa"
        # cp_file(s_tmp1, sf_samp_folder)
        # s_tmp2 = sf_source_folder + "tmp/cns/candidate_sites_all_clip.fq"
        # cp_file(s_tmp2, sf_samp_folder)

        s_tmp1=sf_source_folder+"tmp/cns/temp_disc.sam"
        cp_file(s_tmp1, sf_samp_folder)
        s_tmp2 = sf_source_folder + "tmp/cns/temp_clip.sam"
        cp_file(s_tmp2, sf_samp_folder)
        s_tmp3 = sf_source_folder + "tmp/cns/all_with_polymerphic_flanks.fa"
        cp_file(s_tmp3, sf_samp_folder)
        s_tmp4 = sf_source_folder + "tmp/clip_reads_tmp0"
        cp_file(s_tmp4, sf_samp_folder)
        s_tmp5 = sf_source_folder + "tmp/discordant_reads_tmp0"
        cp_file(s_tmp5, sf_samp_folder)

    #compress the results folder to one file
    sf_compressed=sf_rslts+"results.tar.gz"
    cmd="tar -cvzf {0} -C {1} .".format(sf_compressed, sf_rslts)
    run_cmd(cmd)

def get_flag_by_rep_type(s_type):
    if s_type is REP_TYPE_L1:
        return 1
    elif s_type is REP_TYPE_ALU:
        return 2
    elif s_type is REP_TYPE_SVA:
        return 4
    elif s_type is REP_TYPE_HERV:
        return 8
    elif s_type is REP_TYPE_MIT:
        return 16
    elif s_type is REP_TYPE_MSTA:
        return 32
####
####
def parse_option():
    parser = OptionParser()
    parser.add_option("-D", "--decompress",
                      action="store_true", dest="decompress", default=False,
                      help="Decompress the rep lib and reference file")
    parser.add_option("-M", "--mosaic",
                      action="store_true", dest="mosaic", default=False,
                      help="Calling mosaic events from high coverage data")
    parser.add_option("-U", "--user",
                      action="store_true", dest="user", default=False,
                      help="Use user specific parameters instead of automatically calculated ones")
    parser.add_option("--force",
                      action="store_true", dest="force", default=False,
                      help="Force to start from the very beginning")

    parser.add_option("-i", "--id", dest="id",
                      help="sample id list file ", metavar="FILE")
    parser.add_option("-a", "--par", dest="parameters",
                      help="parameter file ", metavar="FILE")
    parser.add_option("-l", "--lib", dest="lib",
                      help="TE lib config file ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-x", "--x10", dest="x10",
                      help="Input 10X bam file", metavar="FILE")

    parser.add_option("-n", "--cores", dest="cores", type="int", default=8,
                      help="number of cores")
    parser.add_option("-m", "--memory", dest="memory", type="int", default=25,
                      help="Memory limit in GB")
    parser.add_option("-q", "--partition", dest="partition", type="string", default="medium",
                      help="Which queue to run the job")
    parser.add_option("-t", "--time", dest="time", type="string",default="0-15:00",
                      help="Time limit")

    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-r", "--ref", dest="ref", type="string",
                      help="reference genome")
    parser.add_option("-g", "--gene", dest="gene", type="string",
                      help="Gene annotation file")
    parser.add_option("--xtea", dest="xtea", type="string",
                      help="xTEA folder")

    parser.add_option("-f", "--flag", dest="flag", type="int", default=19,
                      help="Flag indicates which step to run (1-clip, 2-disc, 4-barcode, 8-xfilter, 16-filter, 32-asm)")

    parser.add_option("-y", "--reptype", dest="rep_type", type="int",
                      help="Type of repeats working on: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrial")

    parser.add_option("--flklen", dest="flklen", type="int",
                      help="flank region file")
    parser.add_option("--nclip", dest="nclip", type="int",
                      help="cutoff of minimum # of clipped reads")
    parser.add_option("--cr", dest="cliprep", type="int",
                      help="cutoff of minimum # of clipped reads whose mates map in repetitive regions")
    parser.add_option("--nd", dest="ndisc", type="int",
                      help="cutoff of minimum # of discordant pair")
    parser.add_option("--nfclip", dest="nfilterclip", type="int",
                      help="cutoff of minimum # of clipped reads in filtering step")
    parser.add_option("--nfdisc", dest="nfilterdisc", type="int",
                      help="cutoff of minimum # of discordant pair of each sample in filtering step")
    parser.add_option("--teilen", dest="teilen", type="int",
                      help="minimum length of the insertion for future analysis")

    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("--blacklist", dest="blacklist", default="null",
                      help="Reference panel database for filtering, or a blacklist region", metavar="FILE")
    (options, args) = parser.parse_args()
    return (options, args)
####
if __name__ == '__main__':
    (options, args) = parse_option()
    sf_id = options.id
    sf_bams = options.bam ###input is a bam file
    sf_bams_10X = options.x10
    sf_sbatch_sh = options.output  # this is the shell for submitting the jobs
    s_wfolder = options.wfolder
    b_mosaic=options.mosaic
    b_user_par=options.user
    b_force=options.force

    if s_wfolder[-1]!="/":
        s_wfolder+="/"
    if os.path.exists(s_wfolder) == False:
        scmd = "mkdir {0}".format(s_wfolder)
        Popen(scmd, shell=True, stdout=PIPE).communicate()

    if os.path.isfile(sf_bams) == False:
        sf_bams = "null"
    if os.path.isfile(sf_bams_10X) == False:
        sf_bams_10X = "null"
####
    spartition = options.partition
    stime = options.time
    smemory = options.memory
    ncores = options.cores
    sf_folder_rep1 = options.lib  ##this is the lib folder path
    sf_ref1=options.ref ####reference genome
    sf_folder_xtea=options.xtea

    sf_folder_rep=sf_folder_rep1
    sf_ref=sf_ref1
    if options.decompress==True:
        decompress(sf_folder_rep1, s_wfolder)
        decompress(sf_ref1, s_wfolder)
        sf_folder_rep = s_wfolder+"rep_lib_annotation/"
        sf_ref=s_wfolder+"genome.fa"
    sf_gene=options.gene
    sf_black_list=options.blacklist

    i_rep_type=options.rep_type
    l_rep_type = []
    if i_rep_type & 1 != 0:
        l_rep_type.append(REP_TYPE_L1)
    if i_rep_type & 2 != 0:
        l_rep_type.append(REP_TYPE_ALU)
    if i_rep_type & 4 != 0:
        l_rep_type.append(REP_TYPE_SVA)
    if i_rep_type & 8 != 0:
        l_rep_type.append(REP_TYPE_HERV)
    if i_rep_type & 16 != 0:
        l_rep_type.append(REP_TYPE_MIT)
    if i_rep_type & 32 != 0:
        l_rep_type.append(REP_TYPE_MSTA)

    gnrt_running_shell(sf_id, sf_bams, sf_bams_10X, l_rep_type, b_mosaic, b_user_par, b_force, s_wfolder, sf_folder_rep,
                       sf_ref, sf_gene, sf_black_list, sf_folder_xtea, spartition, stime, smemory, ncores, sf_sbatch_sh)


    #run_pipeline(l_rep_type, sample_id, s_wfolder)
    #cp_compress_results(s_wfolder, l_rep_type, sample_id)

####
####