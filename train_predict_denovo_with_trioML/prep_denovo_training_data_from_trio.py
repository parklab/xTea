##06/12/2025
##@@author: Simon (Chong) Chu, The Wistar Institute
##@@contact: chong.simon.chu@gmail.com
####
#this is to prepare the training data for de novo TE insertion calling from trio data
#mainly prepare shuffled images of TP insertions

import os
import sys
from optparse import OptionParser
from x_genotype_classify import *
from x_BamSnap import *


####
def prep_cmds_to_generate_positive_training_images(sf_bam_list, sf_child_folder_list, sf_pa_folder_list, sf_ma_folder_list, sf_wfolder, s_type, sf_ref, i_extnd, sf_tp_manual_inspec, sf_out_cmd):
    """
    Prepare the training data for de novo TE insertion calling from trio data.
    :param sf_child_folder_list: Path to the file containing child folder list.
    :param sf_pa_folder_list: Path to the file containing paternal folder list.
    :param sf_ma_folder_list: Path to the file containing maternal folder list.
    :param sf_out_arff: Output ARFF file path.
    """
    m_proband_folders = load_in_folders(sf_child_folder_list)
    m_pa_folders = load_in_folders(sf_pa_folder_list)
    m_ma_folders = load_in_folders(sf_ma_folder_list)
    m_bam_list = load_in_bam_list(sf_bam_list)
    m_tp_manual_inspec=load_in_manual_inspected_TP_insertions(sf_tp_manual_inspec)
    
    ####
    if len(sf_wfolder)>0 and sf_wfolder[-1]!="/":
        sf_wfolder+="/"
    sf_wfolder_positive=sf_wfolder+"positive/"
    if os.path.exists(sf_wfolder_positive)==False:
        os.mkdir(sf_wfolder_positive)
    sf_wfolder_negative=sf_wfolder+"negative/"
    if os.path.exists(sf_wfolder_negative)==False:
        os.mkdir(sf_wfolder_negative)
    ####
    xbamsnap= XBamSnap()
    print("Finished loading folders and bams")

    #1. Define the high quality heterozygous insertion in proband
    #   1.1 Make sure one location only recruite one insertion (only select one for high population AF insertions)
    m_all_slcted_01={} #{chrm: position: family_id: [features]}
    m_negative_cmd={}
    i_check_slack=30 #no nearby insertions within 30bp
    for s_family in m_proband_folders:
        if (s_family in m_pa_folders) and (s_family in m_ma_folders):
            sf_child=m_proband_folders[s_family]
            sf_pa=m_pa_folders[s_family]
            sf_ma=m_ma_folders[s_family]
            m_slcted_01=select_high_quality_heterozygous_insertions(sf_child, sf_pa, sf_ma, s_type)
            print("working on sample %s to select high quality insertions"%s_family)
            ####
            for s_chrm in m_slcted_01:
                if s_chrm not in m_all_slcted_01:
                    m_all_slcted_01[s_chrm]={}
                for i_pos in m_slcted_01[s_chrm]:
                    b_exist=False
                    for i_tmp in range(i_pos-i_check_slack, i_pos+i_check_slack):#
                        if i_tmp in m_all_slcted_01[s_chrm]:
                            print("Warning: %s_%d already has a insertion, skip this one" % (s_chrm, i_tmp))
                            b_exist=True
                            break
                    if b_exist==False:
                        m_all_slcted_01[s_chrm][i_pos] = s_family  # store family id and features
            ####
            ####get the FPs from xTea "de novo" outputs
            ####
            m_negative_sites=select_fp_from_xtea_raw_denovo_insertions(sf_child, s_type, m_tp_manual_inspec, i_check_slack)
            if s_family not in m_bam_list:
                continue
            sf_pa_bam = m_bam_list[s_family][1]  # paternal BAM
            sf_ma_bam = m_bam_list[s_family][2]  # maternal BAM
            sf_proband_bam = m_bam_list[s_family][0] # proband BAM
            for s_tmp_chrm in m_negative_sites:
                for i_tmp_pos in m_negative_sites[s_tmp_chrm]:
                    s_cmd, s_cmd2=xbamsnap.prepare_bamsnap_scripts_one_denovo_loci(s_family, s_tmp_chrm, i_tmp_pos, i_extnd, sf_pa_bam, sf_ma_bam, sf_proband_bam, sf_wfolder_negative, sf_ref)
                    m_negative_cmd[s_family]=(s_cmd, s_cmd2)####
    
    print("Finished selecting high quality heterozygous insertions")

    # output the negative command
    with open(sf_out_cmd+".negative.sh","w") as fout_cmd:
        for s_family in m_negative_cmd:
            s_cmd, s_cmd2=m_negative_cmd[s_family]
            fout_cmd.write(s_cmd)
            fout_cmd.write(s_cmd2)

    #2. For each heterozygous insertion in proband, find two pseudo "parents" that do not have the insertion
    m_parents_clip=load_in_clip_files_all_trios_parents(m_pa_folders, m_ma_folders, s_type)
    with open(sf_out_cmd, "w") as fout_cmd, open(sf_out_cmd+".log", "w") as fout_log:
        for s_chrm in m_all_slcted_01:
            for i_pos in m_all_slcted_01[s_chrm]:
                #s_id = s_chrm + "_" + str(i_pos)
                s_family = m_all_slcted_01[s_chrm][i_pos]
                for s_tmp_family in m_parents_clip:
                    if s_tmp_family == s_family:
                        continue
                    if (s_tmp_family not in m_bam_list) or (s_family not in m_bam_list):
                        print("Warning: %s not in bam list, skip this family" % s_tmp_family)
                        continue
                    m_pa_clip, m_ma_clip = m_parents_clip[s_tmp_family]
                    b_in_pa_clip, s_pa_clip_chrm, i_pa_clip_pos = has_overlap(s_chrm, i_pos, m_pa_clip, 35)
                    b_in_ma_clip, s_ma_clip_chrm, i_ma_clip_pos = has_overlap(s_chrm, i_pos, m_ma_clip, 35)
                    
                    b_paired=False
                    if b_in_pa_clip == False and b_in_ma_clip == False:
                        b_paired=True
                        #generate the cmd to generate the image
                        sf_pa_bam = m_bam_list[s_tmp_family][1]  # paternal BAM
                        sf_ma_bam = m_bam_list[s_tmp_family][2]  # maternal BAM
                        sf_proband_bam = m_bam_list[s_family][0] # proband BAM
                        #here the proband is the one in this "s_family", and the parents are pseudo parents
                        s_cmd, s_cmd2=xbamsnap.prepare_bamsnap_scripts_one_denovo_loci(s_family, s_chrm, i_pos, i_extnd, sf_pa_bam, sf_ma_bam, sf_proband_bam, sf_wfolder_positive, sf_ref)
                        fout_cmd.write(s_cmd)
                        fout_cmd.write(s_cmd2)
                        fout_log.write(s_chrm+"\t"+str(i_pos)+"\t"+s_family+"\t"+s_tmp_family+"\n")
                        break #each insertion we only need one pair of pseudo parents, so break the loop after finding one
                    if b_paired==False:
                        fout_log.write(s_chrm+"\t"+str(i_pos)+"\t"+s_family+"\tnot_paired\n")

####
####generate false positive ones
def prep_cmds_to_generate_negative_training_images():
    #candidate_disc_filtered_cns_post_filtering_denovo.txt
    pass
####

####
def load_in_bam_list(sf_bam_list):
    """
    Load in the BAM list from a file.
    
    :param sf_bam_list: Path to the BAM list file.
    :return: A dictionary mapping sample IDs to their corresponding BAM files.
    """
    m_bams = {}
    with open(sf_bam_list, 'r') as fin_bams:
        for line in fin_bams:
            fields = line.strip().split()
            if len(fields) < 4:
                print("Wrong bam file list ", line)
                continue
            l_sf_bams = []
            b_exist = True
            for sf_tmp_bam in fields[1:]:
                if not os.path.isfile(sf_tmp_bam):
                    print("BAM file does not exist: %s" % sf_tmp_bam)
                    b_exist = False
                    break
                l_sf_bams.append(sf_tmp_bam)
            if b_exist:
                m_bams[fields[0]] = l_sf_bams
    return m_bams
####

def load_in_folders(sf_list):
    m_folders={}
    with open(sf_list) as fin_folders:
        for line in fin_folders:
            fields=line.rstrip().split()
            s_sample=fields[0]
            sf_folder=fields[1]
            if len(sf_folder) > 0 and sf_folder[-1] != "/":
                sf_folder += "/"
            m_folders[s_sample]=sf_folder
    return m_folders


def load_in_clip_files_all_trios_parents(m_pa_folders, m_ma_folders, s_type):
    """
    Load in all the trio parents from the paternal and maternal folders.
    
    :param m_pa_folders: Dictionary of paternal folders.
    :param m_ma_folders: Dictionary of maternal folders.
    :return: A dictionary containing all the trio parents.
    """
    m_parents_clip = {}
    for s_family in m_pa_folders:
        if s_family in m_ma_folders:
            sf_pa = m_pa_folders[s_family]
            sf_ma = m_ma_folders[s_family]
            m_pa_clip, m_ma_clip = load_in_clip_files_parents_one_trio(sf_pa, sf_ma, s_type)
            if len(m_pa_clip) > 0 and len(m_ma_clip) > 0:
                m_parents_clip[s_family] = (m_pa_clip, m_ma_clip)
            else:
                print("Failed to load clip files for family: %s" % s_family)
    return m_parents_clip


#load in all the clip files from the child, paternal and maternal folders
def load_in_clip_files_parents_one_trio(sf_pa, sf_ma, s_type):
    s_clip_name = "candidate_list_from_clip.txt"
    sf_pa_clip = sf_pa + s_type+"/"+s_clip_name
    sf_ma_clip = sf_ma + s_type+"/"+s_clip_name

    if os.path.isfile(sf_pa_clip) == False or os.path.isfile(sf_ma_clip) == False:
        return {}, {}

    m_pa_clip = load_in_xtea_output(sf_pa_clip)
    m_ma_clip = load_in_xtea_output(sf_ma_clip)

    return m_pa_clip, m_ma_clip
####


#select the high quality heterozygous insertions from child
def select_high_quality_heterozygous_insertions(sf_child, sf_pa, sf_ma, s_type):
    #print("Working on %s"%sf_child)
    s_hq_name="candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt"
    sf_child_hq=sf_child+s_type+"/"+s_hq_name
    sf_pa_hq = sf_pa + s_type+"/"+s_hq_name
    sf_ma_hq = sf_ma + s_type+"/"+s_hq_name
    if os.path.isfile(sf_child_hq)==False or os.path.isfile(sf_pa_hq)==False or os.path.isfile(sf_ma_hq)==False:
        print("%s or %s or %s failed to have high quality results" % (sf_child_hq, sf_pa_hq, sf_ma_hq))
        return {}
    ####
    ####
    #s_clip_name = "candidate_list_from_clip.txt"
    s_clip_name="candidate_list_from_disc.txt"
    sf_pa_clip = sf_pa + s_type+"/"+s_clip_name
    sf_ma_clip = sf_ma + s_type+"/"+s_clip_name
    if os.path.isfile(sf_pa_clip) == False or os.path.isfile(sf_ma_clip) == False:
        print("%s failed to have paternal/maternal clip file"%sf_pa_clip)
        return {}

    m_child_hq_tprt_both=load_in_xtea_output_hq(sf_child_hq, "two_side_tprt")
    m_pa_hq_tprt_both=load_in_xtea_output(sf_pa_hq, "two_side_tprt")
    m_ma_hq_tprt_both = load_in_xtea_output(sf_ma_hq, "two_side_tprt")

    m_pa_clip=load_in_xtea_output(sf_pa_clip)
    m_ma_clip = load_in_xtea_output(sf_ma_clip)

    #l_slcted_01=[]
    m_slcted_01={}
    for s_chrm in m_child_hq_tprt_both:
        for i_pos in m_child_hq_tprt_both[s_chrm]:
            s_id = s_chrm + "_" + str(i_pos)
            l_fields = m_child_hq_tprt_both[s_chrm][i_pos]
            l_features = load_in_features_one_insertion(l_fields)

            b_in_pa_hq, s_pa_hq_chrm, i_pa_hq_pos=has_overlap(s_chrm, i_pos, m_pa_hq_tprt_both, 15)
            b_in_ma_hq, s_ma_hq_chrm, i_ma_hq_pos=has_overlap(s_chrm, i_pos, m_ma_hq_tprt_both, 15)
            ####
            b_in_pa_clip, s_pa_clip_chrm, i_pa_clip_pos=has_overlap(s_chrm, i_pos, m_pa_clip, 25)
            b_in_ma_clip, s_ma_clip_chrm, i_ma_clip_pos = has_overlap(s_chrm, i_pos, m_ma_clip, 25)

            ####whether it is qualified 0/1 in proband
            b_qualified=False
            f_clip_ratio = l_features[9]
            f_disc_ratio = l_features[10]
            b_VAF_in_range=(f_clip_ratio>=0.3 and f_clip_ratio<=0.85) and (f_disc_ratio>=0.3 and f_disc_ratio<=0.85)
            #b_VAF_in_range=True
            if b_in_pa_hq==True and b_in_ma_clip==False and b_VAF_in_range==True:
                b_qualified=True
            elif b_in_ma_hq==True and b_in_pa_clip==False and b_VAF_in_range==True:
                b_qualified=True
            ####
            ####

            if b_qualified==True:
                #l_slcted_01.append(l_features)
                if s_chrm not in m_slcted_01:
                    m_slcted_01[s_chrm]={}
                m_slcted_01[s_chrm][i_pos]=1####here just given a hard value
            
    return m_slcted_01
####

#load in the xTea identified "de novo" raw TE insertions
def select_fp_from_xtea_raw_denovo_insertions(sf_child, s_type, m_tp_manual_insp, i_slack):
    m_slcted_fp={}
    s_denovo_name="candidate_disc_filtered_cns_post_filtering_denovo.txt"
    sf_child_denovo=sf_child+s_type+"/"+s_denovo_name
    
    if os.path.isfile(sf_child_denovo)==False:
        print("%s failed to have raw de novo insertion results" % (sf_child_denovo))
        return {}

    with open(sf_child_denovo) as fin_raw_denovo:
        for line in fin_raw_denovo:
            fields=line.rstrip().split()
            s_chrm=fields[0]
            i_pos=int(fields[1])

            b_tp=False
            if s_chrm in m_tp_manual_insp:
                i_start=i_pos-i_slack
                i_end=i_pos+i_slack
                for i in range(i_start, i_end):
                    if i in m_tp_manual_insp[s_chrm]:
                        b_tp=True
                        break
            if b_tp==True:
                continue
            if s_chrm not in m_slcted_fp:
                m_slcted_fp[s_chrm]={}
            m_slcted_fp[s_chrm][i_pos]=1
    return m_slcted_fp


#load in the manual inspected insertions
def load_in_manual_inspected_TP_insertions(sf_tp):
    m_tp_manual_inspect={}
    with open(sf_tp) as fin:
        for line in fin:
            fields=line.rstrip().split()
            s_chrm=fields[0]
            i_pos=int(fields[1])
            if s_chrm not in m_tp_manual_inspect:
                m_tp_manual_inspect[s_chrm]={}
            m_tp_manual_inspect[s_chrm][i_pos]=1
    return m_tp_manual_inspect


####
def load_in_xtea_output(sf_rslts, s_flag=""):
    m_rslts={}
    with open(sf_rslts) as fin_rslts:
        for line in fin_rslts:
            fields = line.rstrip().split()##
            if len(s_flag)>1:
                if s_flag not in line:
                    continue

            s_chrm=fields[0]
            i_pos=int(fields[1])
            if s_chrm not in m_rslts:
                m_rslts[s_chrm]={}
            if len(s_flag)>1:
                m_rslts[s_chrm][i_pos]=1
            else:
                m_rslts[s_chrm][i_pos] = 1
    return m_rslts


def load_in_xtea_output_hq(sf_rslts, s_flag=""):
    m_rslts={}
    with open(sf_rslts) as fin_rslts:
        for line in fin_rslts:
            fields = line.rstrip().split()##
            if len(s_flag)>1:
                if s_flag not in line:
                    continue

            s_chrm=fields[0]
            i_pos=int(fields[1])
            if s_chrm not in m_rslts:
                m_rslts[s_chrm]={}
            if len(s_flag)>1:
                m_rslts[s_chrm][i_pos]=fields
            else:
                m_rslts[s_chrm][i_pos] = []
    return m_rslts

#
def load_in_features_one_insertion(l_fields):
    gntp_classifier = GntpClassifier_DF21()#
    l_features=gntp_classifier._parser_features(l_fields)
    return l_features

#
def has_overlap(s_chrm, i_pos, m_sites, i_slack=15):
    b_hit=False
    if s_chrm not in m_sites:
        return b_hit, "None", -1
    for i_tmp in range(i_pos-i_slack, i_pos+i_slack):
        if i_tmp in m_sites[s_chrm]:
            return True, s_chrm, i_tmp
    return False, "None", -1


#function to parse parameters
##parse the options
def parse_option():
    parser = OptionParser()

    parser.add_option("--pa", dest="pa_list", default="",
                        help="pa output folder list", metavar="FILE")
    parser.add_option("--ma", dest="ma_list", default="",
                        help="ma output folder list", metavar="FILE")
    parser.add_option("--proband", dest="proband_list", default="",
                        help="proband output folder list", metavar="FILE")
    parser.add_option("--bam", dest="bam_list", default="",
                        help="bam list", metavar="FILE")
    parser.add_option("--tp", dest="tp",
                      help="Manual inspected TPs", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-i", "--input", dest="input", default="",
                    help="input file ", metavar="FILE")
    parser.add_option("--type", dest="type", type="string", default="L1",
                      help="Repeat Type, e.g. L1, Alu, SVA")
    parser.add_option("-r", "--reference", dest="reference",
                      help="The reference file ", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string", default="./",
                      help="Working folder")
    parser.add_option("-e", "--extend", dest="extend", type="int", default=450,
                      help="extend length")
    (options, args) = parser.parse_args()
    return (options, args)

####
if __name__ == '__main__':
    (options, args) = parse_option()

    sf_child_list=options.proband_list
    sf_pa_list=options.pa_list
    sf_ma_list=options.ma_list
    sf_bam_list=options.bam_list
    sf_tp_manual_inspec=options.tp
    s_working_folder = options.wfolder
    sf_ref = options.reference
    i_extnd = options.extend
    s_type=options.type
    if len(s_working_folder) > 0 and s_working_folder[-1] != "/":
        s_working_folder += "/"
    sf_out_cmd = options.output

prep_cmds_to_generate_positive_training_images(sf_bam_list, sf_child_list, sf_pa_list, sf_ma_list, s_working_folder, s_type, sf_ref, i_extnd, sf_tp_manual_inspec, sf_out_cmd)

####