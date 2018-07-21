##06/12/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
####this is the stand alone filtering version

import os
import sys
from optparse import OptionParser
from x_clip_disc_filter import *


def parse_option():
    parser = OptionParser()
    parser.add_option("-R", "--realign",
                      action="store_true", dest="realign", default=False,
                      help="Collect and realign the clip and disc reads")
    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("-c", "--cns", dest="cns",
                      help="The consensus file ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("--ref", dest="reference",
                      help="Reference genome", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")
    parser.add_option("-e", "--extend", dest="extend", type="int",
                      help="extend length")
    parser.add_option("--copies", dest="copies",
                      help="repeat copies", metavar="FILE")
    parser.add_option("--fflank", dest="fflank",
                      help="flank region file", metavar="FILE")
    parser.add_option("--flklen", dest="flklen", type="int",
                      help="flank region file")
    parser.add_option("--nclip", dest="nclip", type="int",
                      help="cutoff of minimum # of clipped parts fall in repeats")
    parser.add_option("--ndisc", dest="ndisc", type="int",
                      help="cutoff of minimum # of discordant pair")
    parser.add_option("--is", dest="insz", type="int",
                      help="Mean insert size +/- 3*derivation")
    (options, args) = parser.parse_args()
    return (options, args)

####
def get_candidate_list_from_TEA_mem_output(sf_input, sf_output):
    m_old_info={}
    with open(sf_input) as fin_list, open(sf_output, "w") as fout_list:
        for line in fin_list:
            fields=line.split()
            if fields[0]=="sample":
                continue
            chrm=fields[1]
            clip_pos=fields[7]#by default, use the right clip position
            nlclip=int(fields[15])
            nrclip=int(fields[16])
            if nlclip>=nrclip:
                clip_pos=fields[6]
            fout_list.write(chrm+"\t"+clip_pos+"\n")
            if chrm not in m_old_info:
                m_old_info[chrm]={}
            m_old_info[chrm][int(clip_pos)]=line.rstrip()
    return m_old_info

def dump_old_info(sf_old_pos, m_old_info, sf_out):
    with open(sf_out,"w") as fout_info:
        with open(sf_old_pos) as fin_pos:
            for line in fin_pos:
                fields=line.split()
                ins_chrm=fields[0]
                pos=int(fields[1])
                if (ins_chrm in m_old_info) and (pos in m_old_info[ins_chrm]):
                    s_info=m_old_info[ins_chrm][pos]
                    fout_info.write(s_info+"\n")

####
if __name__ == '__main__':
    (options, args) = parse_option()

    sf_bam_list = options.bam  ###read in a bam list file
    s_working_folder = options.wfolder
    n_jobs = options.cores
    sf_reference = options.reference

    sf_candidate_list = options.input
    sf_rep_cns = options.cns  ####repeat copies here, with the flank regions
    sf_output = options.output

    iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
    bin_size = 50000000  # block size for parallelization
    bmapped_cutoff = 0.5
    f_concord_ratio = 0.45
    i_concord_dist = options.insz  # this should be the 3*std_derivation
    sf_rep_copies = options.copies  ####repeat copies here, with the flank regions
    sf_flank = options.fflank
    i_flank_lenth = options.flklen
    n_clip_cutoff = options.nclip# this is the sum of left and right clipped reads
    n_disc_cutoff = options.ndisc #each sample should have at least this number of discordant reads

    sf_candidate_list_brkpnt=sf_candidate_list+".brkpnt"
    m_old_info=get_candidate_list_from_TEA_mem_output(sf_candidate_list, sf_candidate_list_brkpnt)

    if os.path.exists(s_working_folder)==False:
        print "Working folder {0} doesn't exist!!!".format(s_working_folder)
        sys.exit()

    x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_reference)
    x_cd_filter.call_MEIs_consensus(sf_candidate_list_brkpnt, iextnd, bin_size, sf_rep_copies, sf_flank, i_flank_lenth,
                                    bmapped_cutoff, i_concord_dist, f_concord_ratio, n_clip_cutoff, n_disc_cutoff,
                                    sf_output)
    sf_old_pos=sf_output+".old_positions"
    sf_tea_mem=sf_output+".tea"
    dump_old_info(sf_old_pos, m_old_info, sf_tea_mem)
####