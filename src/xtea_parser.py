#!/usr/bin/env python

##04/26/2024
##@@author: Corinne Sexton, DBMI, Harvard Medical School
##@@contact: corinne_sexton@hms.harvard.edu
##


import os
from subprocess
import argparse
import ntpath


def make_parser():
    desc = "Detecting TE insertions in short and long-read data."

    #TODO add example
    example = "Example: run_xtea ..."

    parser = argparse.ArgumentParser(prog='run_xtea', description=desc, epilog=example)

    parser.add_argument("-C", "--case_control",
                        action="store_true", dest="case_control", default=False,
                        help="Run in case control mode")
    parser.add_argument("-M", "--mosaic",
                        action="store_true", dest="mosaic", default=False,
                        help="Calling mosaic events from high coverage data")
    parser.add_argument("--denovo",
                        action="store_true", dest="denovo", default=False,
                        help="Run in de novo mode")
    parser.add_argument("-U", "--user",
                        action="store_true", dest="user", default=False,
                        help="Use user specific parameters instead of automatically calculated ones")
    parser.add_argument("--force",
                        action="store_true", dest="force", default=False,
                        help="Force to start from the very beginning")
    parser.add_argument("--hard",
                         action="store_true", dest="hard", default=False,
                        help="This is hard-cut for fitering out coverage abnormal candidates")
    parser.add_argument("--tumor",
                        action="store_true", dest="tumor", default=False,
                        help="Working on tumor samples")
    parser.add_argument("--purity", dest="purity", type="float", default=0.45,  # by default tumor purity set to 45%
                        help="Tumor purity")
    parser.add_argument("--lsf",
                        action="store_true", dest="lsf", default=False,
                        help="Indiates submit to LSF system")
    parser.add_argument("--slurm",
                        action="store_true", dest="slurm", default=False,
                        help="Indiates submit to slurm system")
    parser.add_argument("--resume",
                        action="store_true", dest="resume", default=False,
                        help="Resume the running, which will skip the step if output file already exists!")
    parser.add_argument("-v", "--version",
                        action="store_true", dest="version", default=False,
                        help="Print xTea version")
    
    parser.add_argument("-i", "--id", dest="id",
                      help="sample id list file ", metavar="FILE")
    parser.add_argument("-a", "--par", dest="parameters",
                      help="parameter file ", metavar="FILE")
    parser.add_argument("-l", "--lib", dest="lib",
                      help="TE lib config file ", metavar="FILE")
    parser.add_argument("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_argument("-x", "--x10", dest="x10",
                      help="Input 10X bam file", metavar="FILE")
    
    parser.add_argument("-n", "--cores", dest="cores", type="int", default=1,
                      help="number of cores")
    parser.add_argument("-m", "--memory", dest="memory", type="int", default=16,
                      help="Memory limit in GB")
    parser.add_argument("-q", "--partition", dest="partition", type="string",
                      help="Which queue to run the job")
    parser.add_argument("-t", "--time", dest="time", type="string", default="",
                      help="Time limit")
    
    parser.add_argument("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_argument("-r", "--ref", dest="ref", type="string",
                      help="reference genome")
    parser.add_argument("-g", "--gene", dest="gene", type="string",default="gene.aff3",
                      help="Gene annotation file")
    parser.add_argument("--xtea", dest="xtea", type="string",
                      help="xTEA folder")


    parser.add_argument("-f", "--flag", dest="flag", type="int",
                      help="Flag indicates which step to run (1-clip, 2-disc, 4-barcode, 8-xfilter, 16-filter, 32-asm)")

    parser.add_argument("-y", "--reptype", dest="rep_type", type="int", default=1,
                      help="Type of repeats working on: 1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrial")

    parser.add_argument("--flklen", dest="flklen", type="int", default=3000,
                      help="flank region file")
    parser.add_argument("--nclip", dest="nclip", type="int", default=3,
                      help="cutoff of minimum # of clipped reads")
    parser.add_argument("--cr", dest="cliprep", type="int", default=1,
                      help="cutoff of minimum # of clipped reads whose mates map in repetitive regions")
    parser.add_argument("--nd", dest="ndisc", type="int", default=5,
                      help="cutoff of minimum # of discordant pair")
    parser.add_argument("--nfclip", dest="nfilterclip", type="int", default=3,
                      help="cutoff of minimum # of clipped reads in filtering step")
    parser.add_argument("--nfdisc", dest="nfilterdisc", type="int", default=5,
                      help="cutoff of minimum # of discordant pair of each sample in filtering step")
    parser.add_argument("--teilen", dest="teilen", type="int", default=50,
                      help="minimum length of the insertion for future analysis")

    parser.add_argument("-o", "--output", dest="output", default="submit_calling_jobs_for_samples.sh",
                      help="The output file", metavar="FILE")
    parser.add_argument("--blacklist", dest="blacklist", default="null",
                      help="Reference panel database for filtering, or a blacklist region", metavar="FILE")
    (options, args) = parser.parse_args()

    return parser



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    b_version=options.version
    if b_version==True:
        print(("xTea %s for short and linked reads on hg38\n" % S_VERSION))
    else:
        sf_id = options.id
        sf_bams = options.bam ###input is a bam file
        sf_bams_10X = options.x10
        sf_sbatch_sh = options.output  # this is the shell for submitting the jobs
        s_wfolder = options.wfolder
        b_mosaic = options.mosaic
        b_user_par = options.user
        b_force = options.force
        b_case_control=options.case_control #whether case control mode
        b_tumor = options.tumor
        f_purity = options.purity
        b_hard_cut = options.hard
        b_lsf=options.lsf
        b_slurm=options.slurm ####
        b_resume = options.resume
        b_denovo=options.denovo
    ####

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
        sf_folder_xtea=options.xtea#

        sf_folder_rep=sf_folder_rep1
        sf_ref=sf_ref1
        if options.decompress==True:
            decompress(sf_folder_rep1, s_wfolder)
            decompress(sf_ref1, s_wfolder)
            sf_folder_rep = s_wfolder+"rep_lib_annotation/" #trim tar.gz
            sf_ref=s_wfolder+"genome.fa"
        sf_gene=options.gene
        sf_black_list = options.blacklist

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
        if i_rep_type & 64 != 0:
            l_rep_type.append(REP_TYPE_PSEUDOGENE)

        if b_case_control==False:
            #non-case-ctrl mode (mainly for germline)
            gnrt_running_shell(sf_id, sf_bams, sf_bams_10X, l_rep_type, b_mosaic, b_user_par, b_force, b_tumor, f_purity, b_hard_cut,
                               s_wfolder, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea, spartition, stime,
                               smemory, ncores, sf_sbatch_sh, "null", b_lsf, b_slurm, b_resume)

        elif b_denovo==False:
            #if in case-control mode, then need to separate to two steps
            #thus need to prepare the files separately
            sf_sprt_bams=sf_bams + CASE_LIST_SUFFIX
            sf_control_bams=sf_bams + CTRL_LIST_SUFFIX
            prepare_case_control_bam(sf_bams, sf_sprt_bams, sf_control_bams)
            #first, run the jobs seperately for all the case and control samples
            gnrt_running_shell(sf_id, sf_sprt_bams, sf_bams_10X, l_rep_type, b_mosaic, b_user_par, b_force, b_tumor,
                               f_purity, b_hard_cut, s_wfolder, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea,
                               spartition, stime, smemory, ncores, sf_sbatch_sh, sf_bams, b_lsf, b_slurm, b_resume)
        else:
            # if in de novo mode, then need to separate to two steps
            # thus need to prepare the files separately
            sf_sprt_bams = sf_bams + CASE_LIST_SUFFIX
            sf_control_bams = sf_bams + CTRL_LIST_SUFFIX
            prepare_case_control_bam(sf_bams, sf_sprt_bams, sf_control_bams)
            # first, run the jobs seperately for all the case and control samples
            gnrt_running_shell(sf_id, sf_bams, sf_bams_10X, l_rep_type, b_mosaic, b_user_par, b_force, b_tumor,
                               f_purity, b_hard_cut, s_wfolder, sf_folder_rep, sf_ref, sf_gene, sf_black_list, sf_folder_xtea,
                               spartition, stime, smemory, ncores, sf_sbatch_sh, sf_bams, b_lsf, b_slurm, b_resume, b_denovo)

####
####
####