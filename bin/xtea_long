#!/usr/bin/env python

##11/04/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
from subprocess import *
from optparse import OptionParser
import ntpath
#import global_values

S_VERSION="v0.1"
LRD_PSEUDOGENE_INPUT_SUF=".for_pseudogene.fa"#

####
PUB_CLIP="pub_clip"
REP_TYPE_L1="L1"
REP_TYPE_ALU="Alu"
REP_TYPE_SVA="SVA"
REP_TYPE_HERV="HERV"
REP_TYPE_MIT="Mitochondrion"
REP_TYPE_MSTA="MSTA"

####
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

####
def gnrt_script_head_lsf(spartition, ncores, stime, smemory, s_id):
    s_head = "#!/bin/bash\n\n"
    s_head += "#BSUB -n {0}\n".format(ncores)
    if stime != None and stime != "":
        s_head += "#BSUB -W {0}\n".format(stime)
    s_head += "#BSUB -M {0}G\n".format(smemory)
    s_head += "#BSUB -q {0}\n".format(spartition)
    s_head += "#BSUB -o {0}_%J.out\n".format(s_id)
    s_head += "#BSUB -e {0}_%J.err\n".format(s_id)
    s_head += "#BSUB -R \"rusage[mem={0}G]\"\n".format(int(smemory/ncores))
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

# gnrt pars
def gnrt_parameters(l_pars):
    s_pars = ""
    for rcd in l_pars:
        sid = rcd[0]
        svalue = str(rcd[1])
        sline = sid + "=" + svalue + "\n"
        s_pars += sline
    return s_pars
####
####
####
def gnrt_calling_command(ncores, iflag, i_rep_type, sf_rmsk, sf_rep_cns, i_min_copy_len, i_peak_win, b_complex,
                         b_fast_mode=False, b_asm_cmd=False, b_call_seq=False, b_clean=False):####

    sclip_step = "python ${{XTEA_PATH}}\"l_main.py\" -C -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                 "-o ${{PREFIX}}\"candidate_list_from_clip.txt\"  -n {0} -w {1} \n".format(ncores, i_peak_win)

    s_fast_mode="--fast"
    if b_fast_mode==False:
        s_fast_mode=""
    s_asm_step = "python ${{XTEA_PATH}}\"l_main.py\" -A -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                 "-i ${{PREFIX}}\"candidate_list_from_clip.txt\" -o ${{PREFIX}}\"all_ins_seqs.fa\"" \
                 " --rep ${{REP_LIB}} -n {0} {1}\n".format(ncores, s_fast_mode)
    if b_complex==True:
        s_asm_step = "python ${{XTEA_PATH}}\"l_main.py\" -A --collect_asm -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                     "-i ${{PREFIX}}\"candidate_list_from_clip.txt.left_breakponts\" -o ${{PREFIX}}\"all_ins_seqs_left.fa\"" \
                     " --rep ${{REP_LIB}} -n {0} {1}\n".format(ncores, s_fast_mode)
        s_asm_step += "python ${{XTEA_PATH}}\"l_main.py\" -A --collect_asm -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                     "-i ${{PREFIX}}\"candidate_list_from_clip.txt.right_breakponts\" -o ${{PREFIX}}\"all_ins_seqs_right.fa\"" \
                     " --rep ${{REP_LIB}} -n {0} {1}\n".format(ncores, s_fast_mode)
    if b_asm_cmd==True and b_call_seq==False:
        s_asm_step = "python ${{XTEA_PATH}}\"l_main.py\" -A --cmd -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                     "-i ${{PREFIX}}\"candidate_list_from_clip.txt\" -o ${{PREFIX}}\"all_ins_seqs.fa\"" \
                     " --rep ${{REP_LIB}} -n {0} {1}\n".format(ncores, s_fast_mode)
        if b_complex==True:
            s_asm_step = "python ${{XTEA_PATH}}\"l_main.py\" -A --cmd --collect_asm -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                         "-i ${{PREFIX}}\"candidate_list_from_clip.txt.left_breakponts\" -o ${{PREFIX}}\"all_ins_seqs_left.fa\"" \
                         " --rep ${{REP_LIB}} -n {0} {1}\n".format(ncores, s_fast_mode)
            s_asm_step += "python ${{XTEA_PATH}}\"l_main.py\" -A --cmd --collect_asm -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                         "-i ${{PREFIX}}\"candidate_list_from_clip.txt.right_breakponts\" -o ${{PREFIX}}\"all_ins_seqs_right.fa\"" \
                         " --rep ${{REP_LIB}} -n {0} {1}\n".format(ncores, s_fast_mode)
    elif b_asm_cmd==False and b_call_seq==True:
        s_asm_step = "python ${{XTEA_PATH}}\"l_main.py\" -A --mei_no_asm -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                     "-i ${{PREFIX}}\"candidate_list_from_clip.txt\" -o ${{PREFIX}}\"all_ins_seqs.fa\"" \
                     " --rep ${{REP_LIB}} -n {0} {1}\n".format(ncores, s_fast_mode)

    s_ghost_step="python ${{XTEA_PATH}}\"l_main.py\" -N -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}}\"ghost\" " \
                 "-o ${{PREFIX}}\"ghost_reads.fa\" --rmsk {0} --cns {1} --min {2}" \
               " -n {3}\n".format(sf_rmsk, sf_rep_cns, i_min_copy_len, ncores)
    s_classify_step="python ${{XTEA_PATH}}\"l_main.py\" -Y -i ${{PREFIX}}\"all_ins_seqs.fa\" -r ${{REF}} " \
                    "-p ${{TMP}}\"classification\" --rep ${{REP_LIB}} -y {0} " \
                    "-o ${{PREFIX}}\"classified_results.txt\" -n {1}\n".format(i_rep_type, ncores)
    if b_complex==True:#complex events
        s_classify_step+="python ${{XTEA_PATH}}\"l_main.py\" -X -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                         "-i ${{PREFIX}}\"candidate_list_from_clip.txt.left_breakponts\" " \
                         "--i2 ${{PREFIX}}\"candidate_list_from_clip.txt.right_breakponts\" " \
                         " --iline ${{PREFIX}}\"classified_results.txt.merged_LINE1.txt\" " \
                         "--isva ${{PREFIX}}\"classified_results.txt.merged_SVA.txt\" --rep ${{REP_LIB}} " \
                         "-o ${{PREFIX}}\"complex_TE_SV\" -n {1}\n".format(i_rep_type, ncores)

    s_clean="python ${{XTEA_PATH}}\"l_main.py\" --clean -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                     "-i ${{PREFIX}}\"candidate_list_from_clip.txt\" " \
                     " -n {0}\n".format(ncores)

    sf_ins_seq_for_pseudo="all_ins_seqs.fa"+LRD_PSEUDOGENE_INPUT_SUF
    s_pseudo="python ${{XTEA_PATH}}\"l_main.py\" --pseudo -i ${{PREFIX}}\"{0}\" " \
             "-r ${{REP_LIB}}\"gene_annotation/exon.fa\" " \
             "--gene ${{REP_LIB}}\"gene_annotation/gencode.v28.GRCh38.annotation.gff3\" " \
             "-p ${{TMP}}\"classification\" " \
             "-o ${{PREFIX}}\"pseudogene_insertion_results.txt\" -n {1}\n".format(sf_ins_seq_for_pseudo, ncores)

    s_dimorphic="python ${{XTEA_PATH}}\"l_dimorphic_HERV.py\"  -p ${{TMP}} " \
                "-i ${{PREFIX}}\"candidate_list_from_clip.txt.left_breakponts\" " \
                "--i2 ${{PREFIX}}\"candidate_list_from_clip.txt.right_breakponts\" " \
                "-a ${{REP_LIB}}\"LTR/hg38_LTR.fa.out\" " \
                "-o ${{PREFIX}}\"dimorphic_herv_candidates.txt\" -n {0}\n".format(ncores)

    s_ref_sva="python ${{XTEA_PATH}}\"l_main.py\" --rsva -b ${{BAM_LIST}} -r ${{REF}} -p ${{TMP}} " \
                     "-i ${{SVA_REF_COPY}} -o ${{PREFIX}}\"ref_SVA_copy_seqs.fa\"" \
                     " -n {0}\n".format(ncores)

    s_cmd = ""
    if iflag & 1 == 1:
        s_cmd += sclip_step
    if iflag & 2 == 2:
        s_cmd += s_asm_step
    if iflag & 4 ==4:
        s_cmd += s_ghost_step
    if iflag & 8 !=0:
        s_cmd += s_classify_step
    if iflag & 16 !=0:
        s_cmd += s_clean
    if iflag & 32 !=0:
        s_cmd += s_pseudo
    if iflag & 64 != 0:
        s_cmd += s_dimorphic
    if iflag & 128 !=0:
        s_cmd += s_ref_sva
    return s_cmd
####

####gnrt the whole pipeline
def gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_id, sf_bams, sf_root_folder):
    l_sbatch_files=[]
    sf_working_folder=sf_root_folder
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
            sf_tmp_ghost = sf_folder + "/tmp/ghost"
            if os.path.exists(sf_tmp_ghost) == False:
                cmd = "mkdir {0}".format(sf_tmp_ghost)
                run_cmd(cmd)
            sf_tmp_classify = sf_folder + "/tmp/classification"
            if os.path.exists(sf_tmp_classify) == False:
                cmd = "mkdir {0}".format(sf_tmp_classify)
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
                    fout_bam_list.write(sf_tmp_bam + "\n")

        ####gnrt the pipeline file
        sf_out_sh = sf_folder + "run_xTEA_pipeline.sh"
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


def write_to_config(sf_ref, sf_xtea, s_bl, s_tmp, sf_rep_folder, sf_ref_sva, sf_config):
    with open(sf_config, "w") as fout_L1:
        fout_L1.write(sf_ref)
        fout_L1.write(sf_xtea)
        fout_L1.write(s_bl)
        fout_L1.write(s_tmp)
        fout_L1.write(sf_rep_folder)
        fout_L1.write(sf_ref_sva)

####
#generate library configuration files
def gnrt_lib_config(sf_ref, sf_folder_xtea, sf_config_prefix, sf_rep_folder, sf_ref_sva, sf_config):
    if sf_folder_xtea[-1] != "/":
        sf_folder_xtea += "/"
    if sf_config_prefix[-1] != "/":
        sf_config_prefix += "/"
    s_bl = "BAM_LIST ${PREFIX}\"bam_list.txt\"\n"
    s_tmp = "TMP ${PREFIX}\"tmp/\"\n"
    sf_ref = "REF " + sf_ref + "\n"
    sf_xtea = "XTEA_PATH " + sf_folder_xtea + "\n"
    sf_rep_folder="REP_LIB "+ sf_rep_folder +"\n"
    sf_ref_sva="SVA_REF_COPY "+sf_ref_sva+"\n"
    write_to_config(sf_ref, sf_xtea, s_bl, s_tmp, sf_rep_folder, sf_ref_sva, sf_config)

#
def cp_file(sf_from, sf_to):
    cmd = "cp {0} {1}".format(sf_from, sf_to)
    if os.path.isfile(sf_from)==False:
        return
    run_cmd(cmd)

def run_cmd(cmd):
    print(cmd)
    Popen(cmd, shell=True, stdout=PIPE).communicate()

def get_sample_id(sf_bam):
    fname = ntpath.basename(sf_bam)
    fname_fields = fname.split(".")
    if fname_fields[-1] != "bam" and fname_fields[-1] != "cram":
        print("Alignment is not end with .bam")
        return None
    sample_id = ".".join(fname_fields[:-1])
    return sample_id

####gnrt the running shell
def gnrt_running_shell(sf_ids, sf_bams, s_wfolder, sf_ref, sf_folder_xtea, sf_rep_folder, spartition, stime, smemory,
                       ncores, sf_rmsk, sf_rep_cns, sf_ref_sva, i_min_copy_len, i_peak_win, sf_submit_sh, b_complex,
                       b_fast_mode=False, b_lsf=False, b_slurm=False, b_asm_cmd=False, b_call_seq=False, b_clean=False):
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
    split_bam_list(m_id, sf_bams, s_wfolder)

    l_sh=[]
    for sid_tmp in m_id:
        sf_sample_folder=s_wfolder + sid_tmp + "/"
        sf_config = sf_sample_folder + "xtea.config"
        gnrt_lib_config(sf_ref, sf_folder_xtea, sf_sample_folder, sf_rep_folder, sf_ref_sva, sf_config)

        sf_rep_sample_id = sf_sample_folder + "/sample_id.txt"
        sf_rep_bam = sf_sample_folder + "/bam_list1.txt"
        if os.path.isfile(sf_rep_bam)==False:
            sf_rep_bam="null"

        s_head = ""
        if b_lsf==True:
            s_head = gnrt_script_head_lsf(spartition, ncores, stime, smemory, sid_tmp)
        elif b_slurm==True:
            s_head = gnrt_script_head(spartition, ncores, stime, smemory, sid_tmp)

        l_libs = load_par_config(sf_config)
        s_libs = gnrt_parameters(l_libs)
        iflag = options.flag #which step to run
        i_rep_type=options.rep_type
        s_calling_cmd = gnrt_calling_command( ncores, iflag, i_rep_type, sf_rmsk, sf_rep_cns,
                                              i_min_copy_len, i_peak_win, b_complex, b_fast_mode,
                                              b_asm_cmd,  b_call_seq, b_clean)

        l_tmp_sh=gnrt_pipelines(s_head, s_libs, s_calling_cmd, sf_rep_sample_id, sf_rep_bam,
                       sf_sample_folder)
        for tmp_sh in l_tmp_sh:
            l_sh.append(tmp_sh)
    with open(sf_submit_sh, "w") as fout_submit:
        fout_submit.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            if b_lsf==False:
                fout_submit.write("sbatch < "+s_sh+"\n")
            else:
                fout_submit.write("bsub < " + s_sh + "\n")
####

####Input:
# m_ids: sample id dictionary
def split_bam_list(m_ids, sf_bams, s_wfolder):
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

    for sid_tmp in m_ids:
        sf_sample_folder = s_wfolder + sid_tmp + "/"
        if os.path.exists(sf_sample_folder) == False:
            cmd = "mkdir {0}".format(sf_sample_folder)
            run_cmd(cmd)

        # need to generate two files for each sample
        sf_rep_sample_id = sf_sample_folder + "/sample_id.txt"
        with open(sf_rep_sample_id, "w") as fout_rep_sample_id:
            fout_rep_sample_id.write(sid_tmp)
        if sid_tmp in m_bams:
            sf_rep_bam = sf_sample_folder + "/bam_list1.txt"
            with open(sf_rep_bam, "w") as fout_rep_bams:
                for sf_tmp_bam in m_bams[sid_tmp]:
                    fout_rep_bams.write(sid_tmp+"\t"+sf_tmp_bam+"\n")
#
####run the pipelines
def run_pipeline(l_rep_type, sample_id, s_wfolder):
    for rep_type in l_rep_type:
        sf_sbatch_sh_rep = s_wfolder + sample_id + "/run_{0}.sh".format(rep_type)
        cmd="sh {0}".format(sf_sbatch_sh_rep)
        run_cmd(cmd)

####
def parse_option():
    parser = OptionParser()
    parser.add_option("-i", "--id", dest="id",
                      help="sample id list file ", metavar="FILE")
    parser.add_option("-a", "--par", dest="parameters",
                      help="parameter file ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")

    parser.add_option("-n", "--cores", dest="cores", type="int", default=8,
                      help="number of cores")
    parser.add_option("-m", "--memory", dest="memory", type="int", default=32,
                      help="Memory limit in GB")
    parser.add_option("-q", "--partition", dest="partition", type="string", default="long",
                      help="Which queue to run the job")
    parser.add_option("-t", "--time", dest="time", type="string", default="120:00",
                      help="Time limit")
    parser.add_option("--lsf",
                      action="store_true", dest="lsf", default=False,
                      help="Indiates submit to LSF system")
    parser.add_option("--slurm",
                      action="store_true", dest="slurm", default=False,
                      help="Indiates submit to slurm system")
    parser.add_option("--cmd",
                      action="store_true", dest="asm_cmd", default=False,
                      help="Generate asm command script (for cluster)")
    parser.add_option("--complex",
                      action="store_true", dest="complex", default=False,
                      help="Call complex events (TE promoted SV)")
    parser.add_option("--mei_no_asm",
                      action="store_true", dest="mei_no_asm", default=False,
                      help="Call MEI only without asm")
    parser.add_option("--clean",
                      action="store_true", dest="clean", default=False,
                      help="Clean the intermediate files")#
    parser.add_option("-V", "--version",
                      action="store_true", dest="version", default=False,
                      help="Print xTea version")
    parser.add_option("--fast",
                      action="store_true", dest="fast", default=False,
                      help="This is the fast mode, which may sacrifice the sensitivity")


    parser.add_option("-p", "--path", dest="wfolder", type="string", default="./",
                      help="Working folder")
    parser.add_option("-r", "--ref", dest="ref", type="string",
                      help="reference genome")
    parser.add_option("-g", "--gene", dest="gene", type="string",
                      help="Gene annotation file")
    parser.add_option("-w", "--win", dest="win", type="int", default=75,
                      help="peak window size")
    parser.add_option("--xtea", dest="xtea", type="string",
                      help="xTEA folder")
    parser.add_option("--rep", dest="rep_lib", type="string",
                      help="Repeat library folder")
    parser.add_option("--rmsk", dest="rmsk", type="string", default="null",
                      help="RepeatMasker output on the reference")
    parser.add_option("--cns", dest="consensus", default="null",
                      help="repeat consensus folder", metavar="FILE")
    parser.add_option("--ref_sva", dest="ref_sva", default="null",
                      help="reference SVA copies", metavar="FILE")
    parser.add_option("--min", dest="min_copy_len", type="int", default=4000,
                      help="Minimum copy length for find polymorhpic copies")

    parser.add_option("-f", "--flag", dest="flag", type="int", default=31,
                      help="Flag indicates which step to run (1-clip, 2-asm, 4-ghost, 8-classification, 16-clean)")
    parser.add_option("-y", "--type", dest="rep_type", type="int", default=7,
                      help="Type of repeats working on (1-LINE1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondria)")
    parser.add_option("-o", "--output", dest="output", default="submit_jobs_long_reads.sh",
                      help="The output file", metavar="FILE")
    (options, args) = parser.parse_args()
    return (options, args)

####
if __name__ == '__main__':
    (options, args) = parse_option()
    b_version = options.version
    if b_version==True:#
        print(("xTea %s for long reads on hg38\n" % S_VERSION))
    else:
        sf_id = options.id
        sf_bams = options.bam ###input is a bam file
        sf_sbatch_sh = options.output  # this is the shell for submitting the jobs
        s_wfolder = options.wfolder
        b_lsf = options.lsf
        b_slurm=options.slurm
        b_asm_cmd = options.asm_cmd
        b_complex=options.complex #for complex events
        b_call_seq = options.mei_no_asm
        b_clean=options.clean
        i_peak_win = options.win #maximum window size when try to cluster the breakpoints

        if s_wfolder[-1]!="/":
            s_wfolder+="/"
        if os.path.exists(s_wfolder) == False:
            scmd = "mkdir {0}".format(s_wfolder)
            Popen(scmd, shell=True, stdout=PIPE).communicate()

        if os.path.isfile(sf_bams) == False:
            sf_bams = "null"

    ####
        spartition = options.partition
        stime = options.time
        smemory = options.memory
        ncores = options.cores
        sf_ref=options.ref ####reference genome
        sf_folder_xtea=options.xtea

    ####
        sf_rmsk = options.rmsk
        i_min_copy_len=options.min_copy_len
        sf_rep_cns = options.consensus  # repeat consensus
        sf_rep_folder=options.rep_lib#repeat library folder
        sf_ref_sva=options.ref_sva
        b_fast_mode=options.fast
    ####
        gnrt_running_shell(sf_id, sf_bams, s_wfolder, sf_ref, sf_folder_xtea, sf_rep_folder, spartition, stime, smemory,
                           ncores, sf_rmsk, sf_rep_cns, sf_ref_sva, i_min_copy_len, i_peak_win, sf_sbatch_sh, b_complex,
                           b_fast_mode, b_lsf, b_slurm, b_asm_cmd, b_call_seq, b_clean)
    ####