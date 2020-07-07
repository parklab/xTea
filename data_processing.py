####
##04/16/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

#Data processing:
# 1. given an dbgap sra id, download the data using prefetch and aspc
# 2. parse bam from the sra file
# 3. index the bam
# 4. align the given fastq list by sample

import os
import sys
import ntpath
from subprocess import *
from multiprocessing import Pool
from optparse import OptionParser


def run_cvt(cmd):
    Popen(cmd, shell=True, stdout=PIPE).communicate()

def parse_bam(sf_list, n_jobs, s_wfolder):
    if s_wfolder[-1]!="/":
        s_wfolder+="/"
    with open(sf_list) as fin_list:
        l_cmd=[]
        for line in fin_list:
            fields_ori=line.split('/')
            fields=fields_ori[-1].split(".")
            sf_sra=line.rstrip()
            sf_out=s_wfolder+fields[0]+".bam"
            cmd="sam-dump {0} | samtools view -bS -o {1} - ".format(sf_sra, sf_out)
            print cmd
            l_cmd.append(cmd)

        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()

def parse_bam_from_sra2(sf_runinfo_table, n_jobs):
    with open(sf_runinfo_table) as fin_list:
        l_cmd = []
        for line in fin_list:
            fields = line.split()
            sid = fields[13]
            if sid == "SRA_Sample":
                continue
            #sf_sra = line.rstrip()
            sf_sra=sid+".sra"
            if os.path.isfile(sf_sra)==False:
                print "Cannot find {0}\n".format(sf_sra)
                continue
            sf_out = fields[0] + ".bam"
            cmd = "sam-dump {0} | samtools view -bS - > {1}".format(sf_sra, sf_out)
            l_cmd.append(cmd)

        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()

def parse_fastq_from_sra(sf_list, n_jobs):
    with open(sf_list) as fin_list:
        l_cmd=[]
        for line in fin_list:
            fields=line.split()
            sid=fields[0]
            s_prefix="/n/data1/hms/dbmi/park/simon_chu/projects/data/dbGaP/sra/"
            sf_sra=s_prefix+sid+".sra"
            cmd="fastq-dump --split-files {0}".format(sf_sra)
            l_cmd.append(cmd)

        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()

####
#By default, the sf_runinfo_table is the directly downloaded file from dbGaP
def download_sra(sf_runinfo_table, n_jobs):
    with open(sf_runinfo_table) as fin_list:
        l_cmd=[]
        for line in fin_list:
            fields=line.split()
            sid=fields[8]
            if sid=="SRA_Sample":
                continue
            cmd="prefetch --max-size 500000000 -t ascp -a \"ascp|asperaweb_id_dsa.openssh\" " \
                "--ascp-options \"-k1 -Tr \" {0}".format(sid)
            l_cmd.append(cmd)

            #print cmd

        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()

def download_sra_with_list(sf_list, n_jobs):
    with open(sf_list) as fin_list:
        l_cmd=[]
        for line in fin_list:
            fields=line.split()
            sid=fields[4]
            if sid=="SRA_Sample" or sid=="Run":
                continue
            cmd="prefetch --max-size 500000000 -t ascp -a \"ascp|asperaweb_id_dsa.openssh\" " \
                "--ascp-options \"-k1 -Tr \" {0}".format(sid)
            l_cmd.append(cmd)

        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()

def realign_from_bam(sf_bam, sf_ref, i_mem, sf_out, sf_wfolder, n_jobs):
    sf_tmp=sf_out+"_tmp.bam"
    cmd = "samtools fastq {0} | bwa mem -t {1} {2} - | samtools view -hSb - " \
          "> {3}".format(sf_bam, n_jobs, sf_ref, sf_tmp)
    #run_cvt(cmd)
    cmd2="sambamba sort -m {0}G -o {1} --tmpdir={2} -t {3} {4}".format(i_mem, sf_out, sf_wfolder, n_jobs, sf_tmp)
    #run_cvt(cmd2)
    #remove the tmp files
    cmd3="rm {0}".format(sf_tmp)

    s_all_cmd=cmd+"\n"
    s_all_cmd+=(cmd2+"\n")
    s_all_cmd+=(cmd3+"\n")
    return s_all_cmd

####
def align_from_fq(rcd_fq, sf_ref, i_mem, sf_out_bam, sf_wfolder, n_jobs):
    sf_lfq=rcd_fq[0]
    sf_rfq=rcd_fq[1]
    sf_tmp = sf_out_bam + "_tmp.bam"
    cmd = "bwa mem -t {0} {1} {2} {3} | samtools view -hSb - " \
          "> {4}".format(n_jobs, sf_ref, sf_lfq, sf_rfq, sf_tmp)
    cmd2 = "sambamba sort -m {0}G -o {1} --tmpdir={2} -t {3} {4}".format(i_mem, sf_out_bam, sf_wfolder, n_jobs, sf_tmp)
    # remove the tmp files
    cmd3 = "rm {0}".format(sf_tmp)

    s_all_cmd = cmd + "\n"
    s_all_cmd += (cmd2 + "\n")
    s_all_cmd += (cmd3 + "\n")
    return s_all_cmd
####

####
def gnrt_script_head(spartition, ncores, stime, imemory):
    s_head = "#!/bin/bash\n\n"
    s_head += "#SBATCH -n {0}\n".format(ncores)
    s_head += "#SBATCH -t {0}\n".format(stime)
    s_head += "#SBATCH --mem={0}G\n".format(imemory)
    s_head += "#SBATCH -p {0}\n".format(spartition)
    s_head += "#SBATCH -o hostname_%j.out\n"
    if spartition == "park" or spartition == "priopark":
        s_head += "#SBATCH --account=park_contrib\n\n"
    return s_head

def parse_id_from_path(sf_bam):
    s_name = ntpath.basename(sf_bam)
    s_prefix = os.path.splitext(s_name)
    return s_prefix[0]

####
def realign_from_bam_list(sf_bam_list, sf_ref, spartition, i_mem, stime, n_cores, sf_wfolder, sf_run_sh):
    l_sh=[]
    sf_cmd_folder=sf_wfolder+"cmd/"
    if os.path.exists(sf_cmd_folder)==False:
        cmd="mkdir {0}".format(sf_cmd_folder)
        run_cvt(cmd)

    i_request_mem=i_mem+3
    i_request_cores=n_cores+2

    s_head=gnrt_script_head(spartition, i_request_cores, stime, i_request_mem)

    with open(sf_bam_list) as fin_list:
        for sline in fin_list:
            sf_bam=sline.rstrip()
            s_id = parse_id_from_path(sf_bam)
            sf_out_bam = sf_wfolder + s_id + ".re_align.sorted.bam"
            s_cmd = realign_from_bam(sf_bam, sf_ref, i_mem, sf_out_bam, sf_wfolder, n_cores)


            s_sh=sf_cmd_folder+"algn_{0}.sh".format(s_id)
            with open(s_sh,"w") as fout_tmp_sh:
                fout_tmp_sh.write(s_head)
                fout_tmp_sh.write(s_cmd)
            l_sh.append(s_sh)

    with open(sf_run_sh, "w") as fout_sh:
        fout_sh.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            fout_sh.write("sbatch < " + s_sh + "\n")

def index_check_status(sf_list, n_jobs):
    with open(sf_list) as fin_list:
        l_cmd=[]
        l_bam=[]
        for line in fin_list:
            sf_bam=line.rstrip()
            l_bam.append(sf_bam)
            cmd="samtools index {0}".format(sf_bam)
            l_cmd.append(cmd)

        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()

        l_cmd=[]
        for sf_bam in l_bam:
            sf_bam_status=sf_bam+".stats"
            cmd="samtools idxstats {0} > {1}".format(sf_bam, sf_bam_status)
            l_cmd.append(cmd)
        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()
####
####
def load_sample_fq(sf_fq_list):
    m_samples={}
    with open(sf_fq_list) as fin_fq:
        for line in fin_fq:
            fields=line.split()
            s_sample=fields[0]
            sf_lfq=fields[1]
            sf_rfq=fields[2]
            m_samples[s_sample]=(sf_lfq, sf_rfq)
    return m_samples
####

def parse_sample_id_all_lanes(sf_sample_fq):
    m_fq={}
    m_fq_lanes={}
    with open(sf_sample_fq) as fin_fq:
        for line in fin_fq:
            fields=line.split("/")#first split by "/", to get the file name and path
            s_prefix="/".join(fields[:-1])
            s_file=fields[-1]
            l_name_fields=s_file.split("_")
            s_id="_".join(l_name_fields[:-1])

            if s_id not in m_fq:
                m_fq[s_id]=[]
            m_fq[s_id].append((s_prefix, s_id))

            s_id_all_lanes = "_".join(l_name_fields[:-2])
            if s_id_all_lanes not in m_fq_lanes:
                m_fq_lanes[s_id_all_lanes]=[]
            m_fq_lanes[s_id_all_lanes].append(s_id)
    return m_fq_lanes
####
'''
Examples:
####Two lanes:
s12148_WB_1_USPD16097714-N707-AK428_HKG3YDSXX_L1_1.fastq.gz
s12148_WB_1_USPD16097714-N707-AK428_HKG3YDSXX_L1_2.fastq.gz
s12148_WB_1_USPD16097714-N707-AK428_HKG3YDSXX_L2_1.fastq.gz
s12148_WB_1_USPD16097714-N707-AK428_HKG3YDSXX_L2_2.fastq.gz

or (one lane)
/n/data1/hms/dbmi/park/DATA/Ting/pool1/C202SC19040221/raw_data/LAM11_Pre/LAM11_Pre_USPD16097422_HJN3MDSXX_L1_1.fq.gz
/n/data1/hms/dbmi/park/DATA/Ting/pool1/C202SC19040221/raw_data/LAM11_Pre/LAM11_Pre_USPD16097422_HJN3MDSXX_L1_2.fq.gz
'''
#parse sample_id from fastq file
def parse_sample_id(sf_sample_fq):
    m_fq={}
    b_fastq=True
    m_fq_lanes={}
    with open(sf_sample_fq) as fin_fq:
        for line in fin_fq:
            fields=line.split("/")#first split by "/", to get the file name and path
            s_prefix="/".join(fields[:-1])
            s_file=fields[-1]
            l_name_fields=s_file.split("_")
            s_id="_".join(l_name_fields[:-1])

            if s_id not in m_fq:
                m_fq[s_id]=[]
            if len(l_name_fields[-1])<=8:#x.fq.gz
                b_fastq=False
            m_fq[s_id].append((s_prefix, s_id))

            s_id_all_lanes = "_".join(l_name_fields[:-2])
            if s_id_all_lanes not in m_fq_lanes:
                m_fq_lanes[s_id_all_lanes]=1

    m_rtn={}
    for s_id in m_fq:
        if len(m_fq[s_id])>1:
            sf_lfq=m_fq[s_id][0][0]+"/{0}_1.fq.gz".format(m_fq[s_id][0][1])
            sf_rfq = m_fq[s_id][0][0] + "/{0}_2.fq.gz".format(m_fq[s_id][0][1])
            if b_fastq==True:
                sf_lfq = m_fq[s_id][0][0] + "/{0}_1.fastq.gz".format(m_fq[s_id][0][1])
                sf_rfq = m_fq[s_id][0][0] + "/{0}_2.fastq.gz".format(m_fq[s_id][0][1])
            m_rtn[s_id]=(sf_lfq, sf_rfq)
    return m_rtn
####


####
def merge_bam_from_different_lanes(sf_ref, sf_ori_fq, sf_bam_folder, i_mem, n_jobs, spartition, stime,
                                   sf_wfolder, sf_run_sh):
    m_id_all_lane=parse_sample_id_all_lanes(sf_ori_fq)
    l_sh=[]

    s_head = gnrt_script_head(spartition, n_jobs, stime, i_mem)
    for s_id_all_lane in m_id_all_lane:
        sf_tmp=sf_wfolder+s_id_all_lane+"_bam_list.txt"
        m_processed={}
        with open(sf_tmp,"w") as fout_tmp:
            n_lanes=len(m_id_all_lane[s_id_all_lane])/2
            if n_lanes==1:
                continue
            for s_tmp_id in m_id_all_lane[s_id_all_lane]:
                if s_tmp_id in m_processed:
                    continue
                m_processed[s_tmp_id]=1
                sf_tmp_bam = sf_bam_folder + "{0}.sorted.bam\n".format(s_tmp_id)
                fout_tmp.write(sf_tmp_bam)

        sf_merged_bam=sf_wfolder+s_id_all_lane+".merged.cram"
        s_cmd = "samtools merge -rf  --output-fmt CRAM --reference {0} {1} -b {2}\n".format(sf_ref, sf_merged_bam, sf_tmp)
        s_cmd += "samtools index {0}\n".format(sf_merged_bam)
        s_sh = sf_wfolder + "merge_{0}.sh".format(s_id_all_lane)
        with open(s_sh, "w") as fout_tmp_sh:
            fout_tmp_sh.write(s_head)
            fout_tmp_sh.write(s_cmd)
        l_sh.append(s_sh)

    with open(sf_run_sh, "w") as fout_sh:
        fout_sh.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            fout_sh.write("sbatch < " + s_sh + "\n")

####merge all the crams of the same id
####
def second_stage_merge(sf_ref, sf_ori_fq, sf_bam_folder, i_mem, n_jobs, spartition, stime,
                                   sf_wfolder, sf_run_sh):
    m_id_all_lane = parse_sample_id_all_lanes(sf_ori_fq)
    l_sh = []
    m_merged_per_sample={} #one sample one bam (merge all the lanes and also different read groups)
    s_head = gnrt_script_head(spartition, n_jobs, stime, i_mem)
    for s_id_all_lane in m_id_all_lane:
        n_lanes = len(m_id_all_lane[s_id_all_lane]) / 2
        if n_lanes == 1:
            continue
        # for i in range(1, n_lanes + 1):
        #     sf_tmp_bam = sf_bam_folder + "{0}_L{1}.sorted.bam\n".format(s_id_all_lane, i)
        s_unique_fields=s_id_all_lane.split("_")
        s_unique_id=s_unique_fields[0]
        if s_unique_id not in m_merged_per_sample:
            m_merged_per_sample[s_unique_id]=[]
        sf_merged_bam = sf_wfolder + s_id_all_lane + ".merged.cram"
        m_merged_per_sample[s_unique_id].append(sf_merged_bam)
####
    for s_uid in m_merged_per_sample:
        sf_tmp=sf_wfolder + s_uid + ".second_stage_merged.txt"
        sf_second_merged_bam=sf_wfolder + s_uid + ".second_stage_merged.cram"
        if len(m_merged_per_sample[s_uid])==1:
            if os.path.isfile(m_merged_per_sample[s_uid][0]):
                cmd="mv {0} {1}".format(m_merged_per_sample[s_uid][0], sf_second_merged_bam)
                Popen(cmd, shell=True, stdout=PIPE).communicate()
            continue

        with open(sf_tmp,"w") as fout_tmp:
            for sf_cram in m_merged_per_sample[s_uid]:
                fout_tmp.write(sf_cram+"\n")
        s_cmd = "samtools merge -rf  --output-fmt CRAM --reference {0} {1} -b {2}\n".format(sf_ref, sf_second_merged_bam,
                                                                                           sf_tmp)
        s_cmd += "samtools index {0}\n".format(sf_second_merged_bam)
        s_sh = sf_wfolder + "merge_{0}.sh".format(s_uid)
        with open(s_sh, "w") as fout_tmp_sh:
            fout_tmp_sh.write(s_head)
            fout_tmp_sh.write(s_cmd)
        l_sh.append(s_sh)

    with open(sf_run_sh, "w") as fout_sh:
        fout_sh.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            fout_sh.write("sbatch < " + s_sh + "\n")

####
def algn_sort_bam(sf_fq_list, sf_ref, i_mem, n_jobs, spartition, stime, sf_wfolder, sf_run_sh):
    m_samples = parse_sample_id(sf_fq_list)

    l_sh = []
    sf_cmd_folder = sf_wfolder + "cmd/"
    if os.path.exists(sf_cmd_folder) == False:
        cmd = "mkdir {0}".format(sf_cmd_folder)
        run_cvt(cmd)

    i_request_mem = i_mem + 3
    i_request_cores = n_jobs + 2
    s_head = gnrt_script_head(spartition, i_request_cores, stime, i_request_mem)

    for s_id in m_samples:
        rcd_fq = m_samples[s_id]
        sf_out_bam = sf_wfolder + s_id + ".sorted.bam"
        s_cmd = align_from_fq(rcd_fq, sf_ref, i_mem, sf_out_bam, sf_wfolder, n_jobs)

        s_sh = sf_cmd_folder + "algn_{0}.sh".format(s_id)
        with open(s_sh, "w") as fout_tmp_sh:
            fout_tmp_sh.write(s_head)
            fout_tmp_sh.write(s_cmd)
        l_sh.append(s_sh)
    with open(sf_run_sh, "w") as fout_sh:
        fout_sh.write("#!/bin/bash\n\n")
        for s_sh in l_sh:
            fout_sh.write("sbatch < " + s_sh + "\n")

    sf_run_sh_merge = sf_run_sh + ".merge_bam"
    merge_bam_from_different_lanes(sf_ref, sf_fq_list, sf_wfolder, i_mem, n_jobs, spartition, stime,
                                   sf_wfolder, sf_run_sh_merge)

    sf_run_sh_merge = sf_run_sh + ".second_merge_bam.sh"
    second_stage_merge(sf_ref, sf_fq_list, sf_wfolder, i_mem, n_jobs, spartition, stime,
                       sf_wfolder, sf_run_sh_merge)
####
def parse_option():
    parser = OptionParser()
    parser.add_option("-D", "--download",
                      action="store_true", dest="download", default=False,
                      help="Download the sra fiels from dbGaP (for runinfo_table)")
    parser.add_option("-L", "--download_list",
                      action="store_true", dest="download_list", default=False,
                      help="Download the sra fiels with given list")
    parser.add_option("-P", "--parse",
                      action="store_true", dest="parse", default=False,
                      help="Parse bam files from sra files")
    parser.add_option("-Q", "--fastq",
                      action="store_true", dest="fastq", default=False,
                      help="Download sra and parse out fastq")
    parser.add_option("-A", "--align",
                      action="store_true", dest="align_list", default=False,
                      help="Re-align a list of bams")
    parser.add_option("-I", "--index",
                      action="store_true", dest="index", default=False,
                      help="Index bam files")
    parser.add_option("-R", "--remap",
                      action="store_true", dest="remap", default=False,
                      help="Re-map from the bam")
    parser.add_option("-M", "--map",
                      action="store_true", dest="map", default=False,
                      help="map sort a given list of samples ")

    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="output file ", metavar="FILE")
    parser.add_option("-r", "--ref", dest="ref",
                      help="reference file ", metavar="FILE")
    parser.add_option("-p", "--wfolder", dest="wfolder",
                      help="Working folder", metavar="FILE")
    parser.add_option("-n", "--cores", dest="cores", type="int", default=8,
                      help="Number of cores")
    parser.add_option("-q", "--partition", dest="partition", type="string", default="park",
                      help="Which queue to run the job")
    parser.add_option("-m", "--mem", dest="memory", type="int", default=30,
                      help="Memory size in GB, should be a number")
    parser.add_option("-t", "--time", dest="time", type="string", default="3-0:00:00",
                      help="Time limit")
    (options, args) = parser.parse_args()
    return (options, args)

####MAIN####
if __name__ == '__main__':
    (options, args) = parse_option()

    if options.download:
        sf_list=options.input
        n_jobs=options.cores
        download_sra(sf_list, n_jobs)
    elif options.download_list:
        sf_list = options.input
        n_jobs = options.cores
        download_sra_with_list(sf_list, n_jobs)
    elif options.parse:
        sf_list = options.input
        s_wfolder=options.wfolder
        n_jobs = options.cores
        parse_bam(sf_list, n_jobs, s_wfolder)
    elif options.fastq:
        sf_list = options.input
        n_jobs = options.cores
        parse_fastq_from_sra(sf_list, n_jobs)

    elif options.align_list:#from a list of bams
        sf_bam_list = options.input
        sf_ref = options.ref
        sf_wfolder = options.wfolder
        if sf_wfolder[-1]!="/":
            sf_wfolder+="/"
        sf_run_sh=options.output

        spartition=options.partition
        i_mem=options.memory
        stime=options.time
        n_cores = options.cores

        realign_from_bam_list(sf_bam_list, sf_ref, spartition, i_mem, stime, n_cores, sf_wfolder, sf_run_sh)

    elif options.remap:#from a single bam
        sf_bam = options.input
        sf_ref=options.ref
        s_wfolder = options.wfolder
        n_jobs = options.cores

        if s_wfolder[-1]!="/":
            s_wfolder+="/"
        sf_tmp_folder=s_wfolder+"tmp"
        if os.path.exists(sf_tmp_folder)==False:
            cmd="mkdir {0}".format(sf_tmp_folder)
            run_cvt(cmd)

        s_mem=options.memory
        s_name=ntpath.basename(sf_bam)
        s_prefix=os.path.splitext(s_name)

        sf_out=s_wfolder+s_prefix[0]+".re_align.sorted.bam"
        realign_from_bam(sf_bam, sf_ref, s_mem, sf_out, sf_tmp_folder, n_jobs)

    elif options.map:
        sf_fq_list=options.input
        sf_ref=options.ref
        s_wfolder = options.wfolder
        n_jobs = options.cores
        if s_wfolder[-1]!="/":
            s_wfolder+="/"
        sf_tmp_folder=s_wfolder+"tmp"
        if os.path.exists(sf_tmp_folder)==False:
            cmd="mkdir {0}".format(sf_tmp_folder)

        s_mem = options.memory
        # s_name = ntpath.basename(sf_bam)
        # s_prefix = os.path.splitext(s_name)
        sf_out_sh = options.output
        spartition = options.partition
        stime = options.time
        algn_sort_bam(sf_fq_list, sf_ref, s_mem, n_jobs, spartition, stime, s_wfolder, sf_out_sh)

    elif options.index:#index and check status
        sf_list = options.input
        n_jobs = options.cores
        index_check_status(sf_list, n_jobs)
####