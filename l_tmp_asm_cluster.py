import os
import global_values
from optparse import OptionParser
from cluster_helper.cluster import cluster_view
import subprocess, select
from l_MEI_caller import *

# this version align the flank regions to contig sequence
# Note: in some cases the called out breakpoint is not accuracy, which may cause problem
def call_MEIs_for_sites(sf_bam_list, sf_sites, sf_ref, l_extd_len, wfolder, n_cores, sf_out_fa, sf_out_sites, l_cluster_info=None):
    '''
    :param sf_bam_list: a list of long reads bam/cram files
    :param sf_sites: candidate site list
    :param sf_ref: reference genome
    :param l_extd_len: list of extended flanking region length for asm
    :param sf_out_fa: constructed fa file
    :param sf_out_sites: sites with insertion seq constructed
    :return: None
    '''
    l_tei_caller=L_MEI_Caller(wfolder, n_cores,sf_ref)

    sf_tmp = wfolder + "l_asm_tmp/"
    if os.path.exists(sf_tmp) == False:
        cmd = "mkdir {0}".format(sf_tmp)
        run_cmd_small_output(cmd)
    l_bams = _load_in_bam(sf_bam_list)#
    if os.path.isfile(sf_out_fa) == True:  # as everytime each called seq is save in "append" way
        os.remove(sf_out_fa)

    m_constructed = {}
    m_new_pos={}
    n_round = 1
    for sf_bam in l_bams:  # run for different bams
        for i_ext_len in l_extd_len:  # run several rounds with diff flanking length
            l_sites = []
            sf_tmp_to_asm=sf_sites+".{0}_to_asm_list.txt".format(i_ext_len)
            with open(sf_sites) as fin_sites, open(sf_tmp_to_asm,"w") as fout_for_asm:
                for line in fin_sites:
                    fields = line.split()

                    ins_chrm = fields[0]
                    ins_pos = int(fields[1])
                    # skip already constructed ones
                    if (ins_chrm in m_constructed) and (str(ins_pos) in m_constructed[ins_chrm]):
                        continue
                    rcd = (sf_bam, ins_chrm, ins_pos, sf_tmp, n_cores)
                    l_sites.append(rcd)
                    fout_for_asm.write(line)
            global_values.set_lrd_extnd_len(i_ext_len)

            if len(l_sites)<=0:
                continue
####comment out for temporary test!!!!!!!!
            l_tei_caller.collect_seq_for_sites_in_parallel(l_sites)

            s_scheduler=l_cluster_info[0]
            s_queue=l_cluster_info[1]
            n_parallel_jobs=l_cluster_info[2]
            n_core_per_job=l_cluster_info[3]
            n_mem_GB=l_cluster_info[4]
            i_wait_mins=l_cluster_info[5]
            s_max_time=l_cluster_info[6]
            s_total_mem=l_cluster_info[7]
            asm_seq_for_sites_in_parallel_on_cluster(l_sites, s_scheduler, s_queue, n_parallel_jobs,
                                                          n_core_per_job, n_mem_GB, i_wait_mins, s_max_time, s_total_mem)

            m_tmp, m_refined_pos = l_tei_caller.call_MEIs_from_asm_in_parallel(sf_ref, sf_tmp_to_asm, l_sites, sf_tmp,
                                                                       sf_out_fa)
            for s_site in m_tmp:
                s_tmp_fields = s_site.split(global_values.SEPERATOR)
                if s_tmp_fields[0] not in m_constructed:
                    m_constructed[s_tmp_fields[0]] = {}
                m_constructed[s_tmp_fields[0]][s_tmp_fields[1]] = 1
                if s_tmp_fields[0] not in m_new_pos:
                    m_new_pos[s_tmp_fields[0]]={}
                m_new_pos[s_tmp_fields[0]][s_tmp_fields[1]]=m_refined_pos[s_site]
            print "Round {0}: {1} sites (out of {2}) are constructed\n".format(n_round, len(m_tmp), len(l_sites))
            n_round += 1

    with open(sf_out_sites, "w") as fout_sites:
        for chrm in m_constructed:
            for pos in m_constructed[chrm]:
                new_pos=m_new_pos[chrm][pos]
                s_site_info = "{0}\t{1}\n".format(chrm, new_pos)
                fout_sites.write(s_site_info)

####
def _load_in_bam(sf_bam_list):
    l_bams = []
    with open(sf_bam_list) as fin_bams:
        for line in fin_bams:
            sf_bam = line.rstrip()
            l_bams.append(sf_bam)
    return l_bams

# def _test_func(cmd):
#     import os
#     print "Running command with output: {0}\n".format(cmd)
#     sf_err = sf_std_out + ".err"
#     errcode = None
#     with open(sf_std_out, "w") as f:
#         p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=None)
#         # errcode = p.wait()
#         p.communicate()
#
#     if os.path.isfile(sf_err):
#         os.remove(sf_err)

####submit jobs to the cluster
def asm_seq_for_sites_in_parallel_on_cluster(l_sites, s_scheduler, s_queue, n_jobs,
                                             n_core_per_job, n_mem_GB,
                                             i_wait_mins, s_max_time, s_total_mem):
    s_extra = s_max_time + ";" + s_total_mem
    with cluster_view(scheduler=s_scheduler, queue=s_queue, num_jobs=n_jobs,
                      cores_per_job=n_core_per_job,
                      start_wait=i_wait_mins,
                      extra_params={'mem': n_mem_GB, 'resources': s_extra}) as view:
        results = list(view.map(_asm_collected_reads_one_site, l_sites))
        print results
        ####


# assemble the collected clipped and contained reads
def _asm_collected_reads_one_site(record):
    import os
    sf_bam = record[0]
    chrm = record[1]
    pos = record[2]
    wfolder = record[3]
    ncores=record[4]

    SEPERATOR = '~'
    LRD_EXTND_LEN = 2000
    FAIL_CLP = 'COLLECT_CLIP_FAIL'
    WTDBG2 = "wtdbg2"
    WTPOA = "wtpoa-cns"
    FAIL_ASM = 'ASM_FAIL'
    FAIL_CNS = "POLISH_FAIL"
    SUCCD = "succeed"

    if wfolder[-1] != "/":
        wfolder += "/"
    s_id = "{0}{1}{2}".format(chrm, SEPERATOR, pos)
    i_extend =LRD_EXTND_LEN
    s_prefix = wfolder + s_id + "_{0}".format(i_extend)
    sf_fa = s_prefix + ".fa"
    print sf_fa, "doesnt' exist!!!!"
    if os.path.isfile(sf_fa) == False:
        return FAIL_CLP
    cmd = "{0}  -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 --ctg-min-nodes 1 " \
          "-i {1}.fa  -fo {2}_wtdbg2 -q -t {3}".format(WTDBG2, s_prefix, s_prefix, ncores)
    print cmd
    # with open(s_prefix + ".fa.cmd","w") as fout:
    #     fout.write(cmd)
    sf_std_out=sf_fa+".std_out"

    #run_cmd_to_file(cmd, sf_std_out)
    print "Running command with output: {0}\n".format(cmd)
    sf_err = sf_std_out + ".err"
    errcode = None
    with open(sf_std_out, "w") as f:
        p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=None)
        # errcode = p.wait()
        p.communicate()

    if os.path.isfile(sf_err):
        os.remove(sf_err)

    sf_lay = "{0}_wtdbg2.ctg.lay".format(s_prefix)
    cmd=""
    if os.path.isfile(sf_lay) == False:
        if os.path.isfile(sf_lay+".gz") == True:
            cmd = "{0} -i {1}_wtdbg2.ctg.lay.gz -fo {2}_wtdbg2.ctg.lay.fa -t {3}".format(WTPOA,
                                                                                         s_prefix, s_prefix, ncores)
        else:
            return FAIL_ASM
    else:
        cmd = "{0} -i {1}_wtdbg2.ctg.lay -fo {2}_wtdbg2.ctg.lay.fa -t {3}".format(WTPOA, s_prefix,
                                                                                  s_prefix, ncores)
    print cmd
    if cmd!="":
        print "Running command with output: {0}\n".format(cmd)
        sf_err = sf_std_out + ".err"
        errcode = None
        with open(sf_std_out, "w") as f:
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=None)
            # errcode = p.wait()
            p.communicate()
        if os.path.isfile(sf_err):
            os.remove(sf_err)

    sf_lay_fa = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
    if os.path.isfile(sf_lay_fa) == False:
        return FAIL_CNS
    return SUCCD

####
def run_cmd_small_output(cmd):
    print "Running command: {0}\n".format(cmd)
    subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()

def merge_sites(sf_sites, sf_sites2, i_merge_slack):
    sf_merged=sf_sites
    if os.path.isfile(sf_sites2)==False:
        return sf_merged
    m_ori_sites={}
    with open(sf_sites) as fin_sites:
        for line in fin_sites:
            fields=line.split()
            ins_chrm=fields[0]
            ins_pos=int(fields[1])
            if ins_chrm not in m_ori_sites:
                m_ori_sites[ins_chrm]={}
            if ins_pos not in m_ori_sites[ins_chrm]:
                m_ori_sites[ins_chrm][ins_pos]=line.rstrip()
    sf_merged=sf_sites+".merged_with_external_source.txt"

    with open(sf_sites2) as fin_rslt2:
        for line in fin_rslt2:
            fields=line.split()
            ins_chrm=fields[0]
            ins_pos=int(fields[1])

            if ins_chrm not in m_ori_sites:
                m_ori_sites[ins_chrm]={}
                m_ori_sites[ins_chrm][ins_pos]=line.rstrip()
                continue

            istart=-1*i_merge_slack
            iend=i_merge_slack
            b_exist=False
            for i in range(istart, iend):
                itmp=i+ins_pos
                if itmp in m_ori_sites[ins_chrm]:
                    b_exist=True
                    break
            if b_exist==False:
                m_ori_sites[ins_chrm][ins_pos] = line.rstrip()

    with open(sf_merged, "w") as fout_merged:
        for ins_chrm in m_ori_sites:
            for ins_pos in m_ori_sites[ins_chrm]:
                fout_merged.write(m_ori_sites[ins_chrm][ins_pos]+"\n")
    return sf_merged

####
##parse the options options
def parse_option():
    parser = OptionParser()

    parser.add_option("-U", "--user",
                      action="store_true", dest="user_specific", default=False,
                      help="User specific cutoff")
    parser.add_option("-C", "--lrd_clip",
                      action="store_true", dest="lrd_clip", default=False,
                      help="Collect the potential clip positions")
    parser.add_option("-A", "--assembly",
                      action="store_true", dest="assembly", default=False,
                      help="Do local assembly for collected long reads")
    parser.add_option("-P", "--Pinpoint",
                      action="store_true", dest="pinpoint", default=False,
                      help="Pinpoint the exact postion and call the insertion seq from assembled seqs")
    parser.add_option("-Y", "--classify",
                      action="store_true", dest="classify", default=False,
                      help="Classify the insertion by alignment to templates")
    parser.add_option("-N", "--polymorphic",
                      action="store_true", dest="polymorphic", default=False,
                      help="Call candidate non reference polymorphic rep copies")
    parser.add_option("--line",
                      action="store_true", dest="call_LINE", default=False,
                      help="Call LINE insertions from rmsk output")
    parser.add_option("--sva",
                      action="store_true", dest="call_SVA", default=False,
                      help="Call SVA insertions from rmsk output")
    parser.add_option("-D", "--delete",
                      action="store_true", dest="delete", default=False,
                      help="Remove the intermediate files")
    parser.add_option("-K", "--keep",
                      action="store_true", dest="keep", default=False,
                      help="Keep the intermediate files")
    parser.add_option("--pacbio",
                      action="store_true", dest="pacbio", default=True,
                      help="Pacbio reads")
    parser.add_option("--hg19",
                      action="store_true", dest="hg19", default=True,
                      help="Working on reference genome hg19")
    parser.add_option("--cns", dest="consensus",
                      help="repeat consensus file", metavar="FILE")

    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-r", "--ref", dest="reference",
                      help="Reference genome", metavar="FILE") #/or reference consensus
    parser.add_option("--rep", dest="rep_lib",
                      help="repeat folder", metavar="FILE")
    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("--i2", dest="input2", default="null",
                      help="input file ", metavar="FILE")
    parser.add_option("-m", "--rmsk", dest="rmsk",
                      help="input file ", metavar="FILE")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")
    parser.add_option("--min", dest="min_copy_len", type="int", default=1000,
                      help="Minimum copy length for collecting polymorphic reads")
    parser.add_option("-s", "--slack", dest="slack", type="int", default=250,
                      help="slack value")
    parser.add_option("-w", "--win", dest="win", type="int", default=75,
                      help="peak window size")
    parser.add_option("-d", "--std", dest="std", type="int", default=50,
                      help="Maximum standard derivation of breakpoints")
    parser.add_option("-c", dest="clip", type="int", default=7,
                      help="cutoff of minimum # of clipped/contained reads")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-y", "--type", dest="type", type="int", default=7,
                      help="Type of repeats working on")
    (options, args) = parser.parse_args()
    return (options, args)

####

####
if __name__ == '__main__':
    (options, args) = parse_option()
    sf_bam_list = options.bam
    sf_sites = options.input
    sf_sites2 = options.input2  # another list from a different source
    n_cores = options.cores
    sf_out_fa = options.output
    sf_out_sites = sf_out_fa + ".sites"
    sf_ref = options.reference
    swfolder = options.wfolder
    if swfolder[-1] != "/":
        swfolder += "/"

    i_slack = 150  # will be merged if distance is smaller than this value
    sf_merged = merge_sites(sf_sites, sf_sites2, i_slack)

    # negative means collect whole reads, and align abs(size) flanking region
    # l_extd_len = [3500, 2500, 1500, 500, 350]
    l_extd_len = [3500, 1500, 500, -1000]

    s_scheduler = "lsf"
    s_queue = "research-rh74"
    n_parallel_jobs = 100
    n_core_per_job = n_cores
    n_mem_GB = 10
    i_wait_mins = 200  # wait for minutes
    s_max_running_time = "W 12:00"
    s_peak_mem = "M {0}G".format(n_mem_GB * n_cores)
    l_cluster_info = [s_scheduler, s_queue, n_parallel_jobs, n_core_per_job, n_mem_GB, i_wait_mins,
                      s_max_running_time, s_peak_mem]
    # align flank regions to the contig
    call_MEIs_for_sites(sf_bam_list, sf_merged, sf_ref, l_extd_len, swfolder, n_cores, sf_out_fa, sf_out_sites, l_cluster_info)