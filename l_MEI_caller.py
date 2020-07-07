##01/07/2010
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import sys
import pysam
from multiprocessing import Pool
import multiprocessing
from optparse import OptionParser
from l_asm import *
import cmd_runner
from rmsk_parser import *
import global_values
from l_clip_collect import *
from x_reference import *
from l_contig import *
from l_algnmt_breakpoints import *
from x_intermediate_sites import *
from l_rep_classification import *
from x_alignments import *
from l_complex_sv import *
#from cluster_helper.cluster import cluster_view

####
import copy_reg
import types

def _reduce_method(m):#
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)
copy_reg.pickle(types.MethodType, _reduce_method)
####

####
# def unwrap_self_collect_reads_one_site(arg, **kwarg):
#     return L_MEI_Caller.collect_clip_contain_reads_one_site(*arg, **kwarg)

def unwrap_self_remove_intermediate_files(arg, **kwarg):
    return L_MEI_Caller.clean_intermediate_output_one_sample(*arg, **kwarg)

def unwrap_self_remove_intermediate_flank_files(arg, **kwarg):
    return L_MEI_Caller.clean_intermediate_output_flanks_one_sample(*arg, **kwarg)

def unwrap_self_asm_collected_reads_one_site(arg, **kwarg):
    return L_MEI_Caller._asm_collected_reads_one_site(*arg, **kwarg)

class L_MEI_Caller():
    def __init__(self, s_wfolder, n_jobs, sf_ref):
        if s_wfolder[-1] != "/":
            s_wfolder += "/"
        self.wfolder = s_wfolder
        self.n_jobs = n_jobs
        self.sf_ref = sf_ref
        self.BRKPNT_FOLDER="l_brkpnts_tmp"
        self.ASM_FOLDER="l_asm_tmp"
        self.COMPLEX_SV_FOLDER="l_complex_sv"

    ####
    def collect_MEIs_breakpoints(self, sf_bam_list, sf_ref, i_peak_win, i_cutoff, i_max_cutoff, std_dev_cutoff,
                                 sf_merged_out, sf_merged_out_left, sf_merged_out_right):
        '''
        Function: call candidate insertion candidate sites
        :param sf_bam_list: bam file list
        :param sf_ref:  reference genome
        :param i_peak_win: if two breakpoints within this distance, then they will be merged
        :param i_cutoff: minimum of the supported clipped and contained reads
        :param sf_merged_out: output merged candidate sites
        Revised: 11/14/2019, save breakpoints only have left (right) clipped reads
        '''
        sf_tmp = self.wfolder + self.BRKPNT_FOLDER+"/"
        cmd_runner = CMD_RUNNER()
        if os.path.exists(sf_tmp) == False:
            cmd = "mkdir {0}".format(sf_tmp)
            cmd_runner.run_cmd_small_output(cmd)
        l_bams = self._load_in_bam(sf_bam_list)

        m_merged = {}  # merge the breakpoints gotten from each bam
        m_merged_left={}#
        m_merged_right={}#
        for sf_bam in l_bams:
            sf_base = os.path.basename(sf_bam)
            sf_tmp_out = sf_tmp + sf_base + ".brkpnt_pos"
            sf_tmp_out_l = sf_tmp + sf_base + ".left_brkpnt_pos"#only have right clipped reads
            sf_tmp_out_r = sf_tmp + sf_base + ".right_brkpnt_pos"#only have left clipped reads
            # each bam create a separate folder
            sf_bam_folder = sf_tmp + sf_base
            if os.path.exists(sf_bam_folder) == False:
                cmd = "mkdir {0}".format(sf_bam_folder)
                cmd_runner.run_cmd_small_output(cmd)
            lrd_brkpnts = LRDBreakpoints(self.n_jobs, sf_bam, sf_ref, sf_bam_folder)
            lrd_brkpnts.collect_breakpoints(i_peak_win, i_cutoff, i_max_cutoff, std_dev_cutoff, sf_tmp_out,
                                            sf_tmp_out_l, sf_tmp_out_r)
            self._load_sites_from_file(sf_tmp_out, m_merged)
            self._load_sites_from_file(sf_tmp_out_l, m_merged_left)
            self._load_sites_from_file(sf_tmp_out_r, m_merged_right)

        x_intermediate_sites = XIntemediateSites()
        x_intermediate_sites.output_candidate_sites(m_merged, sf_merged_out)
        x_intermediate_sites.output_candidate_sites(m_merged_left, sf_merged_out_left)
        x_intermediate_sites.output_candidate_sites(m_merged_right, sf_merged_out_right)
####
####
    def _load_sites_from_file(self, sf_tmp_out, m_merged):
        x_intermediate_sites = XIntemediateSites()
        m_bam_sites = x_intermediate_sites.load_in_candidate_list_str_version(sf_tmp_out)
        for tmp_chrm in m_bam_sites:
            if tmp_chrm not in m_merged:
                m_merged[tmp_chrm] = {}
            for tmp_pos in m_bam_sites[tmp_chrm]:
                if tmp_pos not in m_merged[tmp_chrm]:
                    m_merged[tmp_chrm][tmp_pos] = []
                for elmnt in m_bam_sites[tmp_chrm][tmp_pos]:
                    m_merged[tmp_chrm][tmp_pos].append(elmnt)
####
####
    def merge_sites(self, sf_sites, sf_sites2, i_merge_slack):
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
    # this version align the flank regions to contig sequence
    # Note: in some cases the called out breakpoint is not accuracy, which may cause problem
    def call_MEIs_for_sites(self, sf_bam_list, sf_sites, sf_ref, l_extd_len, sf_out_fa, sf_out_sites, l_cluster_info=None):
        '''
        :param sf_bam_list: a list of long reads bam/cram files
        :param sf_sites: candidate site list
        :param sf_ref: reference genome
        :param l_extd_len: list of extended flanking region length for asm
        :param sf_out_fa: constructed fa file
        :param sf_out_sites: sites with insertion seq constructed
        :return: None
        '''
        sf_tmp = self.wfolder + self.ASM_FOLDER+"/"
        if os.path.exists(sf_tmp) == False:
            cmd = "mkdir {0}".format(sf_tmp)
            cmd_runner = CMD_RUNNER()
            cmd_runner.run_cmd_small_output(cmd)
        l_bams = self._load_in_bam(sf_bam_list)#
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
                        rcd = (sf_bam, ins_chrm, ins_pos, sf_tmp)
                        l_sites.append(rcd)
                        fout_for_asm.write(line)
                global_values.set_lrd_extnd_len(i_ext_len)

                if len(l_sites)<=0:
                    continue
                ####comment out for temporary test!!!!!!!!
                self.collect_seq_for_sites_in_parallel(l_sites)
                #if l_cluster_info is None:
                self.asm_seq_for_sites_in_serial(l_sites)
                # else:#this version doesn't work for some reason!!!!
                #     s_scheduler=l_cluster_info[0]
                #     s_queue=l_cluster_info[1]
                #     n_parallel_jobs=l_cluster_info[2]
                #     n_core_per_job=l_cluster_info[3]
                #     n_mem_GB=l_cluster_info[4]
                #     i_wait_mins=l_cluster_info[5]
                #     s_max_time=l_cluster_info[6]
                #     s_total_mem=l_cluster_info[7]
                #     self.asm_seq_for_sites_in_parallel_on_cluster(l_sites, s_scheduler, s_queue, n_parallel_jobs,
                #                                                   n_core_per_job, n_mem_GB, i_wait_mins, s_max_time, s_total_mem)
                m_tmp, m_refined_pos = self.call_MEIs_from_asm_in_parallel(sf_ref, sf_tmp_to_asm, l_sites, sf_tmp, sf_out_fa)

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
    ####
    # this version align the flank regions to contig sequence
    # Note: in some cases the called out breakpoint is not accuracy, which may cause problem
    def call_MEIs_for_sites_with_cmd_version(self, sf_bam_list, sf_sites, sf_ref, l_extd_len, sf_out_fa, sf_out_sites,
                                             b_complex, b_asm_cmd=False, sf_script_prefix="", b_call_seq=True,
                                             b_skip_exist=False):
        '''
        :param sf_bam_list: a list of long reads bam/cram files
        :param sf_sites: candidate site list
        :param sf_ref: reference genome
        :param l_extd_len: list of extended flanking region length for asm
        :param sf_out_fa: constructed fa file
        :param sf_out_sites: sites with insertion seq constructed
        :return: None
        '''
        sf_tmp_prefix = self.wfolder + self.ASM_FOLDER
        if os.path.exists(sf_tmp_prefix) == False:
            cmd = "mkdir {0}".format(sf_tmp_prefix)
            cmd_runner = CMD_RUNNER()
            cmd_runner.run_cmd_small_output(cmd)
        l_bams = self._load_in_bam(sf_bam_list)  #
        if os.path.isfile(sf_out_fa) == True:  # as everytime each called seq is save in "append" way
            os.remove(sf_out_fa)

        m_constructed = {}
        m_new_pos = {}
        n_round = 1
        for sf_bam in l_bams:  # run for different bams
            for i_ext_len in l_extd_len:  # run several rounds with diff flanking length
                l_sites = []
                sf_tmp = sf_tmp_prefix + "/{0}/".format(i_ext_len)
                if i_ext_len<0:
                    sf_tmp=sf_tmp_prefix+"/neg_{0}/".format(abs(i_ext_len))
                if os.path.exists(sf_tmp) == False:
                    cmd = "mkdir {0}".format(sf_tmp)
                    cmd_runner = CMD_RUNNER()
                    cmd_runner.run_cmd_small_output(cmd)

                sf_tmp_to_asm = sf_sites + ".{0}_to_asm_list.txt".format(i_ext_len)
                with open(sf_sites) as fin_sites, open(sf_tmp_to_asm, "w") as fout_for_asm:
                    for line in fin_sites:
                        fields = line.split()
                        ins_chrm = fields[0]
                        ins_pos = int(fields[1])
                        # skip already constructed ones
                        if (ins_chrm in m_constructed) and (str(ins_pos) in m_constructed[ins_chrm]):
                            continue

                        #also check whether already assembled
                        if b_skip_exist==True:
                            l_lasm=L_Local_ASM()
                            if l_lasm.is_wtdbg2_asm_already_exist(ins_chrm, ins_pos, sf_tmp)==True:
                                continue

                        rcd = (sf_bam, ins_chrm, ins_pos, sf_tmp)
                        l_sites.append(rcd)
                        fout_for_asm.write(line)
                global_values.set_lrd_extnd_len(i_ext_len)
                if len(l_sites) <= 0:
                    continue
####
                #collect the clipped and contained reads for sites
                if b_call_seq == False:
                    self.collect_seq_for_sites_in_parallel(l_sites)
                # if l_cluster_info is None:
                if b_asm_cmd==True:
                    sf_script=sf_script_prefix+"_{0}.sh".format(i_ext_len)
                    self.gnrt_asm_cmd_only(l_sites, sf_script)
                    continue

                if b_call_seq==False:#if this is only for call insertion from asm, then no need for asm
                    self.asm_seq_for_sites_in_serial(l_sites)

                # else:#this version doesn't work for some reason!!!!
                #     s_scheduler=l_cluster_info[0]
                #     s_queue=l_cluster_info[1]
                #     n_parallel_jobs=l_cluster_info[2]
                #     n_core_per_job=l_cluster_info[3]
                #     n_mem_GB=l_cluster_info[4]
                #     i_wait_mins=l_cluster_info[5]
                #     s_max_time=l_cluster_info[6]
                #     s_total_mem=l_cluster_info[7]
                #     self.asm_seq_for_sites_in_parallel_on_cluster(l_sites, s_scheduler, s_queue, n_parallel_jobs,
                #                                                   n_core_per_job, n_mem_GB, i_wait_mins, s_max_time, s_total_mem)

                # for comlex events, only collect and assembly, but do not call insertion seqs
                if b_complex==True:
                    continue

                m_tmp, m_refined_pos = self.call_MEIs_from_asm_in_parallel(sf_ref, sf_tmp_to_asm,
                                                                           l_sites, sf_tmp, sf_out_fa)
####
                for s_site in m_tmp:
                    s_tmp_fields = s_site.split(global_values.SEPERATOR)
                    if s_tmp_fields[0] not in m_constructed:
                        m_constructed[s_tmp_fields[0]] = {}
                    m_constructed[s_tmp_fields[0]][s_tmp_fields[1]] = 1
                    if s_tmp_fields[0] not in m_new_pos:
                        m_new_pos[s_tmp_fields[0]] = {}
                    m_new_pos[s_tmp_fields[0]][s_tmp_fields[1]] = m_refined_pos[s_site]
                print "Round {0}: {1} sites (out of {2}) are constructed\n".format(n_round, len(m_tmp),
                                                                                   len(l_sites))
                n_round += 1

        # for comlex events, only collect and assembly, but do not call insertion seqs
        if b_complex==True:
            return

        with open(sf_out_sites, "w") as fout_sites:
            for chrm in m_constructed:
                for pos in m_constructed[chrm]:
                    new_pos = m_new_pos[chrm][pos]
                    s_site_info = "{0}\t{1}\n".format(chrm, new_pos)
                    fout_sites.write(s_site_info)

####
    # this version align the contigs to target reference genome
    def call_MEIs_for_sites_with_contig_2_ref(self, sf_bam_list, sf_sites, sf_ref, l_extd_len, sf_out_fa, sf_out_sites):
        '''
        :param sf_bam_list: a list of long reads bam/cram files
        :param sf_sites: candidate site list
        :param sf_ref: reference genome
        :param l_extd_len: list of extended flanking region length for asm
        :param sf_out_fa: constructed fa file
        :param sf_out_sites: sites with insertion seq constructed
        :return: None
        '''
        sf_tmp = self.wfolder + self.ASM_FOLDER+"/"
        if os.path.exists(sf_tmp) == False:
            cmd = "mkdir {0}".format(sf_tmp)
            cmd_runner = CMD_RUNNER()
            cmd_runner.run_cmd_small_output(cmd)
        l_bams = self._load_in_bam(sf_bam_list)
        if os.path.isfile(sf_out_fa) == True:  # as everytime each called seq is save in "append" way
            os.remove(sf_out_fa)

        m_constructed = {}
        n_round = 1
        for sf_bam in l_bams:  # run for different bams
            for i_ext_len in l_extd_len:  # run several rounds with diff flanking length
                l_sites = []
                with open(sf_sites) as fin_sites:
                    for line in fin_sites:
                        fields = line.split()

                        ins_chrm = fields[0]
                        ins_pos = int(fields[1])
                        s_sites = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                        if s_sites in m_constructed:  # skip already constructed ones
                            continue
                        rcd = (sf_bam, ins_chrm, ins_pos, sf_tmp)
                        l_sites.append(rcd)
                global_values.set_lrd_extnd_len(i_ext_len)

                self.collect_seq_for_sites_in_parallel(l_sites)
                self.asm_seq_for_sites_in_serial(l_sites)
                m_tmp = self.call_MEIs_from_asm_contig_2_ref_in_parallel(sf_ref, sf_sites, l_sites, sf_tmp, sf_out_fa)
                for s_site in m_tmp:
                    s_tmp_fields = s_site.split(global_values.SEPERATOR)
                    if s_tmp_fields[0] not in m_constructed:
                        m_constructed[s_tmp_fields[0]] = {}
                    m_constructed[s_tmp_fields[0]][s_tmp_fields[1]] = 1
                print "Round {0}: {1} sites (out of {2}) are constructed\n".format(n_round, len(m_tmp),
                                                                                   len(l_sites))
                n_round += 1

        with open(sf_out_sites, "w") as fout_sites:
            for chrm in m_constructed:
                for pos in m_constructed[chrm]:
                    s_site_info = "{0}\t{1}\n".format(chrm, pos)
                    fout_sites.write(s_site_info)


    ####classify the assembled contigs
    def classify_asm_insertions(self, sf_fa, l_rep_types, l_sf_rep, sf_out):
        '''
        :param sf_fa: the input assembled insertion seq, with head in format:
        :param sf_rep: this is the repeat consensus seq designed for masking the repeat region
        :param sf_out: this is final output file
        :return:
        '''
        return

    ####
    def _load_in_bam(self, sf_bam_list):
        l_bams = []
        with open(sf_bam_list) as fin_bams:
            for line in fin_bams:
                sf_bam = line.rstrip()
                l_bams.append(sf_bam)
        return l_bams

    ####Note:1. the reads collection step is parallized
    def collect_seq_for_sites_in_parallel(self, l_sites):
        lr_clip_info = LRDClipReadInfo(self.n_jobs, sf_ref=self.sf_ref)
        lr_clip_info.collect_seqs_for_sites_in_parallel(l_sites)

    def collect_seq_cover_region_in_parallel(self, l_sites):
        lr_clip_info = LRDClipReadInfo(self.n_jobs, sf_ref=self.sf_ref)
        lr_clip_info.collect_seqs_cover_regions_in_parallel(l_sites)
####
    def gnrt_asm_cmd_only(self, l_sites, sf_script):
        lasm = L_Local_ASM()
        l_asm_cmd=[]
        l_cns_cmd=[]
        wfolder=""
        for record in l_sites:
            sf_bam = record[0]
            ins_chrm = record[1]
            ins_pos = record[2]
            wfolder = record[3]
            cmd1, cmd2=lasm.construct_short_cns_wtdbg2_cmd_list(ins_chrm, ins_pos, self.n_jobs, wfolder)
            l_asm_cmd.append(cmd1)
            l_cns_cmd.append(cmd2)

        with open(sf_script, "w") as fout_script:
            for asm_cmd in l_asm_cmd:
                cmd1="bsub -n 8 -W 1:00 -M 8G -R \"rusage[mem=1G]\" \"{0}\"\n".format(asm_cmd)
                fout_script.write(cmd1)
            for cns_cmd in l_cns_cmd:
                cmd2="bsub -n 8 -W 1:00 -M 8G -R \"rusage[mem=1G]\" \"{0}\"\n".format(cns_cmd)
                fout_script.write(cmd2)
            # if wfolder!="":
            #     cmd_rm = "find {0} -type f ! -name \'*.fa\' -delete".format(wfolder)
            #     cmd3 = "bsub -n 1 -W 5:00 -M 5G -R \"rusage[mem=5G]\" \"{0}\"\n".format(cmd_rm)
            #     fout_script.write(cmd3)
####

    def gnrt_asm_cmd_covered_region_only(self, l_sites, sf_script):
        lasm = L_Local_ASM()
        l_asm_cmd = []
        l_cns_cmd = []
        wfolder = ""
        for record in l_sites:
            sf_bam = record[0]
            ins_chrm = record[1]
            ins_pos = record[2]#region start
            ins_end=record[3]#region end
            wfolder = record[4]
            cmd1, cmd2 = lasm.construct_short_cns_wtdbg2_cmd_list(ins_chrm, ins_pos, self.n_jobs, wfolder)
            l_asm_cmd.append(cmd1)
            l_cns_cmd.append(cmd2)

        with open(sf_script, "w") as fout_script:
            for asm_cmd in l_asm_cmd:
                cmd1 = "bsub -n 8 -W 1:00 -M 8G -R \"rusage[mem=1G]\" \"{0}\"\n".format(asm_cmd)
                fout_script.write(cmd1)
            for cns_cmd in l_cns_cmd:
                cmd2 = "bsub -n 8 -W 1:00 -M 8G -R \"rusage[mem=1G]\" \"{0}\"\n".format(cns_cmd)
                fout_script.write(cmd2)

    ######## 2. the assembly step will fail if run in parallel !!!!!!!!!!!!!!!!!!!!!!!
    def asm_seq_for_sites_in_serial(self, l_sites):
        for rcd in l_sites:#
            self._asm_collected_reads_one_site(rcd)  # asm the sites

    # ####submit jobs to the cluster
    # def asm_seq_for_sites_in_parallel_on_cluster(self, l_sites, s_scheduler, s_queue, n_jobs, n_core_per_job, n_mem_GB,
    #                                              i_wait_mins, s_max_time, s_total_mem):
    #     s_extra=s_max_time+";"+s_total_mem
    #     with cluster_view(scheduler=s_scheduler, queue=s_queue, num_jobs=n_jobs, cores_per_job=n_core_per_job,
    #                       start_wait=i_wait_mins, extra_params={'mem': n_mem_GB, 'resources':s_extra}) as view:
    #         #view.map(unwrap_self_asm_collected_reads_one_site, zip([self] * len(l_sites), l_sites))
    #         view.map(self._asm_collected_reads_one_site, l_sites)
####
    # assemble the collected clipped and contained reads
    def _asm_collected_reads_one_site(self, record):
        sf_bam = record[0]
        ins_chrm = record[1]
        ins_pos = record[2]
        wfolder = record[3]
        lasm = L_Local_ASM()
        lasm.construct_short_cns_wtdbg2(ins_chrm, ins_pos, self.n_jobs, wfolder)  # by default use 4 cores

####
    ######## 3. call out from asm
    # each rcd of l_sites in format: (sf_bam, ins_chrm, ins_pos, sf_tmp)
    def call_MEIs_from_asm_in_parallel(self, sf_ref, sf_sites, l_sites, s_working_folder, sf_mei):
        # 3.1. collect flanking regions
        # 3.2. align the flanking regions to the assembled contigs
        i_extend = abs(global_values.LRD_EXTND_LEN) #there is a case of negative (to collect whole reads)
        xref = XReference()
        xref.gnrt_flank_region_for_sites(sf_sites, i_extend, self.n_jobs, sf_ref, s_working_folder)
        xctg = LRD_Contig(s_working_folder, self.n_jobs)

        sflank_folder = s_working_folder + global_values.FLANK_FOLDER
        l_rcd = []
        l_prefix_clean = []
        l_site_flanks = []

        for site_rcd in l_sites:
            # each in format: (sf_bam, ins_chrm, ins_pos, sf_tmp)
            sf_bam = site_rcd[0]
            ins_chrm = site_rcd[1]
            ins_pos = site_rcd[2]
            ori_wfolder = site_rcd[3]

            s_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
####Note, s_prefix need to keep the same as the file name in l_asm.py
####This should be improved later
            s_prefix = ori_wfolder + s_id + "_{0}".format(global_values.LRD_EXTND_LEN)
            l_prefix_clean.append(s_prefix)
            sf_asm = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
            sf_flank_fa = sflank_folder + "/" + "{0}{1}{2}_flanks.fa".format(ins_chrm, global_values.SEPERATOR, ins_pos)
            l_site_flanks.append(sf_flank_fa)
            if os.path.isfile(sf_asm) == False or os.path.isfile(sf_flank_fa) == False:
                continue

            sf_algnmt = s_prefix + global_values.LRD_ALNMT_SUFFIX
            # if os.path.isfile(sf_algnmt)==True:
            #     continue
            l_rcd.append((sf_asm, sf_flank_fa, sf_algnmt, ins_chrm, ins_pos))
            # print sf_asm, sf_flank_fa, sf_algnmt#####################################################################
        ####
        # each record in format: (sf_asm, sf_flanks, sf_algnmt)
        xctg.align_flanks_to_contig_2(l_rcd)

        # 3.3. parse out the insertion from asm alignments
        m_constructed, m_refined_pos = xctg.call_MEI_from_long_read_contig_flank_algnmt(l_rcd, sf_mei)

        # 4. clean the intermediate files
        self._clean_intermediate_files_in_parallel(l_prefix_clean, l_site_flanks)
        # clean the folder
        if not os.listdir(sflank_folder):
            os.rmdir(sflank_folder)  # only remove empty folder
        return m_constructed, m_refined_pos
####
####
    ######## 3. call out from asm
    # each rcd of l_sites in format: (sf_bam, ins_chrm, ins_pos, sf_tmp)
    def call_MEIs_from_region_asm_in_parallel(self, sf_ref, sf_sites, l_sites, s_working_folder, sf_mei):
        # 3.1. collect flanking regions
        # 3.2. align the flanking regions to the assembled contigs
        i_extend = abs(global_values.LRD_EXTND_LEN)  # there is a case of negative (to collect whole reads)
        xref = XReference()
        xref.gnrt_flank_region_for_regions(sf_sites, i_extend, self.n_jobs, sf_ref, s_working_folder)
        xctg = LRD_Contig(s_working_folder, self.n_jobs)

        sflank_folder = s_working_folder + global_values.FLANK_FOLDER
        l_rcd = []
        l_prefix_clean = []
        l_site_flanks = []

        for site_rcd in l_sites:
            # each in format: (sf_bam, ins_chrm, ins_pos, sf_tmp)
            sf_bam = site_rcd[0]
            ins_chrm = site_rcd[1]
            ins_pos = site_rcd[2]
            ori_wfolder = site_rcd[3]

            s_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
            ####Note, s_prefix need to keep the same as the file name in l_asm.py
            ####This should be improved later
            s_prefix = ori_wfolder + s_id + "_{0}".format(global_values.LRD_EXTND_LEN)
            l_prefix_clean.append(s_prefix)
            sf_asm = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
            sf_flank_fa = sflank_folder + "/" + "{0}{1}{2}_flanks.fa".format(ins_chrm, global_values.SEPERATOR,
                                                                             ins_pos)
            l_site_flanks.append(sf_flank_fa)
            if os.path.isfile(sf_asm) == False or os.path.isfile(sf_flank_fa) == False:
                continue

            sf_algnmt = s_prefix + global_values.LRD_ALNMT_SUFFIX
            # if os.path.isfile(sf_algnmt)==True:
            #     continue
            l_rcd.append((sf_asm, sf_flank_fa, sf_algnmt, ins_chrm, ins_pos))
            # print sf_asm, sf_flank_fa, sf_algnmt#####################################################################
        ####
        # each record in format: (sf_asm, sf_flanks, sf_algnmt)
        xctg.align_flanks_to_contig_2(l_rcd)

        # 3.3. parse out the insertion from asm alignments
        m_constructed, m_refined_pos = xctg.call_ref_copy_from_contig_flank_algnmt(l_rcd, sf_mei)

        # 4. clean the intermediate files
        self._clean_intermediate_files_in_parallel(l_prefix_clean, l_site_flanks)
        # clean the folder
        if not os.listdir(sflank_folder):
            os.rmdir(sflank_folder)  # only remove empty folder
        return m_constructed, m_refined_pos
####
####
    def _rname_id_to_sid(self, ilen):
        s_sid=0
        if ilen<0:
            s_sid="neg_"+str(abs(ilen))
        else:
            s_sid=str(abs(ilen))
        return s_sid

    ####
    def call_complex_MEIs_from_asm_in_parallel(self, sf_ref, sf_l_sites, sf_r_sites, s_working_folder, i_extend,
                                               f_flk_map_ratio, sf_l1_rslt, sf_sva_rslt, flk_lenth, sf_rep_folder,
                                               b_hg19, sf_out):
        #first generate the left flanks for left sites
        # and right flanks for right sites
        xref = XReference()
        xref.gnrt_flank_region_for_sites(sf_l_sites, i_extend, self.n_jobs, sf_ref, s_working_folder, True, False)
        xref.gnrt_flank_region_for_sites(sf_r_sites, i_extend, self.n_jobs, sf_ref, s_working_folder, False, True)

        xctg = Complex_TE_SV(s_working_folder, self.n_jobs)
        sf_asm_folder=s_working_folder+self.ASM_FOLDER+"/neg_{0}/".format(abs(i_extend))

        l_left_rcds=self._load_in_sites(sf_l_sites, sf_asm_folder)#each record in format (chrm, pos, sf_wfolder)
        l_right_rcds=self._load_in_sites(sf_r_sites, sf_asm_folder)

        sflank_folder = s_working_folder + global_values.FLANK_FOLDER
        #align the left (right) flank to itself's contig, to get the right (left) part
        l_rcd=[]
        for (chrm, pos, sf_wfolder) in (l_left_rcds+l_right_rcds):
            s_id = "{0}{1}{2}".format(chrm, global_values.SEPERATOR, pos)
            sf_flank_fa = sflank_folder + "/{0}_flanks.fa".format(s_id)
            s_prefix = sf_asm_folder + s_id + "_{0}".format(global_values.LRD_EXTND_LEN)
            sf_asm = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
            sf_algnmt = s_prefix + global_values.LRD_ALNMT_SUFFIX
            l_rcd.append((sf_asm, sf_flank_fa, sf_algnmt, chrm, pos))
        xctg.align_flanks_to_contig_2(l_rcd)

        # parse, merge the trimmed right (left) part to two separate files
        sf_l_trimmed=sf_asm_folder+"left_trimmed_contig.fa"#trimmed segments from left flank alignments
        sf_r_trimmed = sf_asm_folder + "right_trimmed_contig.fa"#trimmed segments from right flank alignments
        m_llen, m_rlen=xctg.parse_trim_seq_from_flank_contig_algnmt(l_rcd, f_flk_map_ratio, sf_l_trimmed, sf_r_trimmed)
####
        #align the left (right) flank to the right (left) merged file
        l_rcd2_left=[]
        for (chrm, pos, sf_wfolder) in l_left_rcds:
            s_id = "{0}{1}{2}".format(chrm, global_values.SEPERATOR, pos)
            sf_flank_fa = sflank_folder + "/{0}_flanks.fa".format(s_id)
            s_prefix = sf_asm_folder + s_id + "_{0}".format(global_values.LRD_EXTND_LEN)
            sf_algnmt = s_prefix + "_to_right_trimmed_"+global_values.LRD_ALNMT_SUFFIX
            l_rcd2_left.append((sf_r_trimmed, sf_flank_fa, sf_algnmt, chrm, pos))
        xctg.align_flanks_to_contig_2(l_rcd2_left)

        l_rcd2_right=[]
        for (chrm, pos, sf_wfolder) in l_right_rcds:
            s_id = "{0}{1}{2}".format(chrm, global_values.SEPERATOR, pos)
            sf_flank_fa = sflank_folder + "/{0}_flanks.fa".format(s_id)
            s_prefix = sf_asm_folder + s_id + "_{0}".format(global_values.LRD_EXTND_LEN)
            sf_algnmt = s_prefix + "_to_left_trimmed_"+global_values.LRD_ALNMT_SUFFIX
            l_rcd2_right.append((sf_l_trimmed, sf_flank_fa, sf_algnmt, chrm, pos))
        xctg.align_flanks_to_contig_2(l_rcd2_right)

        #parse the alignment and call out the candidates
        m_l_seq_rcds=xctg.parse_flk_2_trim_algnmt_for_internal_seq(l_rcd2_left, f_flk_map_ratio)
        m_r_seq_rcds=xctg.parse_flk_2_trim_algnmt_for_internal_seq(l_rcd2_right, f_flk_map_ratio)

########
#format: m_seqs[s_id] = (s_contig, i_seq_start, i_seq_end)
        for s_tmp_id in m_l_seq_rcds:
            print "Test1", s_tmp_id, m_l_seq_rcds[s_tmp_id][0], m_l_seq_rcds[s_tmp_id][1], m_l_seq_rcds[s_tmp_id][2]
        for s_tmp_id in m_r_seq_rcds:
            print "Test2", s_tmp_id, m_r_seq_rcds[s_tmp_id][0], m_r_seq_rcds[s_tmp_id][1], m_r_seq_rcds[s_tmp_id][2]
########

        #mask the middle region and classify based on alignment (segment orientation)
        #print m_l_seq_rcds, m_r_seq_rcds
        f_len_diff_ratio=0.25

        #note, m_simple_sv save candidates of "simple deletion" or "simple translocation"
        m_candidates, m_simple_sv=xctg.classify_complex_SV(m_l_seq_rcds, m_r_seq_rcds, sf_l_trimmed, sf_r_trimmed,
                                                           f_len_diff_ratio)

        #output for log usage: m_candidates[s_new_id]=(s_lchrm, s_lpos, s_rchrm, s_rpos, s_l_seq, s_r_seq)
        with open(s_working_folder+self.ASM_FOLDER+"/tmp_candidates.log", "w") as fout_tmp:
            for s_tmp_id in m_candidates:
                fout_tmp.write(s_tmp_id+"\t")
                (s_lchrm, s_lpos, s_rchrm, s_rpos, s_l_seq, s_r_seq)=m_candidates[s_tmp_id]
                fout_tmp.write(s_lchrm+"\t"+str(s_lpos)+"\t"+s_rchrm+"\t"+str(s_rpos)+"\t"+s_l_seq+"\t"+s_r_seq)
####
        f_segmt_map_ratio=f_flk_map_ratio
        f_cover_ratio=0.9 #not used so far

        #complex sv folder
        sf_complex_sv_folder=s_working_folder+self.COMPLEX_SV_FOLDER
        if os.path.exists(sf_complex_sv_folder)==False:
            try:
                os.mkdir(sf_complex_sv_folder)
            except OSError:
                print ("Creation of the directory %s failed" % sf_complex_sv_folder)#

        xctg.classify_events_by_seq_algnmt_mask(m_candidates, sf_ref, f_segmt_map_ratio, f_cover_ratio, sf_l1_rslt,
                                                sf_sva_rslt, flk_lenth, sf_complex_sv_folder, sf_rep_folder, b_hg19,
                                                m_simple_sv, sf_out)
####
####
    def _load_in_sites(self, sf_sites, sf_wfolder):
        l_sites=[]
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                l_sites.append((chrm, pos, sf_wfolder))
        return l_sites

####
    ####
    ####this version will align the assembled contigs to target reference seqs
    def call_MEIs_from_asm_contig_2_ref_in_parallel(self, sf_ref, sf_sites, l_sites, s_working_folder, sf_mei):
        # 3.1. collect flanking regions
        # 3.2. align the flanking regions to the assembled contigs
        ###here extend twice: 1. for extension of the extended reads 2. extension for inprecise breakpoint
        i_extend = global_values.LRD_EXTND_LEN + global_values.BRKPNT_CK_WIN
        xref = XReference()
        xref.gnrt_target_flank_seq_for_sites(sf_sites, i_extend, self.n_jobs, sf_ref, s_working_folder)
        xctg = LRD_Contig(s_working_folder, self.n_jobs)

        sflank_folder = s_working_folder + global_values.FLANK_FOLDER
        l_rcd = []
        l_prefix_clean = []
        l_site_flanks = []
        for site_rcd in l_sites:
            # each in format: (sf_bam, ins_chrm, ins_pos, sf_tmp)
            sf_bam = site_rcd[0]
            ins_chrm = site_rcd[1]
            ins_pos = site_rcd[2]
            ori_wfolder = site_rcd[3]

            s_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
            s_prefix = ori_wfolder + s_id
            l_prefix_clean.append(s_prefix)
            sf_asm = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
            sf_flank_fa = sflank_folder + "/" + "{0}{1}{2}_flanks.fa".format(ins_chrm, global_values.SEPERATOR,
                                                                             ins_pos)
            l_site_flanks.append(sf_flank_fa)
            if os.path.isfile(sf_asm) == False or os.path.isfile(sf_flank_fa) == False:
                continue

            sf_algnmt = ori_wfolder + s_id + global_values.LRD_ALNMT_SUFFIX
            # if os.path.isfile(sf_algnmt)==True:
            #     continue
            l_rcd.append((sf_asm, sf_flank_fa, sf_algnmt, ins_chrm, ins_pos))
            # print sf_asm, sf_flank_fa, sf_algnmt#######################################################################
            ####
        # each record in format: (sf_asm, sf_flanks, sf_algnmt)
        xctg.align_contig_to_target_ref(l_rcd)

        # 3.3. parse out the insertion from asm alignments
        m_constructed = xctg.call_MEI_from_long_read_contig_flank_algnmt(l_rcd, sf_mei)

        # 4. clean the intermediate files
        self._clean_intermediate_files_in_parallel(l_prefix_clean, l_site_flanks)
        # clean the folder
        if not os.listdir(sflank_folder):
            os.rmdir(sflank_folder)  # only remove empty folder
        return m_constructed
#
####
    # clean the files
    def clean_intermediate_files(self, sf_sites):
        sf_tmp = self.wfolder + self.ASM_FOLDER+"/"

        l_sites = []
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields = line.split()

                ins_chrm = fields[0]
                ins_pos = int(fields[1])
                rcd = (ins_chrm, ins_pos, sf_tmp)
                l_sites.append(rcd)
        self.clean_intermediate_files_one_bam(l_sites, sf_tmp)

    def clean_intermediate_files2(self, l_extnd, sf_sites):
        for iextnd in l_extnd:
            sf_tmp = self.wfolder + self.ASM_FOLDER+"/{0}/".format(iextnd)
            if iextnd<0:
                sf_tmp = self.wfolder + self.ASM_FOLDER+"/neg_{0}/".format(abs(iextnd))

            l_sites = []
            with open(sf_sites) as fin_sites:
                for line in fin_sites:
                    fields = line.split()

                    ins_chrm = fields[0]
                    ins_pos = int(fields[1])
                    rcd = (ins_chrm, ins_pos, sf_tmp, iextnd)
                    l_sites.append(rcd)
            self.clean_intermediate_files_one_bam(l_sites, sf_tmp)

    ####Only clean the intermediate files
    def clean_intermediate_files_one_bam(self, l_sites, s_working_folder):
        sflank_folder = s_working_folder + global_values.FLANK_FOLDER
        l_prefix_clean = []
        l_site_flanks = []
        for site_rcd in l_sites:
            ins_chrm = site_rcd[0]
            ins_pos = site_rcd[1]
            ori_wfolder = site_rcd[2]
            iextnd=site_rcd[3]
            s_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
            s_prefix = ori_wfolder + s_id + "_" + str(iextnd)
            l_prefix_clean.append(s_prefix)
            sf_flank_fa = sflank_folder + "/" + "{0}_{1}_flanks.fa".format(ins_chrm, ins_pos)
            #print sf_flank_fa, "flank test"
            l_site_flanks.append(sf_flank_fa)
####
        self._clean_intermediate_files_in_parallel(l_prefix_clean, l_site_flanks)
        # clean the folder
        if (os.path.exists(sflank_folder)) and (not os.listdir(sflank_folder)):
            os.rmdir(sflank_folder)  # only remove empty folder
####
####
    # remove the intermeidate files in parallel
    def _clean_intermediate_files_in_parallel(self, l_prefix, l_flanks):
        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_remove_intermediate_files, zip([self] * len(l_prefix), l_prefix), 1)
        pool.close()
        pool.join()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_remove_intermediate_flank_files, zip([self] * len(l_flanks), l_flanks), 1)
        pool.close()
        pool.join()
####
####
    # remove one group of intermediate files
    def clean_intermediate_output_one_sample(self, sf_prefix):
        sf_asm_1_reads = sf_prefix + "_wtdbg2.1.reads"
        self.___remove_one_file(sf_asm_1_reads)
        sf_asm_1_dot = sf_prefix + "_wtdbg2.1.dot.gz"
        self.___remove_one_file(sf_asm_1_dot)
        sf_asm_1_nodes = sf_prefix + "_wtdbg2.1.nodes"
        self.___remove_one_file(sf_asm_1_nodes)
        sf_asm_2_dot = sf_prefix + "_wtdbg2.2.dot.gz"
        self.___remove_one_file(sf_asm_2_dot)
        sf_asm_3_dot = sf_prefix + "_wtdbg2.3.dot.gz"
        self.___remove_one_file(sf_asm_3_dot)
        sf_algnmts = sf_prefix + "_wtdbg2.alignments.gz"
        self.___remove_one_file(sf_algnmts)
        sf_binkmer = sf_prefix + "_wtdbg2.binkmer"
        self.___remove_one_file(sf_binkmer)
        sf_bins = sf_prefix + "_wtdbg2.closed_bins"
        self.___remove_one_file(sf_bins)
        sf_clps = sf_prefix + "_wtdbg2.clps"
        self.___remove_one_file(sf_clps)
        sf_ctg_dot = sf_prefix + "_wtdbg2.ctg.dot.gz"
        self.___remove_one_file(sf_ctg_dot)
        sf_ctg_lay = sf_prefix + "_wtdbg2.ctg.lay.gz"
        self.___remove_one_file(sf_ctg_lay)
        # sf_ctg_lay_fa = sf_prefix + "_wtdbg2.ctg.lay.fa"
        sf_events = sf_prefix + "_wtdbg2.events"
        self.___remove_one_file(sf_events)
        sf_frg_dot = sf_prefix + "_wtdbg2.frg.dot.gz"
        self.___remove_one_file(sf_frg_dot)
        sf_frg_nodes = sf_prefix + "_wtdbg2.frg.nodes"
        self.___remove_one_file(sf_frg_nodes)
        sf_kmerdep = sf_prefix + "_wtdbg2.kmerdep"
        self.___remove_one_file(sf_kmerdep)
        sf_std_out = sf_prefix + ".fa.std_out"
        self.___remove_one_file(sf_std_out)

        # remove the flanking and asmbled sequence
        sf_ctg_lay_fa = sf_prefix + "_wtdbg2.ctg.lay.fa"
        #self.___remove_one_file(sf_ctg_lay_fa)
        sf_ctg_lay_fai = sf_prefix + "_wtdbg2.ctg.lay.fa.fai"
        self.___remove_one_file(sf_ctg_lay_fai)
        #align flanks to contigs bam files
        sf_flank_2_ctg_bam=sf_prefix+"_aln_flank_2_ctg.bam"
        self.___remove_one_file(sf_flank_2_ctg_bam)
        sf_flank_2_ctg_bai = sf_prefix + "_aln_flank_2_ctg.bam.bai"
        self.___remove_one_file(sf_flank_2_ctg_bai)


####

    def clean_intermediate_output_flanks_one_sample(self, sf_flank):
        self.___remove_one_file(sf_flank)

    def ___remove_one_file(self, sf_file):
        if os.path.isfile(sf_file):
            os.remove(sf_file)

    def get_ASM_folder(self):
        return self.ASM_FOLDER
####
####

class CallFromRMSK():
    def _load_in_sites(self, sf_sites):
        l_sites = []
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields = line.split()

                ins_chrm = fields[0]
                ins_pos = int(fields[1])
                rcd = (ins_chrm, ins_pos)
                l_sites.append(rcd)
        return l_sites

    # if the region is overlap, then merg
    def _get_max_cover(self, l_regions):
        l_max_cov = []
        cur_start = -1
        cur_end = -1
        l_sorted_regions = sorted(l_regions, key=lambda tup: tup[0])
        for (istart, iend) in l_sorted_regions:
            if cur_start == -1:
                cur_start = istart
                cur_end = iend
            elif istart >= cur_start and istart <= cur_end:  # two region has overlap, then merge
                if iend > cur_end:
                    cur_end = iend
            else:  # non overlap
                l_max_cov.append((cur_start, cur_end))
                cur_start = istart
                cur_end = iend
        if cur_start != -1:
            l_max_cov.append((cur_start, cur_end))
        return l_max_cov

    def _split_by_rc(self, l_hits, l_type):
        l_rc = []
        l_not_rc = []
        for rcd in l_hits:
            i_contig_start = rcd[0]
            i_contig_end = rcd[1]
            b_rc = rcd[2]
            s_sub_family = rcd[3]
            s_rep_family = rcd[4]
            i_cns_start = rcd[5]
            i_cns_end = rcd[6]

            b_hit = False
            for s_type in l_type:
                if s_type in s_sub_family:
                    b_hit = True
                    break
            if b_hit == False:
                continue

            rcd_on_cns = (i_cns_start, i_cns_end)
            if b_rc == True:
                l_rc.append(rcd_on_cns)
            else:
                l_not_rc.append(rcd_on_cns)
        return l_rc, l_not_rc

    ####

    # hit_rcd=(i_contig_start, i_contig_end, b_rc, s_sub_family, s_rep_family, i_cns_start, i_cns_end)
    # here we need to consider inversion for L1, and also the expansion for SVA
    def slct_L1_hits(self, l_hits, l_type):
        l_rc, l_not_rc = self._split_by_rc(l_hits, l_type)

        l_rc_max = self._get_max_cover(l_rc)
        l_non_rc_max = self._get_max_cover(l_not_rc)

        # this is an inversion case
        if len(l_rc_max) > 0 and len(l_non_rc_max) > 0:
            return (l_rc_max, l_non_rc_max)

        elif len(l_rc_max) > 0:
            return (l_rc_max, None)
        elif len(l_non_rc_max) > 0:
            return (None, l_non_rc_max)
        else:
            return (None, None)
            #
            # if some region cover another one, and of the same sub-family, then they are merged
            #

            # hit_rcd=(i_contig_start, i_contig_end, b_rc, s_sub_family, s_rep_family, i_cns_start, i_cns_end)

    #
    def _split_by_rc_SVA(self, l_hits, l_type):
        l_rc = []
        l_not_rc = []
        for rcd in l_hits:
            i_contig_start = rcd[0]
            i_contig_end = rcd[1]
            b_rc = rcd[2]
            s_sub_family = rcd[3]
            s_rep_family = rcd[4]
            i_cns_start = rcd[5]
            i_cns_end = rcd[6]

            b_hit = False
            for s_type in l_type:
                if s_type in s_sub_family:
                    b_hit = True
                    break
            if b_hit == False:
                continue

            rcd_on_ctg = (i_contig_start, i_contig_end)
            if b_rc == True:
                l_rc.append(rcd_on_ctg)
            else:
                l_not_rc.append(rcd_on_ctg)
        return l_rc, l_not_rc

    # As SVA is much more sparse distributed, so we count how many bases of the contig is masked.
    # Here l_hits is the list of regions of one contig
    def slct_SVA_hits(self, l_hits, s_type):
        print "slct_SVA_hits starts..."
        l_rc, l_not_rc = self._split_by_rc_SVA(l_hits, s_type)
        l_rc_max = self._get_max_cover(l_rc)
        l_non_rc_max = self._get_max_cover(l_not_rc)

        # this is an inversion case
        if len(l_rc_max) > 0 and len(l_non_rc_max) > 0:
            return (l_rc_max, l_non_rc_max)

        elif len(l_rc_max) > 0:
            return (l_rc_max, None)
        elif len(l_non_rc_max) > 0:
            return (None, l_non_rc_max)
        else:
            return (None, None)

    def _calc_cov_length(self, l_covs):
        if l_covs is None:
            return 0, -1, -1

        i_cov_len = 0
        i_left = -1
        i_right = -1
        for region in l_covs:
            istart = region[0]
            iend = region[1]

            if i_left == -1:
                i_left = istart
                i_right = iend
            else:
                if i_left > istart:
                    i_left = istart
                if i_right < iend:
                    i_right = iend
            i_cov_len += (iend - istart)
        return i_cov_len, i_left, i_right

    ####
    def call_TEI_from_rmsk_L1(self, sf_rmsk, sf_sites, sf_out):
        rmsk_parser = RMSK_Parser(sf_rmsk)
        # m_hits in format: {chrm:{pos:[]}}
        # m_hits_by_contig in format: {contig-id:[]}
        m_hits, m_hits_by_contig = rmsk_parser.parse_rmsk()
        l_sites = self._load_in_sites(sf_sites)
        l_L1_family = ["L1HS", "L1P"]  # we only consider these two type of LINE1
        with open(sf_out, "w") as fout_rslt:
            for rcd in l_sites:
                ins_chrm = rcd[0]
                ins_pos = int(rcd[1])

                s_inv = "not_inversion"
                if (ins_chrm in m_hits) and (ins_pos in m_hits[ins_chrm]):
                    l_hits = m_hits[ins_chrm][ins_pos]
                    l_rc_max, l_non_rc_max = self.slct_L1_hits(l_hits, l_L1_family)

                    i_len_rc, i_rc_left, i_rc_right = self._calc_cov_length(l_rc_max)
                    i_len_nrc, i_nrc_left, i_nrc_right = self._calc_cov_length(l_non_rc_max)

                    sinfo = "{0}\t{1}\t".format(ins_chrm, ins_pos)
                    sinfo1 = "{0}\t{1}\t{2}\t".format(i_rc_left, i_rc_right, i_len_rc)
                    sinfo2 = "{0}\t{1}\t{2}\t".format(i_nrc_left, i_nrc_right, i_len_nrc)
                    all_len = 0
                    if i_len_rc > 0 and i_len_nrc > 0:
                        b_inv, all_len = self._check_inversion(i_rc_left, i_rc_right, i_nrc_left, i_nrc_right)
                        if b_inv == True:
                            s_inv = "inversion"
                    else:
                        all_len = i_len_rc + i_len_nrc
                    sinfo3 = "{0}\t{1}\n".format(s_inv, all_len)
                    fout_rslt.write(sinfo + sinfo1 + sinfo2 + sinfo3)
                #
                ####

    def _check_inversion(self, i_rc_left, i_rc_right, i_nrc_left, i_nrc_right):
        all_len = 0
        b_inv = False
        l_tmp = []
        l_tmp.append((i_rc_left, i_rc_right))
        l_tmp.append((i_nrc_left, i_nrc_right))
        l_tmp_cov = self._get_max_cover(l_tmp)
        for tmp_rcd in l_tmp_cov:
            all_len += (tmp_rcd[1] - tmp_rcd[0])
        if i_nrc_right > i_rc_right and abs(i_nrc_left - i_rc_right) < 100:
            b_inv = True
        elif i_rc_right > i_nrc_right and abs(i_rc_left - i_nrc_right) < 100:
            b_inv = True

        # if contained, then this is not inversion
        if i_nrc_left < i_rc_left and i_nrc_right > i_rc_right:
            b_inv = False
        if i_rc_left < i_nrc_left and i_rc_right > i_nrc_right:
            b_inv = False

        return b_inv, all_len

    ####
    def call_TEI_from_rmsk_SVA(self, sf_rmsk, sf_sites, sf_out):
        rmsk_parser = RMSK_Parser(sf_rmsk)
        # m_hits in format: {chrm:{pos:[]}}
        # m_hits_by_contig in format: {contig-id:[]}
        m_hits, m_hits_by_contig = rmsk_parser.parse_rmsk()
        l_sites = self._load_in_sites(sf_sites)
        l_SVA_family = ["SVA"]

        with open(sf_out, "w") as fout_rslt:
            for rcd in l_sites:
                ins_chrm = rcd[0]
                ins_pos = int(rcd[1])
                if (ins_chrm in m_hits_by_contig) and (ins_pos in m_hits_by_contig[ins_chrm]):
                    i_max_hit_len = 0
                    s_info_max = ""
                    for s_contig in m_hits_by_contig[ins_chrm][ins_pos]:
                        l_hits = m_hits_by_contig[ins_chrm][ins_pos][s_contig]
                        if l_hits is None:
                            continue
                        l_rc_max, l_non_rc_max = self.slct_SVA_hits(l_hits, l_SVA_family)  #
                        print l_rc_max, l_non_rc_max, "info"
                        s_inv = "not_inversion"
                        i_len_rc, i_rc_left, i_rc_right = self._calc_cov_length(l_rc_max)
                        i_len_nrc, i_nrc_left, i_nrc_right = self._calc_cov_length(l_non_rc_max)

                        sinfo = "{0}\t{1}\t".format(ins_chrm, ins_pos)
                        sinfo1 = "{0}\t{1}\t{2}\t".format(i_rc_left, i_rc_right, i_len_rc)
                        sinfo2 = "{0}\t{1}\t{2}\t".format(i_nrc_left, i_nrc_right, i_len_nrc)
                        all_len = 0
                        if i_len_rc > 0 and i_len_nrc > 0:
                            b_inv, all_len = self._check_inversion(i_rc_left, i_rc_right, i_nrc_left, i_nrc_right)
                            if b_inv == True:
                                s_inv = "inversion"
                        else:
                            all_len = i_len_rc + i_len_nrc
                        sinfo3 = "{0}\t{1}\n".format(s_inv, all_len)

                        if i_max_hit_len < all_len:
                            s_info_max = sinfo + sinfo1 + sinfo2 + sinfo3
                            i_max_hit_len = all_len
                    fout_rslt.write(s_info_max)
