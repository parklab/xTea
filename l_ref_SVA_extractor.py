##12/27/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu, chong.simon.chu@gmail.com

import os
import pysam
import global_values
from l_MEI_caller import *
####
class RefSVAExtractor():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1] != "/":
            self.swfolder += "/"
        self.n_jobs = n_jobs

    ####
    def collect_asm_extract_ref_sva_seq(self, sf_sites, sf_ref, sf_bam_list, b_asm_cmd, sf_out_fa):
        #first load in sites
        xchrm=XChromosome()
        l_sites=[]
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields=line.rstrip().split()
                if len(fields)<3:
                    print("Wrong line %s" % line)
                    continue
                chrm=fields[0]#
                if xchrm.is_decoy_contig_chrms(chrm)==True:
                    continue
                istart=int(fields[1])
                iend=int(fields[2])
                l_sites.append((chrm, istart, iend))

        lcaller = L_MEI_Caller(self.swfolder, self.n_jobs, sf_ref)
        sf_tmp_prefix = self.swfolder + lcaller.get_ASM_folder()
        if os.path.exists(sf_tmp_prefix) == False:
            cmd = "mkdir {0}".format(sf_tmp_prefix)
            cmd_runner = CMD_RUNNER()
            cmd_runner.run_cmd_small_output(cmd)

        l_bams = lcaller._load_in_bam(sf_bam_list)
        if os.path.isfile(sf_out_fa) == True:  # as everytime each called seq is save in "append" way
            os.remove(sf_out_fa)

        m_constructed = {}
        m_new_pos = {}
        n_round = 1
        i_ext_len=-1000
        for sf_bam in l_bams:  # run for different bams
            sf_tmp = sf_tmp_prefix + "/{0}/".format(i_ext_len)
            if i_ext_len < 0:
                sf_tmp = sf_tmp_prefix + "/neg_{0}/".format(abs(i_ext_len))
            if os.path.exists(sf_tmp) == False:
                cmd = "mkdir {0}".format(sf_tmp)
                cmd_runner = CMD_RUNNER()
                cmd_runner.run_cmd_small_output(cmd)

            sf_tmp_to_asm = sf_out_fa + ".{0}_to_asm_list.txt".format(i_ext_len)
            l_new_sites=[]
            l_region_start=[]#
            with open(sf_tmp_to_asm, "w") as fout_for_asm:
                for old_rcd in l_sites:
                    ins_chrm = old_rcd[0]
                    i_copy_start = int(old_rcd[1])
                    i_copy_end=int(old_rcd[2])

                    rcd = (sf_bam, ins_chrm, i_copy_start, i_copy_end, sf_tmp)
                    l_new_sites.append(rcd)
                    fout_for_asm.write(ins_chrm+"\t"+str(i_copy_start)+"\t"+str(i_copy_end)+"\n")
                    rcd2 = (sf_bam, ins_chrm, i_copy_start, sf_tmp)
                    l_region_start.append(rcd2)
            global_values.set_lrd_extnd_len(i_ext_len)
            if len(l_new_sites) <= 0:
                continue
            # collect the clipped and contained reads for sites
####
            lcaller.collect_seq_cover_region_in_parallel(l_new_sites)
            # if l_cluster_info is None:#
            if b_asm_cmd == True:#
                sf_script = sf_out_fa + ".asm_{0}.sh".format(i_ext_len)
                lcaller.gnrt_asm_cmd_only(l_region_start, sf_script)
                continue
            else:
                lcaller.asm_seq_for_sites_in_serial(l_region_start)

            m_tmp, m_refined_pos = lcaller.call_MEIs_from_region_asm_in_parallel(sf_ref, sf_tmp_to_asm,
                                                                                 l_region_start, sf_tmp, sf_out_fa)
            ####
            for s_site in m_tmp:
                s_tmp_fields = s_site.split(global_values.SEPERATOR)
                if s_tmp_fields[0] not in m_constructed:
                    m_constructed[s_tmp_fields[0]] = {}
                m_constructed[s_tmp_fields[0]][s_tmp_fields[1]] = 1
                if s_tmp_fields[0] not in m_new_pos:
                    m_new_pos[s_tmp_fields[0]] = {}
                m_new_pos[s_tmp_fields[0]][s_tmp_fields[1]] = m_refined_pos[s_site]
            print "Round {0}: {1} sites (out of {2}) are constructed\n".format(n_round, len(m_tmp), len(l_new_sites))
            n_round += 1

        sf_out_sites=sf_out_fa+".constructed_ref_copy_sites"
        with open(sf_out_sites, "w") as fout_sites:
            for chrm in m_constructed:
                for pos in m_constructed[chrm]:
                    new_pos = m_new_pos[chrm][pos]
                    s_site_info = "{0}\t{1}\n".format(chrm, new_pos)
                    fout_sites.write(s_site_info)
####