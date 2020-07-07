##11/14/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import sys
import pysam
from subprocess import *
import global_values
from cmd_runner import *

#To-do list:
#1. Need to check the cases contained in the cigar!!!!!

#Given a list of candidates, check whether they can be assembled from long reads and repeatmasker the region
#1. wtdbg2 seems work better than minimap2
#2. Performance varied for different bams aligned from different tool

####
class L_Local_ASM():
    def __init__(self):
        self.cmd_runner=CMD_RUNNER()

    #this is to use wtdbg2 to asm the clipped long reads
    def construct_cns_wtdbg2(self, chrm, pos, wfolder):
        if wfolder[-1] != "/":
            wfolder += "/"
        s_id = "{0}_{1}".format(chrm, pos)
        s_prefix = wfolder + s_id

        sf_fa=s_prefix+".fa"
        if os.path.isfile(sf_fa)==False:
            return global_values.FAIL_CLP
        sf_std_out = sf_fa + ".std_out"
        cmd = "{0} -i {1}.fa  -fo {2}_wtdbg2".format(global_values.WTDBG2, s_prefix, s_prefix)
        print cmd
        self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)

        ####
        sf_lay="{0}_wtdbg2.ctg.lay".format(s_prefix)
        if os.path.isfile(sf_lay)==False:
            return global_values.FAIL_ASM
        cmd = "{0} -i {1}_wtdbg2.ctg.lay -fo {2}_wtdbg2.ctg.lay.fa".format(global_values.WTPOA, s_prefix, s_prefix)
        print cmd
        self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)
        sf_lay_fa="{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
        if os.path.isfile(sf_lay_fa)==False:
            return global_values.FAIL_CNS
        return global_values.SUCCD

    ##this is to construct the consensus
    def construct_short_cns_wtdbg2(self, chrm, pos, ncores, wfolder):
        if wfolder[-1] != "/":
            wfolder += "/"
        s_id = "{0}{1}{2}".format(chrm, global_values.SEPERATOR, pos)
        i_extend = global_values.LRD_EXTND_LEN
        s_prefix = wfolder + s_id + "_{0}".format(i_extend)
        sf_fa = s_prefix + ".fa"
        if os.path.isfile(sf_fa) == False:
            return global_values.FAIL_CLP
        cmd = "{0}  -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 --ctg-min-nodes 1 " \
              "-i {1}.fa  -fo {2}_wtdbg2 -q -t {3}".format(global_values.WTDBG2, s_prefix, s_prefix, ncores)
        print cmd
        sf_std_out=sf_fa+".std_out"
        self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)

        sf_lay = "{0}_wtdbg2.ctg.lay".format(s_prefix)
        cmd=""
        if os.path.isfile(sf_lay) == False:
            if os.path.isfile(sf_lay+".gz") == True:
                cmd = "{0} -i {1}_wtdbg2.ctg.lay.gz -fo {2}_wtdbg2.ctg.lay.fa -t {3}".format(global_values.WTPOA,
                                                                                             s_prefix, s_prefix, ncores)
            else:
                return global_values.FAIL_ASM
        else:
            cmd = "{0} -i {1}_wtdbg2.ctg.lay -fo {2}_wtdbg2.ctg.lay.fa -t {3}".format(global_values.WTPOA, s_prefix,
                                                                                      s_prefix, ncores)
        print cmd
        if cmd!="":
            self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)
        sf_lay_fa = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
        if os.path.isfile(sf_lay_fa) == False:
            return global_values.FAIL_CNS
        return global_values.SUCCD

    def get_wtdbg2_asm_contig_file(self, chrm, pos, wfolder):
        if wfolder[-1] != "/":
            wfolder += "/"
        s_id = "{0}{1}{2}".format(chrm, global_values.SEPERATOR, pos)
        i_extend = global_values.LRD_EXTND_LEN
        s_prefix = wfolder + s_id + "_{0}".format(i_extend)
        sf_lay_fa = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
        return sf_lay_fa

    ##this is to construct the consensus
    def construct_short_cns_wtdbg2_cmd_list(self, chrm, pos, ncores, wfolder):
        if wfolder[-1] != "/":
            wfolder += "/"
        s_id = "{0}{1}{2}".format(chrm, global_values.SEPERATOR, pos)
        i_extend = global_values.LRD_EXTND_LEN
        s_prefix = wfolder + s_id + "_{0}".format(i_extend)
        cmd1 = "{0}  -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 --ctg-min-nodes 1 " \
              "-i {1}.fa  -fo {2}_wtdbg2 -q -t {3}".format(global_values.WTDBG2, s_prefix, s_prefix, ncores)
        print cmd1
        cmd2 = "{0} -i {1}_wtdbg2.ctg.lay.gz -fo {2}_wtdbg2.ctg.lay.fa -t {3}".format(global_values.WTPOA,
                                                                                     s_prefix, s_prefix, ncores)
        print cmd2

        return cmd1, cmd2

        ####

    ##this is to construct the consensus
    def construct_short_cns_wtdbg2_given_file(self, sf_fa, ncores, wfolder):
        if wfolder[-1] != "/":
            wfolder += "/"

        if os.path.isfile(sf_fa) == False:
            return global_values.FAIL_CLP
        cmd = "{0}  -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 --ctg-min-nodes 1 " \
              "-i {1}  -fo {2}_wtdbg2 -q -t {3}".format(global_values.WTDBG2, sf_fa, sf_fa, ncores)
        print cmd
        sf_std_out = sf_fa + ".std_out"
        self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)

        sf_lay = "{0}_wtdbg2.ctg.lay".format(sf_fa)
        cmd = ""
        if os.path.isfile(sf_lay) == False:
            if os.path.isfile(sf_lay + ".gz") == True:
                cmd = "{0} -i {1}_wtdbg2.ctg.lay.gz -fo {2}_wtdbg2.ctg.lay.fa -t {3}".format(global_values.WTPOA,
                                                                                             sf_fa, sf_fa, ncores)
            else:
                return global_values.FAIL_ASM
        else:
            cmd = "{0} -i {1}_wtdbg2.ctg.lay -fo {2}_wtdbg2.ctg.lay.fa -t {3}".format(global_values.WTPOA, sf_fa,
                                                                                      sf_fa, ncores)
        print cmd
        if cmd!="":
            self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)
        sf_lay_fa = "{0}_wtdbg2.ctg.lay.fa".format(sf_fa)
        if os.path.isfile(sf_lay_fa) == False:
            return global_values.FAIL_CNS
        return global_values.SUCCD

####
    #this is use minimap2 to assemble the clipped long reads, and then use racon to polish the assembled contigs
    def construct_cns_minimap2(self, chrm, pos, wfolder):
        if wfolder[-1] != "/":
            wfolder += "/"
        s_id = "{0}_{1}".format(chrm, pos)
        sfa_clip = wfolder + "{0}.fa".format(s_id)
        sf_paf = wfolder + "{0}.paf".format(s_id)
        sf_gfa = wfolder + "{0}.gfa".format(s_id)
        sf_utg = wfolder + "{0}.utg".format(s_id)
        sfa_utg = wfolder + "{0}.utg.fa".format(s_id)
        sf_paf2 = wfolder + "{0}_2_utg.paf".format(s_id)
        sf_cns = wfolder + "{0}.utg.cns.fa".format(s_id)

        cmd = "{0} -Hk19 -Xw3 -m40 -g10000 --max-chain-skip 25 {1} {2} > {3}".format(global_values.MINIMAP2, sfa_clip,
                                                                                     sfa_clip, sf_paf)
        self.cmd_runner.run_cmd_small_output(cmd)
        cmd = "{0} -s100 -S4 -c1 -f {1} {2} 1> {3}".format(global_values.MINIMAP2, sfa_clip, sf_paf, sf_gfa)
        self.cmd_runner.run_cmd_small_output(cmd)

        with open(sf_gfa) as fin_gfa, open(sfa_utg, "w") as fout_utg_fa:
            for line in fin_gfa:
                fields = line.split()
                stype = fields[0]
                if stype != "S":
                    continue
                uid = fields[1]
                s_seq = fields[2]
                fout_utg_fa.write(">" + uid + "\n")
                fout_utg_fa.write(s_seq + "\n")

        cmd = "{0} -Hk19 -Xw3 -m40 -g10000 --max-chain-skip 25 {1} {2} > {3}".format(global_values.MINIMAP2, sfa_utg,
                                                                                     sfa_clip, sf_paf2)
        self.cmd_runner.run_cmd_small_output(cmd)
        cmd = "{0} {1} {2} {3} 1>{4}".format(global_values.RACON, sfa_clip, sf_paf2, sfa_utg, sf_cns)
        self.cmd_runner.run_cmd_small_output(cmd)

####
    def is_wtdbg2_asm_already_exist(self, ins_chrm, ins_pos, sf_folder):
        s_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
        ####Note, s_prefix need to keep the same as the file name in l_asm.py
        ####This should be improved later
        s_prefix = sf_folder + s_id + "_{0}".format(global_values.LRD_EXTND_LEN)
        sf_asm = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
        if os.path.isfile(sf_asm)==True:
            return True
        return False
