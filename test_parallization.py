##11/14/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import sys
import subprocess
import global_values
import select
#from multiprocessing import Pool
import multiprocessing as mp
#To-do list:
#1. Need to check the cases contained in the cigar!!!!!

#Given a list of candidates, check whether they can be assembled from long reads and repeatmasker the region
#1. wtdbg2 seems work better than minimap2
#2. Performance varied for different bams aligned from different tool


class CMD_RUNNER():
    ####This is assume the output is small amount of data
    ####Note: one issue for Popen.communicate() is: when output file is large, the process will be hanged, which will
    ####cause deadlock
    def run_cmd_small_output(self, cmd):
        print "Running command: {0}\n".format(cmd)
        subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()

    ####This is try to solve the issue mentioned previously with assumption that stderr is small
    def run_cmd_to_file(self, cmd, sf_out):
        print "Running command with output: {0}\n".format(cmd)
        sf_err=sf_out+".err"
        errcode=None
        with open(sf_out, "w") as f:
            p = subprocess.Popen(cmd, shell=True, stdout=f, stderr=None)
            #errcode = p.wait()
            p.communicate()

        # if errcode:
        #     errmess = p.stderr.read()
        #     log.error('cmd failed <%s>: %s' % (errcode, errmess))
        if os.path.isfile(sf_err):
            os.remove(sf_err)
        return errcode
####

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
    def construct_short_cns_wtdbg2(self, chrm, pos, wfolder):
        if wfolder[-1] != "/":
            wfolder += "/"
        s_id = "{0}_{1}".format(chrm, pos)
        s_prefix = wfolder + s_id
        sf_fa = s_prefix + ".fa"
        if os.path.isfile(sf_fa) == False:
            return global_values.FAIL_CLP
        cmd = "{0}  -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 --ctg-min-nodes 1 " \
              "-i {1}.fa  -fo {2}_wtdbg2 -q".format(global_values.WTDBG2, s_prefix, s_prefix)
        print cmd
        sf_std_out=sf_fa+".std_out"
        self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)

        sf_lay = "{0}_wtdbg2.ctg.lay".format(s_prefix)
        if os.path.isfile(sf_lay) == False:
            return global_values.FAIL_ASM
        cmd = "{0} -i {1}_wtdbg2.ctg.lay -fo {2}_wtdbg2.ctg.lay.fa".format(global_values.WTPOA, s_prefix, s_prefix)
        print cmd
        self.cmd_runner.run_cmd_to_file(cmd, sf_std_out)
        sf_lay_fa = "{0}_wtdbg2.ctg.lay.fa".format(s_prefix)
        if os.path.isfile(sf_lay_fa) == False:
            return global_values.FAIL_CNS
        return global_values.SUCCD

####

######## 2. the assembly step will fail if run in parallel
def asm_seq_for_sites_in_serial(sf_sites, n_jobs, wfolder):
    l_sites=[]
    with open(sf_sites) as fin_sites:
        for line in fin_sites:
            fields=line.split()
            chrm=fields[0]
            pos=fields[1]
            l_sites.append(("none", chrm, pos, wfolder))

    pool = mp.get_context("spawn").Pool(processes = n_jobs)
    pool.map(_asm_collected_reads_one_site, l_sites, 1)
    pool.close()
    pool.join()


# assemble the collected clipped and contained reads
def _asm_collected_reads_one_site(record):
    sf_bam = record[0]
    ins_chrm = record[1]
    ins_pos = record[2]
    wfolder = record[3]
    lasm = L_Local_ASM()
    lasm.construct_short_cns_wtdbg2(ins_chrm, ins_pos, wfolder)

if __name__ == '__main__':
    sf_sites=sys.argv[1]
    n_jobs=int(sys.argv[2])
    wfolder=sys.argv[3]
    asm_seq_for_sites_in_serial(sf_sites, n_jobs, wfolder)

