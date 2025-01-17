##07/14/2021
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com
##This is a class for assemble the MEIs from paired-end short reads

#1. collect the clipped reads and discordant reads
#2. assemble the collected reads
#3. parse out the insertion sequence from the assembled contigs

#Input: this is assume given the candidate insertions, and the bam file;
#Output: then, output the assembled insertions

import os
import pysam
import global_values
from subprocess import *
from multiprocessing import Pool
from x_clip_disc_filter import *
from x_local_assembly import *

def unwrap_self_asm_clip_disc_reads_MEIs(arg, **kwarg):
    return ShortRead_MEI_ASM.asm_clip_disc_reads_one_MEI(*arg, **kwarg)

class ShortRead_MEI_ASM():####
    def __init__(self, sf_input, sf_bam_list, sf_wfolder, n_jobs):#
        self.sf_candidate_list=sf_input
        self.sf_bam_list=sf_bam_list
        self.sf_working_folder=sf_wfolder
        if "/" != self.sf_working_folder[-1]:
            self.sf_working_folder += "/"
        self.n_jobs=n_jobs

    #Here, sf_bam_ref need to be provided if it is cram files
    def asm_MEI_PE_reads(self, sf_bam_ref, extnd, bin_size, sf_out_fa):####
        x_cd_filter = XClipDiscFilter(self.sf_bam_list, self.sf_working_folder, self.n_jobs, sf_bam_ref)
        ####
        # collect the clipped and discordant reads
        # each record in format like: @20~41951715~L~1~41952104~0~0,
        # (chrm, map_pos, global_values.FLAG_LEFT_CLIP, is_rc, insertion_pos, n_cnt_clip, sample_id)
        sf_clip_fq = self.sf_working_folder + "candidate_sites_all_clip.fq"
        sf_disc_fa = self.sf_working_folder + "candidate_sites_all_disc.fa"
        x_cd_filter.collect_clipped_disc_reads(self.sf_candidate_list, extnd, bin_size, sf_clip_fq, sf_disc_fa)##

        #split the reads to small files, and then assemble each one
        m_sid_fa_folder= x_cd_filter.split_merged_reads_by_site(sf_clip_fq, sf_disc_fa,
                                                                                  self.sf_working_folder)
        l_site_reads=[]
        l_site_asm_folders=[]
        for s_id in m_sid_fa_folder:
            (sf_tmp_fa, s_site_folder)=m_sid_fa_folder[s_id]
            l_site_reads.append((s_site_folder, sf_tmp_fa))
            l_site_asm_folders.append((s_id, s_site_folder))

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_asm_clip_disc_reads_MEIs,
                              list(zip([self] * len(l_site_reads), l_site_reads)), 1)
        pool.close()
        pool.join()

        with open(sf_out_fa,"w") as fout_fa:
            #select and merge the results to one file
            for (s_id, sf_site_folder) in l_site_asm_folders:
                s_ins_seq=self._get_asm_sequence(s_id, sf_site_folder)
                fout_fa.write(s_ins_seq)
####
####
####
    #group by sites and save the reads to one file, and them assemble them
    def asm_clip_disc_reads_one_MEI(self, rcd):
        sf_local_working_folder=rcd[0]
        x_local_asm=XLocalAssembly(sf_local_working_folder,"")
        x_local_asm.run_local_assembly(rcd)

    ####get the assembled insertion sequence
    def _get_asm_sequence(self, s_id, sf_folder):
        l_contigs=[]
        sf_contig=sf_folder+"/contig.fa"
        l_contigs.append(sf_contig)
        for i in [100,80,60,40,20]:
            l_contigs.append(sf_folder+"/contig-%d.fa"%i)

        b_saved=False
        s_ins_seq=">"+s_id+"_ins_seqs\n"
        for sf_tmp_contig in l_contigs:
            if b_saved==True:
                break
            if os.path.isfile(sf_tmp_contig)==False:
                continue
            with pysam.FastxFile(sf_tmp_contig) as fin_contig:
                for entry in fin_contig:
                    s_seq=str(entry.sequence)
                    if b_saved==True:
                        s_ins_seq+="NNNNNNNNNN"
                    s_ins_seq+=s_seq
                    b_saved=True
        s_ins_seq+="\n"
        return s_ins_seq
####