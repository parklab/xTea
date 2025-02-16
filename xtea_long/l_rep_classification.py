##03/21/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com
##added CHM13 support --Feb 15, 2025

####Given assembled insertion, mask the insertion sequence to:
####1) Classify the insertion to: LINE1, Alu, SVA, HERV, Psudogene, Mitochondria, Orphan transduction, other insertion
####2) Call out the transduction events for specific type of insertions
####3) Call associated SVs (TBD)

####Bug: in function gnrt_combined_final_result: duplicate insertions in to-be-decided-cases, and confirmed2
########current solution: filter out the duplicate ones at the merging step (not the optimal solution)
####

import pysam
from x_contig import *
from x_polyA import *
from x_reference import *
from cmd_runner import *
from l_LINE1_masker import *
from l_SVA_masker import *
from l_Alu_masker import *
from l_rep_masker import *
from l_Mitochondria_masker import *
from l_HERV_masker import *
from l_polyA_masker import *
from l_combined_masker import *

####
class LRepClassification():#
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1]!="/":
            self.swfolder+="/"
        self.n_jobs = n_jobs
        self._3mer="3mer"
        self._5mer="5mer"
        self._both_transdct="both"
        self._PolyA="polyA"#decoy polyA sequence, which is used to filter out the polyA only sequences
        self._LINE1="LINE1"
        self._Alu="ALU"
        self._SVA="SVA"
        self._HERV="HERV"
        self._HERV_K="HERV-K"
        self._Mitochondria="Mitochondrion"
        self._ORPHAN="ORPHAN" #orphan transduction
        self._PSUDOGENE = "Psudogene"
        self._COMBINED="Combined"
        self.cmd_runner = CMD_RUNNER()
        self.l_type=[]#save the repeat type
        self.m_sf_flank={}#save the flanking files
        self.l_sf_cns=[]#save the consensus files
        self.l_sf_cns_kmer=[]#this is the consensus files for kmer generating
        self.sva_rmsk=""#

    ####
    def get_sf_cns_kmer_cnt(self):
        l_cns=[]
        for rcd in self.l_sf_cns_kmer:
            l_cns.append(rcd)
        return l_cns
####
    ####i_type: indicates which type of repeats
    def set_rep_configuration(self, i_type, s_rep_folder, b_hg19=False, b_chm13=False):##
        self.l_type=self.set_rep_type(i_type)
        if (s_rep_folder is None) or (len(s_rep_folder)==0):
            print("repeat folder is not set!!!")
            return
        if s_rep_folder[-1]!="/":
            s_rep_folder+="/"

        if b_hg19==False and b_chm13==False:
            for s_type in self.l_type:
                s_tmp_cns=s_rep_folder+"consensus_mask_lrd/"+s_type+".fa"
                self.l_sf_cns.append(s_tmp_cns)#
                if s_type is self._LINE1:
                    s_tmp_flank=s_rep_folder+ "consensus_mask_lrd/hg38_FL_L1_flanks_3k.fa"
                    self.m_sf_flank[s_type]=s_tmp_flank
                elif s_type is self._SVA:
                    s_tmp_flank = s_rep_folder + "consensus_mask_lrd/hg38_FL_SVA_flanks_3k.fa"
                    self.m_sf_flank[s_type] = s_tmp_flank
                    self.sva_rmsk=s_rep_folder + "consensus_mask_lrd/hg38_SVA.out"
                    s_tmp_cns = s_rep_folder + "consensus_mask_lrd/" + s_type + "_ori.fa"

                if s_type is not self._COMBINED:
                    self.l_sf_cns_kmer.append(s_tmp_cns)
        elif b_hg19==False and b_chm13==True:
            for s_type in self.l_type:
                s_tmp_cns=s_rep_folder+"consensus_mask_lrd/"+s_type+".fa"
                self.l_sf_cns.append(s_tmp_cns)#
                if s_type is self._LINE1:
                    s_tmp_flank=s_rep_folder+ "consensus_mask_lrd/chm13_FL_L1_flanks_3k.fa"
                    self.m_sf_flank[s_type]=s_tmp_flank
                elif s_type is self._SVA:
                    s_tmp_flank = s_rep_folder + "consensus_mask_lrd/chm13_FL_SVA_flanks_3k.fa"
                    self.m_sf_flank[s_type] = s_tmp_flank
                    self.sva_rmsk=s_rep_folder + "consensus_mask_lrd/chm13_SVA.out"
                    s_tmp_cns = s_rep_folder + "consensus_mask_lrd/" + s_type + "_ori.fa"

                if s_type is not self._COMBINED:
                    self.l_sf_cns_kmer.append(s_tmp_cns)
        else:
            for s_type in self.l_type:
                s_tmp_cns=s_rep_folder+"consensus_mask_lrd/"+s_type+".fa"
                self.l_sf_cns.append(s_tmp_cns)
                if s_type is self._LINE1:
                    s_tmp_flank=s_rep_folder+ "consensus_mask_lrd/hg19_FL_L1_flanks_3k.fa"
                    self.m_sf_flank[s_type]=s_tmp_flank
                elif s_type is self._SVA:
                    s_tmp_flank = s_rep_folder + "consensus_mask_lrd/hg19_FL_SVA_flanks_3k.fa"
                    self.m_sf_flank[s_type] = s_tmp_flank
                    self.sva_rmsk=s_rep_folder + "consensus_mask_lrd/hg19_SVA.out"
                    s_tmp_cns = s_rep_folder + "consensus_mask_lrd/" + s_type + "_ori.fa"

                if s_type is not self._COMBINED:
                    self.l_sf_cns_kmer.append(s_tmp_cns)##
####
    ####
    def classify_from_ref_algnmt(self, sf_ref, sf_rep_ins, sf_rslt):
        ####align all the sequence to reference genome
        ####get:1) unique fully mapped ones --> for further checking whether candidate orphan transductions
                ####saved in a file named "all_asm_seqs_2_ref_unique_aligned_for_further_checking.sam"
                ####insertion itself has polyA, TSD, and 3' of a full length copy (by repeatmasker)
        ####    2) unmapped ones --> for further checking whether will be aligned to the non-reference transductions
                ####saved in a file named "all_asm_seqs_2_ref_unmapped_for_further_checking.sam"
        ####    3) check whether fall in known psuedogene region --> psuedogene insertion
        ####    4) spliced alignment --> novel psuedogene insertion
        xtea_contig = XTEContig(self.swfolder, self.n_jobs)
        sf_algnmt=self.swfolder + "all_tei_seq_2_ref.bam"
        xtea_contig.align_contigs_2_reference_genome(sf_ref, sf_rep_ins, self.n_jobs, sf_algnmt)
        return


    ####classify the insertion type by aligning to repeat consensus
    def classify_ins_seqs(self, sf_rep_ins, sf_ref, flk_lenth, sf_rslt):
        #first classify by align all the seqs to reference genome
        self.classify_from_ref_algnmt(sf_ref, sf_rep_ins, sf_rslt)
        xtea_contig = XTEContig(self.swfolder, self.n_jobs)
        sf_rep_ins_tmp = sf_rep_ins  # repeat insertion temporary file
############
        for s_type, s_rep_cns in zip(self.l_type, self.l_sf_cns):
            # s_rep_cns = sf_rep_cns_folder + s_type + ".fa"
            # 1. first align the insertion seqs to consensus
            sf_algnmt = self.swfolder + s_type + "_cns.bam"
            # if s_type == self._SVA:  # for SVA, align the consensus segments (special processing) to the assembled contigs
            #     xtea_contig.align_short_contigs_2_cns_minimap2(sf_rep_ins_tmp, s_rep_cns, self.n_jobs, sf_algnmt)
            # else:# otherwise, align the contigs to repeat consensus
############
            print("Working on {0} with contigs {1} and consensus {2}\n".format(s_type, sf_rep_ins_tmp, s_rep_cns))
            xtea_contig.align_short_contigs_2_cns_minimap2(s_rep_cns, sf_rep_ins_tmp, self.n_jobs, sf_algnmt)
            sf_tmp_out=""
            if s_type == self._PolyA:#mask out those polyA only (like to be homopolymer) cases
                sf_tmp_out = sf_rslt + "." + self._PolyA + ".txt"
                polyA_masker = PolyAMasker(self.swfolder, self.n_jobs)
                polyA_masker.parse_polyA_from_algnmt(sf_algnmt, sf_ref, sf_tmp_out)
            elif s_type == self._Alu:
                sf_tmp_out = sf_rslt + "." + self._Alu + ".txt"
                alu_masker=AluMasker(self.swfolder, self.n_jobs)
                alu_masker.parse_Alu_from_algnmt(sf_algnmt, sf_ref, global_values.LRD_PRI_SUP_MAX_CLIP, sf_tmp_out)
            elif s_type == self._HERV:
                sf_tmp_out = sf_rslt + "." +self._HERV + ".txt"
                herv_masker = HERVMasker(self.swfolder, self.n_jobs)
                herv_masker.parse_HERV_from_algnmt(sf_algnmt, sf_ref, global_values.LRD_PRI_SUP_MAX_CLIP, sf_tmp_out)
            elif s_type == self._Mitochondria:
                sf_tmp_out = sf_rslt + "." + self._Mitochondria + ".txt"
                mit_masker=MitochondriaMasker(self.swfolder, self.n_jobs)
                mit_masker.parse_Mit_from_algnmt(sf_algnmt, sf_ref, global_values.LRD_PRI_SUP_MAX_CLIP, sf_tmp_out)
            elif s_type==self._LINE1:
                sf_tmp_out=sf_rslt + "." + self._LINE1 + ".txt"#skip the mid sgmt ones
                sf_tmp_out_tmp=sf_rslt + "." + self._LINE1 + ".txt.tmp" #keep the small middle sgmt ones
                sf_ref_fl_flank=self.m_sf_flank[self._LINE1]
                line1_masker = LINE1masker(self.swfolder, self.n_jobs)
                line1_masker.parse_LINE1_from_algnmt(sf_algnmt, global_values.LRD_PRI_SUP_MAX_CLIP, sf_ref,
                                                     sf_ref_fl_flank, flk_lenth, sf_tmp_out_tmp, sf_tmp_out)
                sf_tmp_out=sf_tmp_out_tmp
################
            elif s_type == self._SVA:
                sf_tmp_out = sf_rslt + "." + self._SVA + ".txt"
                sf_ref_fl_flank = self.m_sf_flank[self._SVA]
                sva_masker = SVAmasker(self.swfolder, self.n_jobs)
                sva_masker.parse_SVA_from_algnmt(sf_algnmt, self.sva_rmsk, global_values.LRD_PRI_SUP_MAX_CLIP,
                                                     sf_ref, sf_ref_fl_flank, flk_lenth, sf_tmp_out)
            elif s_type==self._COMBINED:
                sf_tmp_out = sf_rslt + "." + self._COMBINED + ".txt"
                cmb_masker = CombinedMasker(self.swfolder, self.n_jobs)
                cmb_masker.parse_combined_from_algnmt(
                    sf_algnmt, sf_ref, global_values.LRD_PRI_SUP_MAX_CLIP, sf_tmp_out)
            #################
            # 2. The previous step will generate a new "sf_rep_ins" file, and pass to the other rep type
            sf_new_tmp = sf_rep_ins_tmp + ".after_" + s_type + "_round.fa"
            self.get_unmasked_seqs(sf_rep_ins_tmp, sf_tmp_out, sf_new_tmp)
            sf_rep_ins_tmp = sf_new_tmp
        ####
        #rename the final un-classified insertion-seqs (for pseudogene insertion purpose)
        os.rename(sf_rep_ins_tmp, sf_rep_ins+global_values.LRD_PSEUDOGENE_INPUT_SUF)
        ####generate the merged results
        self.gnrt_combined_final_result(sf_rslt)
####
####
    ####
    def gnrt_combined_final_result(self, sf_rslt):
        l_alu = []
        l_line1 = []
        l_sva = []
        l_herv = []
        l_mit = []
        for s_type in self.l_type:
            if s_type is self._Alu:
                sf_alu = sf_rslt + "." + self._Alu + ".txt"
                sf_alu_tbd= sf_alu+"_to_be_decided_cases.txt"
                l_alu.append(sf_alu)
                l_alu.append(sf_alu_tbd)
            elif s_type is self._LINE1:
                sf_line1 = sf_rslt + "." + self._LINE1 + ".txt"
                sf_line1_tbd=sf_line1+"_to_be_decided_cases.txt"
                sf_line1_tbd2=sf_rslt + "." + self._COMBINED + ".txt"+"_" + self._LINE1 + "_to_be_further_confirmed2.txt"
                l_line1.append(sf_line1)
                l_line1.append(sf_line1_tbd)
                l_line1.append(sf_line1_tbd2)
            elif s_type is self._SVA:
                sf_sva=sf_rslt + "." + self._SVA + ".txt"
                sf_sva_tbd = sf_sva + "_to_be_decided_cases.txt"
                sf_sva_tbd2 = sf_rslt + "." + self._COMBINED + ".txt" + "_" + self._SVA + "_to_be_further_confirmed2.txt"
                #classified_results.txt.SVA.txt_non_transduction_within_rep_copy.txt
                sf_in_ref_copy=sf_sva+"_non_transduction_within_rep_copy.txt"
                sf_sva_tbd_no_ref_copy = sf_sva + "_to_be_decided_cases_no_ref_copy.txt"
####Hard code here!!!!!!!!!!!!!
                b_with_chr=True
                self.sprt_SVA_fall_in_rep_copy(self.sva_rmsk, b_with_chr, sf_sva_tbd2, sf_in_ref_copy, sf_sva_tbd_no_ref_copy)
                l_sva.append(sf_sva)
                l_sva.append(sf_sva_tbd)
                l_sva.append(sf_sva_tbd_no_ref_copy)
            elif s_type is self._HERV:
                sf_herv = sf_rslt + "." + self._HERV + ".txt"
                sf_herv_tbd2 = sf_rslt + "." + self._COMBINED + ".txt" + "_" + self._HERV_K + "_to_be_further_confirmed2.txt"
                l_herv.append(sf_herv)
                l_herv.append(sf_herv_tbd2)
####
        ####merge LINE1
        sf_merged_L1=sf_rslt+".merged_" + self._LINE1 + ".txt"
        self.merge_rslt_without_dup_for_list(l_line1, sf_merged_L1)
        ####merge Alu
        sf_merged_Alu = sf_rslt + ".merged_" + self._Alu + ".txt"
        self.merge_rslt_for_list(l_alu, sf_merged_Alu)
        ####merge SVA
        sf_merged_SVA = sf_rslt + ".merged_" + self._SVA + ".txt"
        self.merge_rslt_without_dup_for_list(l_sva, sf_merged_SVA)
        ####merge HERV
        sf_merged_HERV = sf_rslt + ".merged_" + self._HERV + ".txt"
        self.merge_rslt_for_list(l_herv, sf_merged_HERV)

    #get merge results list
    def get_merged_results(self, sf_prefix, i_type):
        l_reps=self.set_rep_type(i_type)
        l_rslts=[]
        for s_type in l_reps[:-1]:
            l_rslts.append((s_type, sf_prefix+".merged_" + s_type + ".txt"))
        return l_rslts

    ####
    def merge_rslt_for_list(self, l_rslt_list, sf_out):
        with open(sf_out, "w") as fout_merged:
            for sf_rslt in l_rslt_list:
                with open(sf_rslt) as fin_rslt:
                    for line in fin_rslt:
                        fout_merged.write(line)

    def merge_rslt_without_dup_for_list(self, l_rslt_list, sf_out):
        m_exist={}
        with open(sf_out, "w") as fout_merged:
            for sf_rslt in l_rslt_list:
                with open(sf_rslt) as fin_rslt:
                    for line in fin_rslt:
                        fields=line.split()
                        chrm=fields[0]
                        pos=int(fields[1])
                        if chrm not in m_exist:
                            m_exist[chrm]={}
                        i_start=pos-global_values.NEARBY_REGION
                        i_end=pos+global_values.NEARBY_REGION
                        for i_tmp_pos in range(i_start, i_end):
                            if i_tmp_pos in m_exist[chrm]:
                                continue
                        m_exist[chrm][pos]=1
                        fout_merged.write(line)
####
    def sprt_SVA_fall_in_rep_copy(self, sf_rmsk, b_with_chr, sf_un_sprted, sf_in_ref_copy, sf_ins):
        sva_masker = SVAmasker(self.swfolder, self.n_jobs)
        xannotation=sva_masker.construct_interval_tree(sf_rmsk, b_with_chr)
        with open(sf_un_sprted) as fin_ori, open(sf_in_ref_copy,"a") as fappnd, open(sf_ins, "w") as fout_ins:
            for line in fin_ori:
                fields=line.split()
                ins_chrm=fields[0]
                ins_pos=int(fields[1])
                b_in_ref_sva, rep_start_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
                if b_in_ref_sva is True:
                    fappnd.write(line)
                else:
                    fout_ins.write(line)

####
    def get_unmasked_seqs(self, sf_ori, sf_slcted, sf_new):
        m_slcted={}
        with open(sf_slcted) as fin_slcted:
            for line in fin_slcted:
                fields=line.split()
                ins_chrm=fields[0]
                ins_pos=fields[1]
                s_id="{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                m_slcted[s_id]=1

        #with open(sf_new, "w") as fout_new, open(sf_ori) as fin_ori:
        with pysam.FastxFile(sf_ori) as fin_ori, open(sf_new, "w") as fout_new:
            for entry in fin_ori:
                fields=entry.name.split(global_values.SEPERATOR)
                s_id=global_values.SEPERATOR.join(fields[:-1])
                if s_id not in m_slcted:
                    fout_new.write(">"+entry.name+"\n")
                    fout_new.write(entry.sequence+"\n")

####
    def set_rep_type(self, i_type):
        l_rep = []
        l_rep.append(self._PolyA)
        if i_type & 1 != 0:
            l_rep.append(self._LINE1)
        if i_type & 2 != 0:
            l_rep.append(self._Alu)
        if i_type & 4 != 0:
            l_rep.append(self._HERV) ####Note the order for HERV and SVA
        if i_type & 8 != 0:
            l_rep.append(self._SVA)
        if i_type & 16 != 0:
            l_rep.append(self._Mitochondria)

        if i_type & 32 != 0:
            l_rep.append(self._ORPHAN)
        if i_type & 64 != 0:
            l_rep.append(self._PSUDOGENE)
        l_rep.append(self._COMBINED)
        return l_rep
#
####
    def get_mask_complex_sv_rep_type(self):
        return 1+2+4+8+16

    def get_output_list_by_type(self, sf_prefix, i_type):
        l_types=self.set_rep_type(i_type)
        l_rslts=[]
        for s_type in l_types:
            sf_merged = sf_prefix + ".merged_" + s_type + ".txt"
            l_rslts.append(sf_merged)
        return l_rslts
####