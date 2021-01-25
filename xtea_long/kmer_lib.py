##12/15/2020
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu or chong.simon.chu@gmail.com

####To do list:
'''
0. separate by subfamily (done)
1. whether generate the TE lib with m-distance kmers???? (done, added separate kmer one)
2. how to set cutoff???
3. check the frequency of kmers???
'''
import os
from kmer_operator import *####
try:
    import cPickle as pickle
except ImportError:  # Python 3.x
    import pickle

class KmerLib():
    def __init__(self):
        self.k_len=15
        self.k_len_polyA=11 #kmer length for polyA
        self.l_kmer_lib=[]#this is the list of m_kmer_libs
        self.l_kmer_lib_rc=[]#this is for the list of m_kmer_libs (reverse complementary)####

        #this polyA kmer is used to remove the polyAs from the kmer lib
        b_rc=True#
        self.m_polyA_kmers=self.init_polyA_T_kmers("A", self.k_len, b_rc)

        b_rc=False
        #this polyA kmer is of a smaller k, and seperate to two: RC vs non-RC
        self.m_polyA_kmers2 = self.init_polyA_T_kmers("A", self.k_len_polyA, b_rc)
        self.m_polyA_kmers2_rc = self.init_polyA_T_kmers("T", self.k_len_polyA, b_rc)####

    ####
    def init_polyA_T_kmers(self, c_A_T, i_k, b_rc):
        kmerop = KmerOperator(i_k)
        s_polyA = ""
        for i in range(i_k):
            s_polyA += c_A_T
        m_polyA_T_kmers = kmerop.gnrt_one_substitution_seqs(s_polyA, b_rc)
        return m_polyA_T_kmers

    def set_k(self, i_k):
        self.k_len=i_k

    def _gnrt_kmer_from_seq(self, s_seq, m_local_kmer=None):##
        m_local={}
        i_n_kmer=len(s_seq)-self.k_len+1
        if i_n_kmer<=0:##
            return m_local
        else:#
            for i in range(i_n_kmer):#
                s_kmer=s_seq[i:i+self.k_len]
                if s_kmer not in m_local:
                    m_local[s_kmer]=1
                else:
                    m_local[s_kmer] += 1

                if m_local_kmer is None:
                    continue
                if s_kmer in self.m_polyA_kmers:#skip polyA
                    continue
                if s_kmer not in m_local_kmer:
                    m_local_kmer[s_kmer]=1
                else:
                    m_local_kmer[s_kmer] += 1
        return m_local

    ####
    def _gnrt_kmer_from_seq2(self, s_seq, i_k):##
        m_local={}
        i_n_kmer=len(s_seq)-i_k+1
        if i_n_kmer<=0:##
            return m_local
        else:#
            for i in range(i_n_kmer):#
                s_kmer=s_seq[i:i+i_k]
                if s_kmer not in m_local:
                    m_local[s_kmer]=1
                else:
                    m_local[s_kmer] += 1
        return m_local

    ####
    def add_to_kmer_lib(self, sf_fa, m_local_kmer, b_rc=False):#
        kmerop=KmerOperator(self.k_len)
        with open(sf_fa) as fin_fa:
            for line in fin_fa:##
                if (len(line)>0 and line[0]==">"):
                    continue
                #1.generate kmers for inital seqs
                s_seq_ori=line.rstrip()
                s_seq=kmerop.cvt_to_upper_case_seq(s_seq_ori)#all are converted to uppercase format
                if b_rc==True:
                    s_seq_rc = kmerop.gnrt_reverse_complementary(s_seq)
                    self._gnrt_kmer_from_seq(s_seq_rc, m_local_kmer)
                else:
                    self._gnrt_kmer_from_seq(s_seq, m_local_kmer)

                #2. generate one substitution seqs and add to the kmer lib
                m_seq_subs=kmerop.gnrt_one_substitution_seqs(s_seq, False)#
                for s_new_seq in m_seq_subs:
                    if b_rc==False:
                        self._gnrt_kmer_from_seq(s_new_seq, m_local_kmer)
                    else:
                        s_new_seq_rc=kmerop.gnrt_reverse_complementary(s_new_seq)
                        self._gnrt_kmer_from_seq(s_new_seq_rc, m_local_kmer)

                #3. generate one insertion seqs and add to the kmer lib
                m_seq_ins=kmerop.gnrt_one_insertion_distance_seqs(s_seq)#
                for s_new_seq in m_seq_ins:
                    if b_rc == False:
                        self._gnrt_kmer_from_seq(s_new_seq, m_local_kmer)
                    else:
                        s_new_seq_rc = kmerop.gnrt_reverse_complementary(s_new_seq)
                        self._gnrt_kmer_from_seq(s_new_seq_rc, m_local_kmer)

                #4. generate one deletion seqs and add to the kmer lib
                m_seq_del=kmerop.gnrt_one_deletion_distance_seqs(s_seq)#
                for s_new_seq in m_seq_del:
                    if b_rc == False:
                        self._gnrt_kmer_from_seq(s_new_seq, m_local_kmer)
                    else:
                        s_new_seq_rc = kmerop.gnrt_reverse_complementary(s_new_seq)##
                        self._gnrt_kmer_from_seq(s_new_seq_rc, m_local_kmer)

    ####Given a list of files of different repeat subfamilies, construct the TE kmer library
    #First check whether the pickle file exist or not, if not, then create; otherwise, directly load from file
    def construct_kmer_lib2(self, l_sf_lib):##
        for sf_fa in l_sf_lib:
            m_local_kmer = {}  #
            sf_pickle=sf_fa+".p"
            if os.path.isfile(sf_pickle)==True:
                m_local_kmer=self._load_from_pickle(sf_pickle)
            else:
                self.add_to_kmer_lib(sf_fa, m_local_kmer)
                self._save_to_pickle(m_local_kmer, sf_pickle)
            self.l_kmer_lib.append(m_local_kmer)

            #this is for reverse complementary kmers
            b_rc=True
            m_local_kmer_rc={}
            sf_pickle_rc = sf_fa + ".rc.p"
            if os.path.isfile(sf_pickle_rc)==True:
                m_local_kmer_rc=self._load_from_pickle(sf_pickle_rc)
            else:
                self.add_to_kmer_lib(sf_fa, m_local_kmer_rc, b_rc)
                self._save_to_pickle(m_local_kmer_rc, sf_pickle_rc)
            self.l_kmer_lib_rc.append(m_local_kmer_rc)

    def is_seq_contain_TE_kmer_rc_non_rc(self, s_seq, n_non_polyA_cutoff, n_polyA_cutoff, b_rc=False):
        if b_rc==False:
            b_contained=self.is_seq_contain_TE_kmer(s_seq, n_non_polyA_cutoff, n_polyA_cutoff, self.l_kmer_lib,
                                                    self.m_polyA_kmers2)
            return b_contained
        else:
            b_contained_rc = self.is_seq_contain_TE_kmer(s_seq, n_non_polyA_cutoff, n_polyA_cutoff, self.l_kmer_lib_rc,
                                                      self.m_polyA_kmers2_rc)
            return b_contained_rc

    ####
    def is_seq_contain_TE_kmer(self, s_seq, n_non_polyA_cutoff, n_polyA_cutoff, l_kmer_lib, m_polyA_kmers):
        m_local=self._gnrt_kmer_from_seq(s_seq)
        m_local_polyA=self._gnrt_kmer_from_seq2(s_seq, self.k_len_polyA)
        b_contained=False
        for m_type_kmer in l_kmer_lib:####
            n_cnt_non_polyA=0
            n_cnt_polyA=0
            b_non_polyA_part = False
            b_polyA = False
            for s_kmer in m_local:
                if s_kmer in m_type_kmer:
                    n_cnt_non_polyA+=1
            for s_kmer_polyA in m_local_polyA:
                if s_kmer_polyA in m_polyA_kmers:
                    n_cnt_polyA+=1
            if n_cnt_non_polyA>=n_non_polyA_cutoff:
                b_non_polyA_part=True
            if n_cnt_polyA>=n_polyA_cutoff:
                b_polyA=True
            #print(n_cnt_non_polyA, n_cnt_polyA)##################################################################
            b_contained=(b_non_polyA_part and b_polyA)
            if b_contained==True:
                return True
        return b_contained


    def _save_to_pickle(self, m_kmer, sf_out_pickle):
        with open(sf_out_pickle, 'wb') as fp:
            pickle.dump(m_kmer, fp, protocol=pickle.HIGHEST_PROTOCOL)

    def _load_from_pickle(self, sf_pickle):
        with open(sf_pickle, 'rb') as fp:
            m_dict = pickle.load(fp)
            return m_dict
####
####
# sf_test="/n/data1/hms/dbmi/park/simon_chu/projects/xTEA_long_reads/results2/HG002/HG002_ccs/tmp/l_asm_tmp/3500/chr2~116215660_3500.fa"
# ####
# sf_lib="/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/rep_lib_annotation/rep_lib_annotation/"
# ##
# l_cns=[]####
# l_cns.append(sf_lib+"consensus_mask_lrd/LINE1.fa")#LINE1
# l_cns.append(sf_lib+"consensus_mask_lrd/ALU.fa")#Alu
# l_cns.append(sf_lib+"consensus_mask_lrd/SVA.fa")#SVA
#
# klib=KmerLib()
# klib.construct_kmer_lib2(l_cns)
#
# n_hit_cnt=0
# n_min_non_polyA_kmer_hit=1
# n_min_polyA_hit=1
# n_cutoff=1
# with open(sf_test) as fin_fa:
#     for line in fin_fa:
#         if len(line) > 0 and line[0] == ">":  #
#             continue
#         s_fa = line.rstrip()
#         b_seq_hit = klib.is_seq_contain_TE_kmer(s_fa, n_min_non_polyA_kmer_hit, n_min_polyA_hit)
#         if b_seq_hit == True:
#             n_hit_cnt += 1
#             if n_hit_cnt >= n_cutoff:  ######
#                 b_hit = True
#                 # break
#     print(n_hit_cnt)

####