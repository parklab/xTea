import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool

'''

'''

####

class XRefCopy():

    ####Call out the snp from repeat copy alignments: From the MD field
    # The bam file is process by "samtools calmd -e " first.
    def call_snp_indel_from_rep_copy_algnmt(self, sf_bam, sf_ref):

        bamfile = pysam.AlignmentFile(sf_bam, "r", reference_filename=sf_ref)

        m_all_mut = {}  # this is position based mutation
        m_copy_mut = {}  # this save the mutations of each copy
        for algnmt in bamfile.fetch():
            if algnmt.has_tag("MD") == False:
                continue
            qname = algnmt.query_name

            map_pos = algnmt.reference_start
            q_seq = algnmt.query_sequence
            s_md = algnmt.get_tag("MD")
            l_cigar = algnmt.cigar

            if (l_cigar[0][0] == 4 and l_cigar[0][1] > 200) or (
                    l_cigar[-1][0] == 4 and l_cigar[-1][1] > 200):  # skip the soft-clip cases
                continue

            m_mut = self.call_snp_indel(l_cigar, map_pos, s_md, q_seq)
            m_copy_mut[qname] = m_mut

        bamfile.close()

        return m_copy_mut
####
    ####
    # Input:
    # 1. l_cigar is a list of cigar strs. 0-M, 1-I, 2-D, 4-S, 5-H, e.g. [(0, 72), (1, 1), (0, 602), (2, 4), (0, 2301), (2, 1), (0, 3053)],
    # 2. s_md is the MD field.
    # In MD field, ^ indicates the start of deletion.
    # e.g. 1G58G98C70T22T97G25G104G67C123^CAGA29T63G80G145G198G38C87T139T69C125T154T146T6A405C188C45G105C262^C350T100C10T510C152A17C82C43G14G0T288T154T5A2T176A107A86C66T230C116C26T2G173A52C65A99A102

    # 3. q_seq is the query seq #this is processed by "samtools calcmd -e". Insertion is showed in seq, but deletion is not.
    # 4. r_seq is the reference seq #this is the original "ref consensus" seq
    # all the postion will be re-calc based on the position on the original ref
    def call_snp_indel(self, l_cigar, ref_start, s_md, q_seq):
        m_mut = {}
        m_mut["M"] = {}  # mismatch
        m_mut["I"] = {}  # insertion
        m_mut["D"] = {}  # deletion

        # parse out deletions and mismatches by order
        l_del, l_mismatch_query = self.call_mismatch_del_from_md(s_md)

        idx_mismatch = 0  # index for l_mismatch_query which save all the mismatches
        idx_del = 0  # index for l_del which save all the mismatches
        ref_pos = ref_start  # global coordinate on the ref, start from 0

        q_start = 0
        q_end = 0
        for (itype, ilen) in l_cigar:
            if itype == 0:  # match or mismatch
                # get the query segment
                q_start = q_end
                q_end += ilen
                s_tmp_q = q_seq[q_start:q_end]

                #             print s_tmp_q
                #             print idx_mismatch, "before"
                #             print ref_pos

                l_tmp_mismatch, idx_mismatch = self.parse_mismatch(s_tmp_q, l_mismatch_query, idx_mismatch, ref_pos)
                #             print idx_mismatch, "after"
                #             print l_tmp_mismatch

                # save the mismatch to dict
                for (ipos, c_r, c_q) in l_tmp_mismatch:
                    m_mut["M"][ipos] = (c_r, c_q)
                ref_pos += ilen
            elif itype == 1:  # insertion
                q_start = q_end
                q_end += ilen
                s_tmp_q = q_seq[q_start:q_end]
                m_mut["I"][ref_pos] = ("", s_tmp_q)
            elif itype == 2:  # deletion
                m_mut["D"][ref_pos] = (l_del[idx_del], "")
                idx_del += 1
                ref_pos += ilen
            elif itype == 4 or itype == 5:  # soft or hard clip
                q_start = q_end
                q_end += ilen
                ref_pos += ilen
            else:
                print "Cigar contain unprocessed type!", l_cigar

        return m_mut

    ####parse out the mismatch for a given region
    def parse_mismatch(self, q_seq, l_mismatch_query, idx_mismatch, ref_pos):
        l_mismatch = []
        istart = ref_pos
        for idx, val in enumerate(q_seq):
            if val == "=":
                istart += 1
                continue
            else:
                istart += 1
                l_mismatch.append((istart, l_mismatch_query[idx_mismatch], q_seq[idx]))  # save the mismatch pair
                idx_mismatch += 1
        return l_mismatch, idx_mismatch

    #
    ####
    def call_mismatch_del_from_md(self, s_md):
        l_del = []
        l_mismatch = []
        b_on = False
        s_del = ""
        for s in s_md:
            if s == '^':
                b_on = True
                continue
            if b_on == True:
                if s >= 'A' and s <= 'Z':
                    s_del += s
                else:
                    l_del.append(s_del)
                    b_on = False
                    s_del = ""
            else:
                if s >= 'A' and s <= 'Z':
                    l_mismatch.append(s)

        # for the last one if possible
        if s_del != "":
            l_del.append(s_del)
        return l_del, l_mismatch

    ####with the mutations of each copy, generate the singleton dict for each copy
    # Input: m_copy_mut in format: {sid:{'M':{ipos:{r , q})}}
    def call_singleton_per_copy(self, m_copy_mut):
        m_singleton = {}

        m_pos_mut_M = {}
        m_pos_mut_I = {}
        m_pos_mut_D = {}

        for sid in m_copy_mut:
            for ipos in m_copy_mut[sid]['M']:  # for mismatch
                if ipos not in m_pos_mut_M:
                    m_pos_mut_M[ipos] = {}
                s_chg = "{0}_{1}".format(m_copy_mut[sid]['M'][ipos][0], m_copy_mut[sid]['M'][ipos][1])
                if s_chg not in m_pos_mut_M[ipos]:
                    m_pos_mut_M[ipos][s_chg] = []
                m_pos_mut_M[ipos][s_chg].append(sid)

            for ipos in m_copy_mut[sid]['I']:  # for insertion
                if ipos not in m_pos_mut_I:
                    m_pos_mut_I[ipos] = {}
                s_chg = "{0}_{1}".format(m_copy_mut[sid]['I'][ipos][0], m_copy_mut[sid]['I'][ipos][1])
                if s_chg not in m_pos_mut_I[ipos]:
                    m_pos_mut_I[ipos][s_chg] = []
                m_pos_mut_I[ipos][s_chg].append(sid)

            for ipos in m_copy_mut[sid]['D']:  # for deletion
                if ipos not in m_pos_mut_D:
                    m_pos_mut_D[ipos] = {}
                s_chg = "{0}_{1}".format(m_copy_mut[sid]['D'][ipos][0], m_copy_mut[sid]['D'][ipos][1])
                if s_chg not in m_pos_mut_D[ipos]:
                    m_pos_mut_D[ipos][s_chg] = []
                m_pos_mut_D[ipos][s_chg].append(sid)

        # find out singleton mismatches
        for ipos in m_pos_mut_M:
            for s_chg in m_pos_mut_M[ipos]:
                if len(m_pos_mut_M[ipos][s_chg]) == 1:  # singleton
                    sid = m_pos_mut_M[ipos][s_chg][0]
                    if sid not in m_singleton:
                        m_singleton[sid] = []
                    s_info = "{0}_{1}".format(ipos, s_chg)
                    m_singleton[sid].append(s_info)

        for ipos in m_pos_mut_I:  ####insertion
            for s_chg in m_pos_mut_I[ipos]:
                if len(m_pos_mut_I[ipos][s_chg]) == 1:  # singleton
                    sid = m_pos_mut_I[ipos][s_chg][0]
                    if sid not in m_singleton:
                        m_singleton[sid] = []
                    s_info = "{0}_{1}".format(ipos, s_chg)
                    m_singleton[sid].append(s_info)

        for ipos in m_pos_mut_D:  # deletion
            for s_chg in m_pos_mut_D[ipos]:
                if len(m_pos_mut_D[ipos][s_chg]) == 1:  # singleton
                    sid = m_pos_mut_D[ipos][s_chg][0]
                    if sid not in m_singleton:
                        m_singleton[sid] = []
                    s_info = "{0}_{1}".format(ipos, s_chg)
                    m_singleton[sid].append(s_info)
        return m_singleton