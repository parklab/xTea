from l_rep_masker import *
from l_SVA_masker import *

####This is one step after the single type masker step.
####Realign the uncharacterized ones to the consensus again, and pick those:
####1. most part is aligned to one type
####2. has polyA or TSD
class CombinedMasker(LRepMasker):
    def __init__(self, swfolder, n_jobs):
        LRepMasker.__init__(self, swfolder, n_jobs)
        self._alu="ALU" ###these are the ids saved in the "combined"
        self._line1="LINE1"
        self._herv="HERV-K"
        self._mit="NC_012920.1"
        self._sva="SVA"
        self.m_type={}
        self.m_type[self._alu]="Alu"
        self.m_type[self._line1] = "LINE1"
        self.m_type[self._herv] = "HERV-K"
        self.m_type[self._mit] = "Mitochondria"
        #others is SVA copies
        self.min_map_ratio=0.75
        self._full_l1=5000
        self._full_sva=350#suppose at least a
        self._l1_3mer_tail=5500
        self._polyA_start_pos = 6015
        self._sva_flank_ovlp_slack=25
        self.polyA_dom_ratio=0.6
        self.polyA_dom_ratio2=0.4
        self.small_sgmt_len=100

####
    #Here i_max_clip is maximum allowed TSD
    def parse_combined_from_algnmt(self, sf_algnmt, sf_ref, i_max_clip, sf_out):
        s_open_fmt = "rb"
        m_succeed={}
        m_succeed[self._alu]=[]
        m_succeed[self._line1] = []
        m_succeed[self._herv] = []
        m_succeed[self._mit] = []
        m_succeed[self._sva] = []

        self.set_TSD(sf_ref)
        sva_masker = SVAmasker(self.swfolder, self.n_jobs)
        i_sva_flank_len = sva_masker.get_flank_length()-500 #as sometimes the benchmark is not good enough
        samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
        m_refs = sva_masker.get_cns_copy_length(samfile)
        xpolyA = PolyA()
        for algnmt in samfile.fetch():  # check each alignment
            if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                continue
            if algnmt.is_unmapped:  ##unmapped
                continue
            if algnmt.is_duplicate:  ##duplication
                continue

            l_cigar = algnmt.cigar
            i_map_pos = algnmt.reference_start
            i_seq_start = 0
            if l_cigar[0][0] == 4:
                i_seq_start = l_cigar[0][1]

            i_map_end, n_mapped, i_seq_end = self._get_map_interval(l_cigar, i_map_pos)

            s_tmp_seq = algnmt.query_sequence
            if xpolyA.is_dominant_polyA(s_tmp_seq, self.polyA_dom_ratio):#if >60% of polyA, then skip
                continue
            elif n_mapped<self.small_sgmt_len and xpolyA.is_dominant_polyA(s_tmp_seq, self.polyA_dom_ratio2):
                #if small segment, and >40% polyA, then skip
                print(algnmt.query_name, "Small segments with polyA dominant!")
                continue

            i_max_all_clip = 2 * i_max_clip
            b_fully_map, i_lclip, i_rclip = self.is_fully_map(algnmt, i_max_all_clip)
            b_lclip = True  ###not used
            b_rclip = True  ###not used
            ####check the small clip region agains the flanking region for TSD
            b_lclip, b_rclip = self.check_TSD_for_site(algnmt, i_lclip, i_rclip, b_lclip, b_rclip)
            s_mapped_ref = algnmt.reference_name
            #print "combined {0} {1} {2} {3} {4}".format(i_map_end, i_map_pos, algnmt.query_name, s_mapped_ref, self._l1_3mer_tail)
            if s_mapped_ref not in self.m_type:
                s_mapped_ref = self._sva

            ####
            b_l1_pass=False
            if (s_mapped_ref == self._line1) and (n_mapped>self._full_l1):
                b_l1_pass=True
            b_sva_pass=False
            if (s_mapped_ref == self._sva) and (n_mapped>self._full_sva):
                b_sva_pass=True

            i_seq_len=len(algnmt.query_sequence)
            if self._is_qualified_algnmt(i_seq_len, n_mapped)==False and b_l1_pass==False and b_sva_pass==False:
                print(algnmt.query_name, "not qualified alignmnt!")
                continue

            # for SVA, if within the flank region, then filter out
            if s_mapped_ref == self._sva:
                b_in_flank, bbtmp, n_lclip, n_rclip=sva_masker.is_fully_map_with_flank(
                    m_refs, algnmt, i_max_clip, i_sva_flank_len, self._sva_flank_ovlp_slack)
                if b_in_flank==True:
                    print(algnmt.query_name, "Fall in flanking region!")
                    continue
                #here i_map_end is the last mappable position on the "reference"
                if sva_masker.hit_front_end_only(i_map_end)==True:#if only hit front end, then skip
                    print(algnmt.query_name, "Only hit front end!")
                    continue
####
            #print "combined {0} {1} {2} {3} {4}".format(i_map_end, i_map_pos, algnmt.query_name, s_mapped_ref, self._line1, self._l1_3mer_tail)
            if s_mapped_ref == self._line1:####for line1, skip the middle ones
                if int(i_map_end) < self._l1_3mer_tail:
                    #print "combined2 {0} {1} {2}".format(i_map_end, i_map_pos, algnmt.query_name)
                    continue
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
####Need further checking!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            i_seq_start=0
            i_seq_end=len(s_tmp_seq)
            tmp_rcd = (i_map_pos, i_map_end, "+", i_seq_start, i_seq_end)#
            if algnmt.is_reverse == True:
                tmp_rcd = (i_map_pos, i_map_end, "-", i_seq_start, i_seq_end)

            rslt_rcd = self.get_ins_information_from_algnmt_ori_seq(algnmt, s_mapped_ref, tmp_rcd, None)
            m_succeed[s_mapped_ref].append(rslt_rcd)

        samfile.close()
        self.clean_TSD()
        self.save_non_transduct_ins_info_to_file_multi_types(m_succeed, global_values.LRD_MIN_INTERNAL_DEL, sf_out)
####
####
    ## most of part is aligned,
    # or it convers the whole repeat consensus
    def _is_qualified_algnmt(self, i_len, i_mapped):
        b_qualified=True
        if i_len<=0:
            return False
        if float(i_mapped)/float(i_len) < self.min_map_ratio:
            return False
        return b_qualified

####
    def save_non_transduct_ins_info_to_file_multi_types(self, m_passed, i_del_slack, sf_out, i_ofst=0):
        with open(sf_out, "w") as fout_merged:
            for s_type in m_passed:
                sf_out_type=sf_out+ "_" + s_type+"_to_be_further_confirmed2.txt"
                with open(sf_out_type, "w") as fout_ins:
                    lis = LInternalStructure()
                    xpolyA=PolyA()
                    l_succeed=m_passed[s_type]
                    for rcd in l_succeed:
                        s_rg1=rcd[4]#region1
                        s_rg2=rcd[5]#region2

                        (ins_chrm, ins_pos, rep_type, structure_type, s_region1, s_region2, s_TSD, transduct_src,
                         s_seq, transduct_seq)=rcd
                        s_structure=lis.get_internal_structure(s_rg1, s_rg2, i_del_slack)
                        s_with_polyA = self.chk_polyA(xpolyA, s_seq)

                        # if s_with_polyA==False:
                        #     s_with_polyA = self.chk_polyA_large_region(xpolyA, s_seq)

                        # b_rpolyA=xpolyA.is_consecutive_polyA_T(seq[:10])
                        # b_lpolyA=xpolyA.is_consecutive_polyA_T(seq[-10:])
                        # if b_lpolyA or b_rpolyA:
                        #     s_with_polyA="with_polyA"
                        s_id = "{0}~{1}".format(ins_chrm, ins_pos)
                        if s_id in self.m_tsd:
                            s_TSD=self.m_tsd[s_id]

                        if len(s_seq)<global_values.LRD_MIN_INS_LTH:
                            print(ins_chrm, ins_pos, "is smaller than minimal length!")
                            continue

                        if rep_type is self._sva:#here require polyA !!!!
                            if s_with_polyA is self._s_no_polyA:
                                print(ins_chrm, ins_pos, "don't have polyA detected!")
                                continue
                        # if self.is_no_polA_no_TSD(s_with_polyA, s_TSD):
                        #     continue

                        if i_ofst != 0:
                            s_region1 = self.adjust_pos(s_region1, i_ofst)
                            s_region2 = self.adjust_pos(s_region2, i_ofst)
                        sinfo = ins_chrm + "\t" + str(
                            ins_pos) + "\t" + rep_type + "\t" + s_structure + "\t" + s_region1 + "\t" + s_region2 + "\t" + \
                                s_TSD + "\t" + transduct_src + "\t" + transduct_src + "\t" + s_seq + "\t" + transduct_seq + \
                                "\t" + transduct_seq + "\t" + s_with_polyA + "\n"
                        fout_ins.write(sinfo)
                        fout_merged.write(sinfo)
####