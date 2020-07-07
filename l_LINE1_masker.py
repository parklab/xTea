##03/21/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

####1) Mask the LINE1 insertions;
####2) Call out the transductions.

from l_rep_masker import *
from l_transduction import *

####
class LINE1masker(LRepMasker):
    def __init__(self, swfolder, n_jobs):
        LRepMasker.__init__(self, swfolder, n_jobs)
        self.transdct = LTransduction(swfolder, n_jobs)
        self.i_min_copy_len = 5500
        self.i_chk_win = 150
        self._LINE1 = "LINE1"
        self._min_map_ratio=0.5
        self._polyA_start_pos=6015
        self._min_ins_lenth_for_chk=150


    ####sf_algnmt: alignments of aligning all the contigs to the consensus
    # 1. fully aligned copy (full length or truncated)
    # 2. fully aligned copy with transduction as clip (pay attention to 3' and 5')
    # 3. inversion case (first part aligned, and the other part aligned in reverse-complementary way)
    # 4. with deletion inside the insertion
    # 4. how about insertion with deletion?????????????????????????????????
    # alignment is: align the contig to the consensus
    def parse_LINE1_from_algnmt(self, sf_algnmt, i_max_clip, sf_ref, sf_ref_fl_flank, flk_lenth, sf_out_tmp, sf_out):
        ###each record in format: t_segmt_pri, t_segmt_sup, algnmt, s_transduction_type, i_5mer_chk_len, i_3mer_chk_len
        ########t_segmt_pri in format (start, end, direction[+/-])
        l_transduct = []  # save the algnmts for further checking
        l_succeed = []
        l_succeed_mid_sgmt=[]

        self.set_TSD(sf_ref)
        samfile = pysam.AlignmentFile(sf_algnmt, "rb")  # read in the sam file
        for algnmt in samfile.fetch():  # check each alignment
            if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                continue
            if algnmt.is_unmapped == True:  ##unmapped
                continue
            if algnmt.is_duplicate == True:  ##duplciation
                continue
            # if algnmt.mapping_quality<30:
            #     continue
####
            print "Working on alignment {0}".format(algnmt.query_name) ######################################################
####
            i_max_all_clip = 2 * i_max_clip
            b_fully_map, i_lclip, i_rclip = self.is_fully_map(algnmt, i_max_all_clip)
            i_map_pos = algnmt.reference_start   ####
            if i_map_pos > self._polyA_start_pos:#purely polyA
                print "Purly polyA"
                continue

            b_lclip = False
            b_rclip = False
            if i_lclip > i_max_clip:
                b_lclip = True
            if i_rclip > i_max_clip:
                b_rclip = True
####
            ####check the small clip region agains the flanking region for TSD
            ####some clip part are aligned, which will change b_lclip or b_rclip
            b_lclip, b_rclip=self.check_TSD_for_site(algnmt, i_lclip, i_rclip, b_lclip, b_rclip)
            if (b_fully_map == True) or ((b_lclip is False) and (b_rclip is False)):
                # m_fully_map[query_name] = query_seq
                l_cigar = algnmt.cigar
                i_seq_start = 0
                if l_cigar[0][0] == 4:
                    i_seq_start = l_cigar[0][1]
                i_map_end, n_mapped, i_seq_end = self._get_map_interval(l_cigar, i_map_pos)
                #skip the short aligned but doesn't reach polyA ones
                if n_mapped<self._min_ins_lenth_for_chk and i_map_end<self._polyA_start_pos:
                    continue
                tmp_rcd = (i_map_pos, i_map_end, "+", i_seq_start, i_seq_end)
                if algnmt.is_reverse == True:
                    tmp_rcd = (i_map_pos, i_map_end, "-", i_seq_start, i_seq_end)
                rslt_rcd = self.get_ins_information_from_algnmt_ori_seq(algnmt, self._LINE1, tmp_rcd, None)
                if i_map_end<self.i_min_copy_len:#doesn't reach the end, then pass
                    l_succeed_mid_sgmt.append(rslt_rcd)
                    continue
                l_succeed.append(rslt_rcd)
            else:  # not fully mapped, and one side clipped then check
                ####collect transduction information, and also other fully mapped cases from Supplementary algnmt
                print "With supplementary alignment"
                self.classify_with_supplementary_algnmt(algnmt, b_lclip, b_rclip, self._LINE1,
                                                        l_transduct, l_succeed, self.i_min_copy_len, l_succeed_mid_sgmt)
        samfile.close()
        self.clean_TSD()
####
        ####output the non transudction ones
        sf_rslts_non_transduction = sf_out + "_non_transduction.txt"
        self.save_non_transduct_ins_info_to_file(l_succeed, global_values.LRD_MIN_INTERNAL_DEL,
                                                 sf_rslts_non_transduction)
        sf_rslts_non_transduction_mid_sgmt = sf_out + "_non_transduction_only_middle_sgmt.txt"
        self.save_non_transduct_ins_info_to_file(l_succeed_mid_sgmt, global_values.LRD_MIN_INTERNAL_DEL,
                                                 sf_rslts_non_transduction_mid_sgmt)


        m_candidates, m_with_polya, m_transduct_tsd, m_tbd = self.transdct.call_transduction_from_candidates(
            l_transduct, l_succeed, self.i_min_copy_len, self.i_chk_win, flk_lenth, sf_ref, sf_ref_fl_flank)
########
        #merge the transduction tsd records to existing ones
        for s_tsd_id in m_transduct_tsd:
            if s_tsd_id in self.m_tsd:
                print "May Error: {0} TSD already exists".format(s_tsd_id)
            else:
                self.m_tsd[s_tsd_id]=m_transduct_tsd[s_tsd_id]
        sf_trsdct = sf_out + "_transductions.txt"
        m_not_slct=self.slct_save_transduction_to_file(l_transduct, m_candidates, m_with_polya, self._LINE1, sf_trsdct)

        # to-be-decided cases
        ####save the m_tbd records to file
        ####save the m_not_slct to file
        sf_rslts_tbd = sf_out + "_to_be_decided_cases.txt"
        self.save_tbd_cases_to_file(m_not_slct, sf_rslts_tbd)

        ####merge the outputs
        ####skip the middle size segments
        self.merge_outputs(sf_rslts_non_transduction, sf_trsdct, sf_out)
        ####Keep all
        self.merge_all_outputs(sf_rslts_non_transduction, sf_rslts_non_transduction_mid_sgmt,
                               sf_trsdct, sf_rslts_tbd, sf_out_tmp)

####
    def merge_all_outputs(self, sf_non_transduct1, sf_non_transduct2, sf_transduct, sf_tbd, sf_merged):
        with open(sf_merged, "w") as fout_merged, open(sf_non_transduct1) as fin_non_td, \
                open(sf_non_transduct2) as fin_non_td2, open(sf_transduct) as fin_td, open(sf_tbd) as fin_tbd:
            for line in fin_non_td:
                fout_merged.write(line)
            for line in fin_non_td2:
                fout_merged.write(line)
            for line in fin_td:
                fout_merged.write(line)
            for line in fin_tbd:
                fout_merged.write(line)
####
####