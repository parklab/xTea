##03/21/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

#To-do-list:
#1) TSD is not decided (function get_ins_information_from_algnmt) [solved]
#2)TSD may help to refine the exact position, need more debug on this
#3) some module need to be moved to alignment related class, some duplicate code is used in l_transduction.py
####

from x_contig import *
from x_polyA import *
from x_reference import *
from cmd_runner import *
from global_values import *
from l_output_fmt_parser import *
from l_TSD import *

####
class LRepMasker():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1]!="/":
            self.swfolder+="/"
        self.n_jobs = n_jobs
        self.lbinfo=LRegionNameBasicInfo()
        self._3mer=self.lbinfo.get_s_3mer()
        self._5mer=self.lbinfo.get_s_5mer()
        self._both_transdct=self.lbinfo.get_s_both_side()
        self._small_polyA_max_clip_len=8 #if hit the end side only with larger clip, then skip
        self.ltsd=None
        self.m_tsd={}
        self._s_no_polyA="no_polyA_found"
        self._s_with_polyA="with_polyA"
        self._f_min_map_ratio=0.66 #for transduction, minimum ratio for all mapped bases (including transduction part)

    ####
    def set_TSD(self,  sf_ref):
        self.ltsd=LTSD(sf_ref)
        self.ltsd.open_ref()

    def clean_TSD(self):
        self.ltsd.close_ref()

    def check_TSD_for_site(self, algnmt, i_lclip, i_rclip, b_lclip, b_rclip):
        # check TSD matching
        b_lclip1=b_lclip
        b_rclip1=b_rclip
        ins_chrm, ins_pos = self.get_ins_pos_from_query_name(algnmt)
        s_lregion, s_rregion = self.ltsd.get_site_short_flanks(ins_chrm, int(ins_pos))#

        b_tsd = False
        if i_lclip < global_values.LRD_MAX_TSD_LEN and i_lclip>=global_values.LRD_MIN_TSD_LEN:
            s_clip = algnmt.query_sequence[:i_lclip]
            s_lregion1=s_lregion[-2*i_lclip:]
            s_rregion1=s_rregion[:2*i_lclip]

            print(s_lregion1, s_rregion1, s_clip, "left")
            b_tsd, s_TSD, i_refined_pos = self.ltsd.check_tsd(s_lregion1, s_rregion1, s_clip)
            if b_tsd == True:
                b_lclip1 = False
                s_tmp_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                self.m_tsd[s_tmp_id] = s_clip
        if b_tsd == False and i_rclip < global_values.LRD_MAX_TSD_LEN and i_rclip>=global_values.LRD_MIN_TSD_LEN:
            s_clip = algnmt.query_sequence[-1 * i_rclip:]
            s_lregion1 = s_lregion[-2 * i_rclip:]
            s_rregion1 = s_rregion[:2 * i_rclip]
            print(s_lregion1, s_rregion1, s_clip, "right")
            b_tsd, s_TSD, i_refined_pos = self.ltsd.check_tsd(s_lregion1, s_rregion1, s_clip)
            if b_tsd == True:
                b_rclip1 = False
                s_tmp_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                self.m_tsd[s_tmp_id] = s_clip
        return b_lclip1, b_rclip1

    def print_tsd(self):
        print(self.m_tsd)

    ####whether this is a qualified transduction
    ####1. hit the 3' end
    ####2. Or contain a full length
    def _is_qualified_transduction(self, i_map_end, i_cns_len):
        b_qualified=True
        if i_map_end < i_cns_len:
            return False
        return b_qualified

    ####if it is both side clipped, then it most contain a full length insertion
    def _is_qualfied_both_side_transduction(self, i_mapped, i_cns_len):
        if i_mapped<i_cns_len:
            return False
        else:
            return True

    ####
    # Input is the algnmt that only have one side clip (the other side either small or no clip)
    # Rule-1:
    # if right clipped, then supplementary aligned position should downstream the clip_pos
    # if left clipped, then supplementary aligned position should upstream the clip-pos
    # Rule-2:
    # supplementary aligned part should from the clipped part (allow very small size overlap)
    # Rule-3:
    # The supplementary part should have major portion aligned
    # 1.For internal deletion
    ####chr14~66616290~0	16	LINE1	3453	60	2001S111M5I6M4D69M3D294M8D304M5D10M2D2M3I42M1D130M
    ####SA:Z:LINE1,1311,-,353S1512M34D1112S,60,273;
    # 2. For inversion
    ####chr3~88973180~0	0	LINE1	1986	60	1745S590M5D990M1I764M
    ####SA:Z:LINE1,4335,-,2345S1728M1I16S,60,15;
    # 3. For inversion with internal deletion:
    # chr4~18161065~0	0	LINE1	3627	60	604S886M8S
    # SA:Z:LINE1,5436,-,891S606M1S,60,2;
    # 4. How about transduction happen together ??????????????????????
    ####
    # l_transduct is used to save potential transduction events that require further checking
    # l_succeed is used to save the succeeded ones in format: [[seg1_start, seg1_end, dir],[...]]
    def classify_with_supplementary_algnmt(self, algnmt, b_lclip, b_rclip, s_rep_type, l_transduct, l_succeed,
                                           i_cns_end=0, l_non_tprt=None):
        # For primary alignment, we parse out: 1. the (start, end) on both reference and the contig
        # 2. clip position, reverse complementary
        l_cigar = algnmt.cigar
        map_pos = algnmt.reference_start  # the mapping position
        query_seq = algnmt.query_sequence
        query_seq_len=len(query_seq)
        query_name = algnmt.query_name
        b_rc = algnmt.is_reverse

        if algnmt.mapping_quality < global_values.LRD_MIN_MAPQ:
            print("Mapping quality is low")
            return
        ####Hard clip are not considered
        #because we use the sequence length in the following steps
        if l_cigar[0][0]==5 or l_cigar[-1][0]==5:
            print("Hard clip is skipped")
            return
####
        # the primary aligned segment
        # find the mapped interval and mapped length
        # mapped interval: [map_pos, i_map_end]
        # mapped length: n_mapped
        i_seq_start = 0
        if l_cigar[0][0] == 4:
            i_seq_start = l_cigar[0][1]
        i_map_end, n_mapped, i_seq_end = self._get_map_interval(l_cigar, map_pos)
        s_rc = "+"
        if b_rc == True:
            s_rc = "-"
        t_segmt_pri = (map_pos, i_map_end, s_rc, i_seq_start, i_seq_end)
        t_segmt_sup = None
        if algnmt.has_tag("SA") == False:  # no supplementary alignment
            print("No supplementary alignment found")
            if self._is_qualified_transduction(i_map_end, i_cns_end)==False:
                return
            if b_lclip==True and b_rclip==True:
                ##here assume only full copy is qualified for both-side-transduction
                if self._is_qualfied_both_side_transduction(n_mapped, i_cns_end)==False:
                    return
            s_transduction_type, i_5mer_chk_len, i_3mer_chk_len=self._get_transduct_info_from_primary_algnmt(
                b_lclip, b_rclip, b_rc, l_cigar)
            l_transduct.append(
                (t_segmt_pri, t_segmt_sup, algnmt, s_transduction_type, i_5mer_chk_len, i_3mer_chk_len))
            return
####
        ####return value in format: "LINE1,6067,+,1382S23M1942S,0,0;"
        ####rname ,pos ,strand ,CIGAR ,mapQ ,NM;
        s_tag = algnmt.get_tag("SA")
        rcd_tag = self._parse_SA_fields(s_tag)
        i_sa_pos = rcd_tag[0]
        b_sa_rc = rcd_tag[1]
        s_sa_cigar = rcd_tag[2]
        i_sa_mapq = rcd_tag[3]  # this will not be checked
        l_sa_cigar = self._cvt_scigar_lcigar(s_sa_cigar)

        ####
        s_sa_rc = "+"
        if b_sa_rc == True:
            s_sa_rc = "-"
        # find the mapped interval of the supplementary region
        # sa mapped interval: [i_sa_pos, sa_ref_pos]
        # sa mapped bases: n_sa_mapped

        i_sa_seq_start = 0
        if l_sa_cigar[0][0] == 4:
            i_sa_seq_start = l_sa_cigar[0][1]
        sa_ref_pos, n_sa_mapped, i_sa_seq_end = self._get_map_interval(l_sa_cigar, i_sa_pos)
        if b_rc is not b_sa_rc:#different direction
            print("Different orientation")
            i_sa_seq_start = 0
            if l_sa_cigar[-1][0]==4:
                i_sa_seq_start = l_sa_cigar[-1][1]
            n_sa_mapped2, i_sa_seq_end=self._get_seq_interval_reverse_order(l_sa_cigar)
####
        # the regions that cover the repeat consensus by supplementary alignment
        t_segmt_sup = (i_sa_pos, sa_ref_pos, s_sa_rc, i_sa_seq_start, i_sa_seq_end)
        print(algnmt.query_name, map_pos, i_sa_pos, b_lclip, b_rclip,  i_seq_start, i_seq_end, i_sa_seq_start, i_sa_seq_end)
        i_pri_sup_max_end=i_map_end
        if sa_ref_pos >  i_pri_sup_max_end:
            i_pri_sup_max_end=sa_ref_pos
        ####
        b_same_ori= (b_rc is b_sa_rc)
        b_sa_qualified=self._is_qualified_SA(query_seq, map_pos, i_sa_pos, b_same_ori, b_lclip, b_rclip,
                                             i_seq_start, i_seq_end, i_sa_seq_start, i_sa_seq_end)

        if b_sa_qualified==False: # first check whether it is qualified "SA"
            print("Un-qualified SA")
            ####save the two flank region to flank for transduction checking
            s_transduction_type, i_5mer_chk_len, i_3mer_chk_len = self._get_transduct_info_from_primary_algnmt(
                b_lclip, b_rclip, b_rc, l_cigar)
            print(algnmt.query_name, i_5mer_chk_len, i_3mer_chk_len)
            if self._is_qualified_transduction(i_pri_sup_max_end, i_cns_end)==True:
                if b_lclip == True and b_rclip == True:
                    if self._is_qualfied_both_side_transduction(n_mapped, i_cns_end) == False:
                        return
                l_transduct.append(
                    (t_segmt_pri, None, algnmt, s_transduction_type, i_5mer_chk_len, i_3mer_chk_len))
        else:#SA is qualified
            print("Qualified SA")
            print(algnmt.query_name, map_pos, i_sa_pos, t_segmt_pri, t_segmt_sup)
            #find out the new left and right clip seq
            i_new_start=i_seq_start
            i_new_end=i_sa_seq_end
            if i_seq_start>i_sa_seq_start:#sa at upstream
                i_new_start=i_sa_seq_start
                i_new_end=i_seq_end
            i_new_5mer_clip=i_new_start #5mer clip length
            i_new_3mer_clip=query_seq_len-i_new_end #3mer clip length

            ####
            b_3mer_clip=True
            if i_new_3mer_clip<=global_values.LRD_PRI_SUP_MAX_CLIP:
                b_3mer_clip=False
            b_5mer_clip=True
            if i_new_5mer_clip <= global_values.LRD_PRI_SUP_MAX_CLIP:
                b_5mer_clip=False

            ####check for TSD
            self.check_TSD_for_site(algnmt, i_new_5mer_clip, i_new_3mer_clip, b_5mer_clip, b_3mer_clip)
            if b_3mer_clip==False and b_5mer_clip==False:#then this is a good case
                rslt_rcd = self.get_ins_information_from_algnmt_ori_seq(algnmt, s_rep_type, t_segmt_pri, t_segmt_sup)
                if map_pos > i_sa_pos:
                    rslt_rcd = self.get_ins_information_from_algnmt_ori_seq(algnmt, s_rep_type, t_segmt_sup, t_segmt_pri)
                #here check whether hit the end of consensus
                if i_pri_sup_max_end<i_cns_end:#doesn't reach the end, then filter out
                    l_non_tprt.append(rslt_rcd)
                else:
                    l_succeed.append(rslt_rcd)
            elif b_5mer_clip==True and b_3mer_clip==False:#then 5'-side clip
                s_tmp_type = self._5mer
                i_5mer_chk_len = i_new_5mer_clip
                i_3mer_chk_len = 0
                if self._is_qualified_transduction(i_pri_sup_max_end, i_cns_end) == False:
                    return
                if map_pos > i_sa_pos:
                    l_transduct.append(
                        (t_segmt_sup, t_segmt_pri, algnmt, s_tmp_type, i_5mer_chk_len, i_3mer_chk_len))
                else:
                    l_transduct.append(
                        (t_segmt_pri, t_segmt_sup, algnmt, s_tmp_type, i_5mer_chk_len, i_3mer_chk_len))
            elif b_5mer_clip==False and b_3mer_clip==True:#then 5'-side clip
                s_tmp_type = self._3mer
                i_5mer_chk_len = 0
                i_3mer_chk_len = i_new_3mer_clip
                if self._is_qualified_transduction(i_pri_sup_max_end, i_cns_end) == False:
                    return
                if map_pos > i_sa_pos:
                    l_transduct.append(
                        (t_segmt_sup, t_segmt_pri, algnmt, s_tmp_type, i_5mer_chk_len, i_3mer_chk_len))
                else:
                    l_transduct.append(
                        (t_segmt_pri, t_segmt_sup, algnmt, s_tmp_type, i_5mer_chk_len, i_3mer_chk_len))
            else:####both side clip
                s_tmp_type = self._both_transdct
                i_5mer_chk_len = i_new_5mer_clip
                i_3mer_chk_len = i_new_3mer_clip
                if self._is_qualified_transduction(i_pri_sup_max_end, i_cns_end) == False:
                    return
                if map_pos > i_sa_pos:
                    l_transduct.append(
                        (t_segmt_sup, t_segmt_pri, algnmt, s_tmp_type, i_5mer_chk_len, i_3mer_chk_len))
                else:
                    l_transduct.append(
                        (t_segmt_pri, t_segmt_sup, algnmt, s_tmp_type, i_5mer_chk_len, i_3mer_chk_len))
####
    ####this function is called in the l_transduction.py for second-level transduction calling
    def parse_sa_for_second_level_transduct(self, algnmt):
        l_cigar = algnmt.cigar
        map_pos = algnmt.reference_start  # the mapping position
        query_seq = algnmt.query_sequence
        query_seq_len = len(query_seq)
        query_name = algnmt.query_name
        b_rc = algnmt.is_reverse
        s_first_source=algnmt.reference_name
####
        ####Hard clip are not considered
        # because we use the sequence length in the following steps
        if l_cigar[0][0] == 5 or l_cigar[-1][0] == 5:
            return

        s_tag = algnmt.get_tag("SA")
        rcd_tag = self._parse_SA_fields(s_tag)
        i_sa_pos = rcd_tag[0]
        b_sa_rc = rcd_tag[1]
        s_sa_cigar = rcd_tag[2]
        i_sa_mapq = rcd_tag[3]  # this will not be checked
        s_second_source=rcd_tag[4]
        l_sa_cigar = self._cvt_scigar_lcigar(s_sa_cigar)
        s_sa_rc = "+"
        if b_sa_rc == True:
            s_sa_rc = "-"

        #primary interval
        i_seq_start = 0
        if l_cigar[0][0] == 4:
            i_seq_start = l_cigar[0][1]
        i_map_end, n_mapped, i_seq_end = self._get_map_interval(l_cigar, map_pos)
        s_rc = "+"
        if b_rc == True:
            s_rc = "-"
        t_segmt_pri = (s_first_source, map_pos, i_map_end, s_rc, i_seq_start, i_seq_end)
        b_src_same_ori=self._two_source_same_orientation(s_first_source, s_second_source)
        #sa interval
        ####Here make sure the two sources are of the same direction (reverse complementary???)
        i_sa_seq_start = 0
        if l_sa_cigar[0][0] == 4:
            i_sa_seq_start = l_sa_cigar[0][1]
        sa_ref_pos, n_sa_mapped, i_sa_seq_end = self._get_map_interval(l_sa_cigar, i_sa_pos)
        b_diff_dir=b_rc ^ b_sa_rc ^ b_src_same_ori
        if b_diff_dir:  # different direction
            i_sa_seq_start = 0
            if l_sa_cigar[-1][0] == 4:
                i_sa_seq_start = l_sa_cigar[-1][1]
            n_sa_mapped2, i_sa_seq_end = self._get_seq_interval_reverse_order(l_sa_cigar)#
        t_segmt_sup = (s_second_source, i_sa_pos, sa_ref_pos, s_sa_rc, i_sa_seq_start, i_sa_seq_end)
        return t_segmt_pri, t_segmt_sup
        ####make sure the two segements are not largely overlapped on the seq side

        ####make sure the overall alignment is qualified

    ####check TSD?
    def _two_source_same_orientation(self, s_src1, s_src2):
        fields1 = s_src1.split(global_values.SEPERATOR)
        b_rc1=False
        if fields1[-1][0]=="1":
            b_rc1=True
        fields2 = s_src2.split(global_values.SEPERATOR)
        b_rc2 = False
        if fields2[-1][0] == "1":
            b_rc2 = True
        return b_rc1 is b_rc2

    ####
    def _get_transduct_info_from_primary_algnmt(self, b_lclip, b_rclip, b_rc, l_cigar):
        i_5mer_chk_len = 0
        i_3mer_chk_len = 0
        if b_lclip == True:
            i_5mer_chk_len = l_cigar[0][1]
        if b_rclip == True:
            i_3mer_chk_len = l_cigar[-1][1]
        # if b_rc == True:  # it is reverse complementary
        #     if b_lclip == True:
        #         i_3mer_chk_len = l_cigar[0][1]
        #     if b_rclip == True:
        #         i_5mer_chk_len = l_cigar[-1][1]
####
        s_transduction_type = self._3mer
        if b_lclip == True and b_rclip == True:
            s_transduction_type = self._both_transdct
        elif i_5mer_chk_len > 0:
            s_transduction_type = self._5mer
        return s_transduction_type, i_5mer_chk_len, i_3mer_chk_len

    ####
    def _is_A_dominant(self, s_seq):
        pa=PolyA()
        if pa.is_dominant_A(s_seq, global_values.LRD_PRI_SUP_DOMINANT_POLYA_RATIO):
            return True
        return False

    def _is_qualified_SA(self, s_seq, map_pos, i_sa_pos, b_same_ori, b_lclip, b_rclip,  i_seq_start, i_seq_end,
                         i_sa_start, i_sa_end):
        b_qualfied=True
        ####
        # check rule-1: if right (left) clipped, then supplementary aligned at downstream (upstream)
        #1745S590M5D990M1I764M
        #SA:Z:LINE1,4335,-,2345S1732M1I12S
        # b_sa_downstream = True
        # if i_sa_pos < map_pos:
        #     b_sa_downstream = False
        # if b_lclip==True and b_rclip==False:
        #     if b_sa_downstream == True and b_same_ori==True:
        #         #print "debug 1"
        #         return False
        #     elif b_sa_downstream == False and b_same_ori==False:
        #         #print "debug 1.1"
        #         return False ####
        # if b_rclip == True and b_lclip==False:
        #     if b_sa_downstream == False and b_same_ori==True:
        #         #print "debug 2"
        #         return False
        #     elif b_sa_downstream == True and b_same_ori==False:
        #         #print "debug 2.1"
        #         return False
####
        #the two segments should not be overlap, or overlap region is short
        #if contained, then not allowd
        if i_seq_start<i_sa_start and i_seq_end> i_sa_end:
            print("Debug 3: One is contained in another: {0}".format(map_pos))
            return False
        if i_sa_start < i_seq_start and i_sa_end>i_seq_end:
            print("Debug 4: One is contained in another: {0}".format(map_pos))
            return False
        #if overlap, then overlap region is short:
        if i_seq_end>i_sa_start and i_seq_end<i_sa_end:#
            if (i_seq_end-i_sa_start)>global_values.LRD_PRI_SUP_MAX_OVRLAP:
                print("Debug 5: Large overlap between the two segment: {0}".format(map_pos))
                return False
        if i_sa_end>i_seq_start and i_sa_end<i_seq_end:
            if (i_sa_end-i_seq_start)>global_values.LRD_PRI_SUP_MAX_OVRLAP:
                print("Debug 6: Large overlap between the two segment: {0}".format(map_pos))
                return False

        #if there is a large gap between them, then it is also not qualified
        if (i_seq_start<i_sa_start) and (abs(i_sa_start-i_seq_end)>global_values.LRD_PRI_SUP_MAX_GAP_IN_SEQ):
            b_a_dominant=False
            if (i_seq_end<i_sa_start) and (self._is_A_dominant(s_seq[i_seq_end:i_sa_start])):
                b_a_dominant=True
            if b_a_dominant==False:
                print("Debug 7: Large gap between the two segment: {0}".format(map_pos))
                return False
            #elif abs(i_sa_start-i_seq_end)>global_values.LRD_PRI_SUP_MAX_GAP_IN_SEQ_WITH_A_DOMNT:
            elif self._is_sa_most_part_mapped(abs(i_sa_start-i_seq_end), (i_sa_end-i_sa_start))==False:
                return False
        elif (i_seq_start>i_sa_start) and (abs(i_seq_start-i_sa_end)>global_values.LRD_PRI_SUP_MAX_GAP_IN_SEQ):
            b_a_dominant = False
            if (i_seq_start>i_sa_end) and (self._is_A_dominant(s_seq[i_sa_end:i_seq_start])):
                b_a_dominant=True
            if b_a_dominant==False:
                print("Debug 8: Large gap between the two segment: {0}".format(map_pos))
                return False
            #elif abs(i_seq_start-i_sa_end)>global_values.LRD_PRI_SUP_MAX_GAP_IN_SEQ_WITH_A_DOMNT:
            elif self._is_sa_most_part_mapped(abs(i_seq_start-i_sa_end), (i_sa_end-i_sa_start))==False:
                return False
        return b_qualfied

####
    def _is_sa_most_part_mapped(self, i_gap, i_mapped):
        if float(i_gap)/float(i_mapped) < global_values.LRD_PRI_SUP_MAX_GAP_MAP_RATIO:
            return True
        return False
    ####save the non transduction inseriton events to a file
    ####each record in format:
    # ins_chrm, ins_pos, rep_type, structure_type, s_region1, s_region2, TSD, transduct_src, s_seq, transduct_seq
    ####i_del_slack: if the interval distance > i_del_slack, then view as an internal insertion
    def save_non_transduct_ins_info_to_file(self, l_succeed, i_del_slack, sf_out, i_ofst=0):
        with open(sf_out, "w") as fout_ins:
            lis = LInternalStructure()
            xpolyA=PolyA()
            for rcd in l_succeed:
                s_rg1=rcd[4]#region1
                s_rg2=rcd[5]#region2

                (ins_chrm, ins_pos, rep_type, structure_type, s_region1, s_region2, s_TSD, transduct_src,
                 s_seq, transduct_seq)=rcd
                s_structure=lis.get_internal_structure(s_rg1, s_rg2, i_del_slack)
                s_with_polyA = self.chk_polyA(xpolyA, s_seq)
                # b_rpolyA=xpolyA.is_consecutive_polyA_T(seq[:10])
                # b_lpolyA=xpolyA.is_consecutive_polyA_T(seq[-10:])
                # if b_lpolyA or b_rpolyA:
                #     s_with_polyA="with_polyA"
                s_id = "{0}~{1}".format(ins_chrm, ins_pos)
                if s_id in self.m_tsd:
                    s_TSD=self.m_tsd[s_id]

                if len(s_seq)<global_values.LRD_MIN_INS_LTH:
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
####
    def is_no_polA_no_TSD(self, s_polyA, s_TSD):
        if (s_polyA == self._s_no_polyA) and (s_TSD == "None"):
            return True
        return False

    def chk_polyA(self, xpolyA, seq):
        s_with_polyA = self._s_no_polyA
        if (seq is None) or (len(seq)<10):
            return s_with_polyA
        b_rpolyA = xpolyA.is_consecutive_polyA_T(seq[:25])
        b_lpolyA = xpolyA.is_consecutive_polyA_T(seq[-25:])
        if b_lpolyA or b_rpolyA:
            s_with_polyA = self._s_with_polyA
        return s_with_polyA

    def chk_polyA_large_region(self, xpolyA, seq):
        s_with_polyA = self._s_no_polyA
        if (seq is None) or (len(seq)<10):
            return s_with_polyA
        b_rpolyA = xpolyA.is_consecutive_polyA_T2(seq[:100])#here request 6A or 6T
        b_lpolyA = xpolyA.is_consecutive_polyA_T2(seq[-100:])
        if b_lpolyA or b_rpolyA:
            s_with_polyA = self._s_with_polyA
        return s_with_polyA
####
    # select and save the transduction information to file
####Note: 5' or 3' transduction should adjust accordingly based on whether it is reverse complementary or not
    def slct_save_transduction_to_file(self, l_transduct, m_candidates, m_with_polya, rep_type, sf_out):
        m_not_confirmed={}
        with open(sf_out, "w") as fout_rslt:
            for rcd in l_transduct:
                algnmt = rcd[2]
                s_seq = algnmt.query_sequence  ###this is already converted back to consensus
                s_transduction_type = rcd[3]
                (ins_chrm, ins_pos) = self.lbinfo.get_ins_info_from_qname(algnmt.query_name)
                sinfo_5mer = ""
                sinfo_3mer = ""
                b_rc=algnmt.is_reverse

                if s_transduction_type == self._5mer or s_transduction_type == self._both_transdct:
                    sinfo_5mer = "{0}~{1}~{2}".format(ins_chrm, ins_pos, self._5mer)
                if s_transduction_type == self._3mer or s_transduction_type == self._both_transdct:
                    sinfo_3mer = "{0}~{1}~{2}".format(ins_chrm, ins_pos, self._3mer)
                s_id = "{0}~{1}".format(ins_chrm, ins_pos)

                s_5mer_source = "None"
                s_3mer_source = "None"
                s_TSD = "None"
                if s_id in self.m_tsd:
                    s_TSD=self.m_tsd[s_id]
                b_with_polyA = False
                s_5mer_seq = "None"
                s_3mer_seq = "None"

                ins_info_rcd = self.get_ins_information_from_algnmt_ori_seq(algnmt, rep_type, rcd[0], rcd[1])
                if (sinfo_5mer not in m_candidates) and (sinfo_3mer not in m_candidates): # not 3' or 5' transduction
                    tmp_rcd = (ins_info_rcd, s_5mer_source, s_3mer_source, s_TSD,
                               b_with_polyA, s_seq, s_5mer_seq, s_3mer_seq)
                    if s_id not in m_with_polya:
                        print("Error happen at [_slct_save_transduction_to_file]: " \
                              "{0} doesn't have hit in m_with_polya\n".format(s_id))
                        m_not_confirmed[s_id] = tmp_rcd
                    if m_with_polya[s_id] == True:  # no hit on flanking regions, but has polyA signal
                        m_not_confirmed[s_id]=tmp_rcd
                        # self._save_transduction_info_to_file(ins_info_rcd, s_5mer_source, s_3mer_source, s_TSD,
                        #                                      b_with_polyA, s_seq, s_5mer_seq, s_3mer_seq, fout_rslt)#
                    continue
                if (sinfo_5mer in m_candidates) and (sinfo_3mer in m_candidates):  # both 5'mer and 3'mer transduction
                    if s_id not in m_with_polya:
                        print("Error happens at [_slct_save_transduction_to_file]: " \
                              "{0} doesn't have hit in m_with_polya\n".format(s_id))
                    if m_with_polya[s_id] == True:
                        b_with_polyA = True

                    s_5mer_source = m_candidates[sinfo_5mer][0]
                    s_3mer_source = m_candidates[sinfo_3mer][0]
###################
                    if s_5mer_source is not s_3mer_source:
                        print("{0} has two transduction, but has different source: {1}  and {2}".format(s_id, s_5mer_source, s_3mer_source))
                    s_5mer_seq=m_candidates[sinfo_5mer][1]
                    s_3mer_seq=m_candidates[sinfo_3mer][1]
                elif sinfo_5mer in m_candidates:
                    s_5mer_source = m_candidates[sinfo_5mer][0]
                    s_5mer_seq = m_candidates[sinfo_5mer][1]
                elif sinfo_3mer in m_candidates:
                    s_3mer_source = m_candidates[sinfo_3mer][0]
                    s_3mer_seq = m_candidates[sinfo_3mer][1]

                if len(s_seq) < global_values.LRD_MIN_INS_LTH:
                    continue

                i_5mer_map_len=0
                if s_5mer_seq!="None":
                    i_5mer_map_len=len(s_5mer_seq)
                i_3mer_map_len = 0
                if s_3mer_seq != "None":
                    i_3mer_map_len = len(s_3mer_seq)
                if self._is_qualified_transduction_af_realign(
                        len(s_seq), rcd[4], rcd[5], i_5mer_map_len, i_3mer_map_len, self._f_min_map_ratio)==False:
                    continue

                if b_rc==False:
                    self._save_transduction_info_to_file(ins_info_rcd, s_5mer_source, s_3mer_source, s_TSD,
                                                         b_with_polyA, s_seq, s_5mer_seq, s_3mer_seq,  fout_rslt)
                else:
                    self._save_transduction_info_to_file(ins_info_rcd, s_3mer_source, s_5mer_source, s_TSD,
                                                         b_with_polyA, s_seq, s_3mer_seq, s_5mer_seq, fout_rslt)
        return m_not_confirmed
####

    def _is_qualified_transduction_af_realign(self, n_all_len, n_lclip, n_rclip, n_ltd_map, n_rtd_map, f_min_ratio):
        if float(n_lclip+n_rclip-n_ltd_map-n_rtd_map)/float(n_all_len) > f_min_ratio:
            return False
        else:
            return True
    ####
    # save a record to file (already opened)
    # each output record in format:
    # chrm pos rep_type structure_type region1 region2 TSD transduction-source, seq, transduction-seq
    ##transdct_rcd in format: (t_segmt_pri, t_segmt_sup, algnmt, s_transduction_type, i_5mer_chk_len, i_3mer_chk_len)
    def _save_transduction_info_to_file(self, ins_info_rcd, s_5mer_source, s_3mer_source, s_TSD,
                                        b_with_polyA, s_seq, s_5mer_seq, s_3mer_seq, fout, i_ofst=0):
        ins_chrm=ins_info_rcd[0]
        ins_pos=ins_info_rcd[1]
        s_rep_type=ins_info_rcd[2]
        s_rgn1=ins_info_rcd[4]
        s_rgn2 = ins_info_rcd[5]

        lis = LInternalStructure()
        s_structure = lis.get_internal_structure(s_rgn1, s_rgn2, global_values.LRD_MIN_INTERNAL_DEL)
        s_with_polyA = "clip_seq_with_polyA"
        if b_with_polyA == True:
            s_with_polyA = "clip_seq_no_polyA"
        if i_ofst!=0:
            s_rgn1=self.adjust_pos(s_rgn1, i_ofst)
            s_rgn2 = self.adjust_pos(s_rgn2, i_ofst)

        sinfo = ins_chrm + "\t" + str(
            ins_pos) + "\t" + s_rep_type + "\t" + s_structure + "\t" + s_rgn1 + "\t" + s_rgn2 + "\t" + \
                s_TSD + "\t" + s_5mer_source + "\t" + s_3mer_source + "\t" + s_seq + "\t" + s_5mer_seq + \
                "\t" + s_3mer_seq + "\t" + s_with_polyA + "\n"
        fout.write(sinfo)

####
    def save_tbd_cases_to_file(self, m_not_slct, sf_out):
        with open(sf_out,"w") as fout:
            for s_id in m_not_slct:
                tmp_rcd=m_not_slct[s_id]
                (ins_info_rcd, s_5mer_source, s_3mer_source, s_TSD, b_with_polyA, s_seq, s_5mer_seq, s_3mer_seq) = tmp_rcd
                self._save_transduction_info_to_file(ins_info_rcd, s_5mer_source, s_3mer_source, s_TSD,
                                                         b_with_polyA, s_seq, s_5mer_seq, s_3mer_seq, fout)
####
    def get_ins_pos_from_query_name(self, algnmt):
        query_name = algnmt.query_name
        name_fields = query_name.split(global_values.SEPERATOR)
        ins_chrm = name_fields[0]
        ins_pos = name_fields[1]
        return ins_chrm, ins_pos

####
    # (3964, 5520, '-', 555, 2111)
    # (5498, 6054, '+', 0, 556)
    # assume there are at most 2 segments
    def get_ins_information_from_algnmt_ori_seq(self, algnmt, rep_type, t_seg1, t_seg2):
        query_seq = algnmt.query_sequence
        query_name = algnmt.query_name
        name_fields = query_name.split(global_values.SEPERATOR)
        ins_chrm = name_fields[0]
        ins_pos = name_fields[1]

        i_start = -1
        i_end = -1
        s_region1 = self.lbinfo.get_s_NONE()
        s_region2 = self.lbinfo.get_s_NONE()
        if t_seg1 is not None:
            i_start = t_seg1[3]
            i_end = t_seg1[4]
            s_region1 = "{0}:{1}:{2}".format(t_seg1[0], t_seg1[1], t_seg1[2])
        if t_seg2 is not None:
            s_region2 = "{0}:{1}:{2}".format(t_seg2[0], t_seg2[1], t_seg2[2])
            # if i_start == -1:
            #     i_start = t_seg2[3]
            # i_end = t_seg2[4]
            if (t_seg2[3]<i_start) or (i_start==-1):
                i_start=t_seg2[3]
            if t_seg2[4]>i_end:
                i_end=t_seg2[4]
####
        s_seq = query_seq[i_start:i_end]
        structure_type = self.lbinfo.get_s_NONE()
        TSD = self.lbinfo.get_s_NONE()
        transduct_src = global_values.NOT_TRANSDUCTION
        transduct_seq = self.lbinfo.get_s_NONE()
        t_rslt = (
            ins_chrm, ins_pos, rep_type, structure_type, s_region1, s_region2, TSD, transduct_src, s_seq,
            transduct_seq)
        return t_rslt
####
####
    ####check whether it is fully mapped
    ####If both side are fully mapped, then save to "l_both_transduct"
    def is_fully_map(self, algnmt, i_max_clip_len):
        l_cigar = algnmt.cigar
        n_lclip = 0
        n_rclip = 0
        if l_cigar[0][0] == 4:  # left-clip
            n_lclip = l_cigar[0][1]
        if l_cigar[-1][0] == 4:  # right clipped
            n_rclip = l_cigar[-1][1]

        if l_cigar[0][0] == 5 and l_cigar[-1][0] == 5:#both side hard clip
            return False, l_cigar[0][1], l_cigar[-1][1]

        if (n_lclip + n_rclip) < i_max_clip_len:
            return True, n_lclip, n_rclip
        else:
            return False, n_lclip, n_rclip

    ####return value in format: "SA:Z:LINE1,6067,+,1382S23M1942S,0,0;"
    ####rname ,pos ,strand ,CIGAR ,mapQ ,NM;
    def _parse_SA_fields(self, s_tag):
        s_tag_fields = s_tag.split(",")
        s_ref=s_tag_fields[0]
        i_sa_pos = int(s_tag_fields[1])
        b_rc = True
        if s_tag_fields[2] == "+":
            b_rc = False
        s_sa_cigar = s_tag_fields[3]
        i_sa_mapq = int(s_tag_fields[4])
        rcd = (i_sa_pos, b_rc, s_sa_cigar, i_sa_mapq, s_ref)
        return rcd

    # get the mapped interval
    def _get_map_interval(self, l_cigar, i_map_pos):
        ref_pos = i_map_pos
        seq_pos = 0
        n_mapped = 0
        for (opn, lth) in l_cigar:
            if opn == 0:  # alignment match (M)
                ref_pos += lth
                seq_pos += lth
                n_mapped += lth  # accumulate the mapped length
            elif opn == 1:  # insertion
                seq_pos += lth
            elif opn == 2 or opn == 3:  # deletion (D) or skipped region (N)
                ref_pos += lth
            elif opn == 4:  # soft-clip (S)
                seq_pos += lth
            elif opn == 5 or opn == 6:  # hard-clip (H) or padding (P)
                seq_pos += lth
            elif opn == 7 or opn == 8:  # sequence match (=) or sequence mismatch (X)
                ref_pos += lth
                seq_pos += lth
                n_mapped += lth
        if len(l_cigar)>1 and (l_cigar[-1][0]==4 or l_cigar[-1][0]==5):
            seq_pos-=l_cigar[-1][1]
        # map position
        return ref_pos, n_mapped, seq_pos

    # get the mapped interval
    def _get_seq_interval_reverse_order(self, l_cigar):
        seq_pos = 0
        n_mapped = 0
        for (opn, lth) in reversed(l_cigar):
            if opn == 0:  # alignment match (M)
                seq_pos += lth
                n_mapped += lth  # accumulate the mapped length
            elif opn == 1:  # insertion
                seq_pos += lth
            elif opn == 4:  # soft-clip (S)
                seq_pos += lth
            elif opn == 5 or opn == 6:  # hard-clip (H) or padding (P)
                seq_pos += lth
            elif opn == 7 or opn == 8:  # sequence match (=) or sequence mismatch (X)
                seq_pos += lth
                n_mapped += lth
        if len(l_cigar)>1 and (l_cigar[0][0]==4 or l_cigar[0][0]==5):
            seq_pos-=l_cigar[0][1]
        # map position
        return n_mapped, seq_pos

####
    ####
    def _cvt_scigar_lcigar(self, s_cigar):
        print(s_cigar)
        m_op = {'S': 4, 'H': 5, 'M': 0, 'D': 2, 'I': 1, 'P': 6, 'X': 8, '=': 7, 'N': 3, 'B': 9}
        s_tmp = ""
        l_op = []
        for ch in s_cigar:
            if ch in m_op:
                i_op = m_op[ch]
                i_len = int(s_tmp)
                l_op.append((i_op, i_len))
                s_tmp = ''
            else:
                s_tmp += ch
        return l_op
        ####
####
    def merge_outputs(self, sf_non_transduct, sf_transduct, sf_merged):
        with open(sf_merged, "w") as fout_merged, open(sf_non_transduct) as fin_non_td, open(sf_transduct) as fin_td:
            for line in fin_non_td:
                fout_merged.write(line)
            for line in fin_td:
                fout_merged.write(line)

    def adjust_pos(self, s_rg, i_ofst):
        if s_rg is None:
            return s_rg
        if s_rg == "":
            return s_rg
        fields=s_rg.split(":")
        if len(fields)<3:
            return s_rg
        i_start=int(fields[0])-i_ofst
        i_end=int(fields[1])-i_ofst
        s_new="{0}:{1}:{2}".format(i_start, i_end, fields[2])
        return s_new
####