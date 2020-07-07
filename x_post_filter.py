####
##04/16/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

####1. for tprt both (with polyA and TSD) events, we filter out those:
        ####1.1 polyA dominant cases, this module will mainly used on Alu
                # (as for L1 and SVA, some middle region or transduction will contain the polyA part)
####2. false positive ones comes from low reads quality part caused clip
####
#For SVA: if fall in a reference SVA copy and only one-side-clip feature,
# then it's more likely because of the VNTR expansion caused FP, or at the boundary of SVA copy

#One side purly polyA clip, but has the other side clipped reads, but none aligned
#Two-side, purely polyA ones are filtered out
#Collect the transduction events into the final output
####

import global_values
from x_annotation import *

class XPostFilter():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1] != "/":
            self.swfolder += "/"
        self.n_jobs = n_jobs
        self.REP_TYPE_L1="LINE1"
        self.REP_TYPE_ALU="ALU"
        self.REP_TYPE_SVA="SVA"
        self.REP_TYPE_HERV="HERV"
        self.REP_TYPE_MIT="Mitochondria"
        self.REP_TYPE_MSTA="Msta"
        self.REP_ALU_MIN_LEN=200
        self.REP_ALU_POLYA_START = 255
        self.REP_LINE_MIN_LEN=300
        self.REP_LINE_POLYA_START = 5950
        self.REP_SVA_MIN_LEN = 300
        self.REP_SVA_CNS_HEAD=400
        self.REP_SVA_POLYA_START=1900
        self.REP_HERV_MIN_LEN = 300
        self.REP_MIT_MIN_LEN = 300
        self.REP_MSTA_MIN_LEN = 300
        self.b_rslt_with_chr=True #####if rslt has "chr" while rmsk doesn't (vice versa), then it will be a problem
        self._two_side="two_side"
        self._one_half_side = "one_half"
        self._one_side="one_side"
        self._other="other"
        self._boundary_extnd=100 #for rmsk copy extend length
########Hard code here!!!!!
        self.nclip_half_cutoff=5 #one-side is polyA, which the other side has more than cutoff unmapped clipped reads!
        self.sva_clip_cluster_diff_cutoff=200
        self.indel_reads_max_ratio=0.3
        self.abnormal_cov_times=2

####
    def post_processing_SVA(self, l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation, m_cutoff, f_cov, sf_new_out):
        i_dft_cluster_diff_cutoff=global_values.TWO_CLIP_CLUSTER_DIFF_CUTOFF
        global_values.set_two_clip_cluster_diff_cutoff(self.sva_clip_cluster_diff_cutoff)
        with open(sf_new_out, "w") as fout_new:
            for rcd in l_old_rcd:
                s_rep_supt_type=xtea_parser.get_ins_sub_type(rcd)
                ins_chrm = rcd[0]
                ins_pos = rcd[1]
                b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
                div_rate, sub_family, family, pos_start, pos_end=xannotation.get_div_subfamily(ins_chrm, i_pos)
                b_with_polyA=xtprt_filter.has_polyA_signal(rcd)
                if s_rep_supt_type is self._two_side:####two sides
                    if xtprt_filter.is_polyA_dominant_two_side_sva(rcd, self.REP_SVA_CNS_HEAD)==True:
                        continue
                    if af_filter.is_qualified_rcd(rcd[-1], m_cutoff) == False:
                        continue
                    if b_with_polyA==False:#require polyA
                        continue
                elif s_rep_supt_type is self._one_half_side:##one and half side
                    #check whether fall in repetitive region of the same type
                    if (b_in_rep is True) and (div_rate<global_values.REP_DIVERGENT_CUTOFF):
                        continue
                    if b_with_polyA==False:#require polyA
                        continue
                    if xtprt_filter.is_polyA_dominant_two_side_sva(rcd, self.REP_SVA_CNS_HEAD) is True:
                        continue
                elif s_rep_supt_type is self._one_side:###one side
                    if (b_in_rep is True) and (div_rate < global_values.REP_DIVERGENT_CUTOFF):
                        continue
                    if b_with_polyA == False:
                        continue
                    # if xtprt_filter.is_polyA_dominant_two_side_sva(rcd) is True:
                    #     continue
                else:#for other type, just skip
                    continue
####
                #if two-side-clip, and none-in-polyA side, then filter out
                if b_in_rep and xtprt_filter.is_two_side_clip_both_polyA_sva(rcd, self.REP_SVA_POLYA_START):
                    continue
                #if fall in rep region, and none of the clipped part hit polyA
                if b_in_rep and xtprt_filter.is_two_side_clip_both_non_polyA_sva(rcd, self.REP_SVA_POLYA_START):
                    continue

                #save the passed ones
                s_in_rep = "not_in_SVA_copy"
                if b_in_rep is True:
                    s_in_rep = "Fall_in_SVA_copy_"+str(div_rate)
                s_pass_info = rcd[-1].rstrip() + "\t" + s_in_rep + "\n"
                fout_new.write(s_pass_info)
        global_values.set_two_clip_cluster_diff_cutoff(i_dft_cluster_diff_cutoff)
####
####
    def post_processing_Alu(self, l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation, m_cutoff, f_cov, sf_new_out):
        with open(sf_new_out, "w") as fout_new:
            for rcd in l_old_rcd:
                s_rep_supt_type = xtea_parser.get_ins_sub_type_alu(rcd) #here if both un-clip, then view as "two-side"
                ins_chrm = rcd[0]
                ins_pos = rcd[1]
                b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
                div_rate, sub_family, family, pos_start, pos_end = xannotation.get_div_subfamily(ins_chrm, i_pos)
                b_with_polyA = xtprt_filter.has_polyA_signal(rcd)

                if s_rep_supt_type is self._two_side:  ####two sides
                    if xtprt_filter.is_polyA_dominant_two_side(rcd) == True:
                        continue
                    if af_filter.is_qualified_rcd(rcd[-1], m_cutoff) == False:
                        continue
                    if (b_with_polyA is False) and \
                                    xtprt_filter.is_two_side_tprt_and_with_polyA(rcd, self.REP_ALU_MIN_LEN)==False:
                        continue
                elif s_rep_supt_type is self._one_half_side:  ##one and half side
                    # check whether fall in repetitive region of the same type
                    #b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
                    if (b_in_rep is True) and (div_rate < global_values.REP_DIVERGENT_CUTOFF):
                        continue
                    if xtprt_filter.is_polyA_dominant_one_side(rcd, self.nclip_half_cutoff) is True:
                        continue
                    if (b_in_rep is True) and (b_with_polyA==False):
                        continue
                elif s_rep_supt_type is self._one_side:  ###one side
                    #b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
                    if (b_in_rep is True) and (div_rate < global_values.REP_DIVERGENT_CUTOFF):
                        continue
                    if xtprt_filter.is_polyA_dominant_one_side(rcd, self.nclip_half_cutoff) is True:
                        continue
                    if (b_in_rep is True) and (b_with_polyA==False):
                        continue
                else:  # for other type, just skip
                    continue

                # if have both side clipped reads, and fall in Alu copy,
                # and both side clip reads aligned to end of consensus, then skip
                if (b_in_rep is True) and \
                                xtprt_filter.is_two_side_clip_and_both_hit_end(rcd, self.REP_ALU_POLYA_START)==True:
                    continue

                # if the nearby region has lots of indel reads, then skip
                b_cov_abnormal=xtprt_filter.cov_is_abnormal(rcd, f_cov*self.abnormal_cov_times)
                n_indel_cutoff=int(f_cov*self.indel_reads_max_ratio)
                if xtprt_filter.has_enough_indel_reads(rcd, n_indel_cutoff)==True \
                        and n_indel_cutoff>0 and b_cov_abnormal==True:
                    continue
####
                # save the passed ones
                # l_slct.append(rcd)
                s_in_rep="not_in_Alu_copy"
                if b_in_rep is True:
                    s_in_rep = "Fall_in_Alu_copy_"+str(div_rate)
                s_pass_info = rcd[-1].rstrip() + "\t" + s_in_rep +"\n"
                fout_new.write(s_pass_info)
####
####
    def post_processing_L1(self, l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation, m_cutoff, f_cov, sf_new_out):
        with open(sf_new_out, "w") as fout_new:
            for rcd in l_old_rcd:
                s_rep_supt_type = xtea_parser.get_ins_sub_type(rcd)
                ins_chrm = rcd[0]
                ins_pos = rcd[1]
                b_in_rep, i_pos = xannotation.is_within_repeat_region_interval_tree(ins_chrm, int(ins_pos))
                div_rate, sub_family, family, pos_start, pos_end = xannotation.get_div_subfamily(ins_chrm, i_pos)
                b_with_polyA = xtprt_filter.has_polyA_signal(rcd)
                if b_with_polyA is False:
                    continue
                if s_rep_supt_type is self._two_side:####two sides
                    if xtprt_filter.is_polyA_dominant_two_side(rcd) == True:
                        continue
                    if xtprt_filter.hit_end_of_cns(rcd)==False:
                        continue
                    if af_filter.is_qualified_rcd(rcd[-1], m_cutoff) == False:
                        continue
                elif s_rep_supt_type is self._one_half_side:  ##one and half side
                    # check whether fall in repetitive region of the same type
                    if (b_in_rep is True) and (div_rate < global_values.REP_DIVERGENT_CUTOFF):
                        continue
                        # if xtprt_filter.is_polyA_dominant_one_side(rcd, self.nclip_half_cutoff) is True:
                        #     continue
                elif s_rep_supt_type is self._one_side:  ###one side
                    if (b_in_rep is True) and (div_rate < global_values.REP_DIVERGENT_CUTOFF):
                        continue
                        # if xtprt_filter.is_polyA_dominant_one_side(rcd, self.nclip_half_cutoff) is True:
                        #     continue
                else:# for other type, just skip
                    continue

                #if coverage is abnormal (>2 times coverage) or within rep region, but two sides clipped are from polyA
                b_cov_abnormal = xtprt_filter.cov_is_abnormal(rcd, f_cov * self.abnormal_cov_times)
                if ((b_in_rep is True) or (b_cov_abnormal is True)) and \
                                xtprt_filter.is_two_side_clip_and_both_hit_end(rcd, self.REP_LINE_POLYA_START)==True:
                    continue

                #if within L1 repeat, but also neither clip side hit the end, then skip
                #this is initially added for 10X
                if (b_in_rep is True) and (xtprt_filter.hit_L1_tail(rcd, self.REP_LINE_POLYA_START)==False):
                    continue

                # save the passed ones
                # l_slct.append(rcd)
                s_in_rep = "not_in_LINE1_copy"
                if b_in_rep is True:
                    s_in_rep = "Fall_in_LINE1_copy_"+str(div_rate)
                s_pass_info = rcd[-1].rstrip() + "\t" + s_in_rep + "\n"
                fout_new.write(s_pass_info)
    ####
    # def load_basic_info_samples(self, working_folder):
    #     sf_basic_info = working_folder + global_values.BASIC_INFO_FILE
    #     if os.path.isfile(sf_basic_info)==False:
    #         return None
    #     m_basic_info=self.load_basic_info_from_file(sf_basic_info)
    #     return m_basic_info
    #
#####
    def run_post_filtering(self, sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, f_cov, sf_new_out):
        #1. for two side cases, filter out the polyA dominant ones
        #2. for single-end cases, filter out those fall in same type of reference copies
        xtea_parser=XTEARsltParser()
        xtea_parser.set_rep_support_type(self._two_side, self._one_half_side, self._one_side, self._other)
        l_old_rcd=xtea_parser.load_in_xTEA_rslt(sf_xtea_rslt)
        xannotation = self.construct_interval_tree(sf_rmsk, i_min_copy_len, self.b_rslt_with_chr)
        xtprt_filter=XTPRTFilter(self.swfolder, self.n_jobs)
        af_filter = AFConflictFilter(self.swfolder, self.n_jobs)
        l_types = af_filter.get_rep_type()
        m_cutoff = af_filter.get_cutoff_by_type(l_types)

####Hard code here !!!!!!!!
        if i_rep_type & 1 is not 0:
            self.post_processing_L1(l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation,
                                    m_cutoff, f_cov, sf_new_out)
        elif i_rep_type & 4 is not 0:
            self.post_processing_SVA(l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation,
                                     m_cutoff, f_cov, sf_new_out)
        else:
            self.post_processing_Alu(l_old_rcd, xtea_parser, xtprt_filter, af_filter, xannotation,
                                     m_cutoff, f_cov, sf_new_out)
####
####
    ####parse out the specific repeat type
    def parse_rep_type(self, i_rep_type):
        l_rep_type = []
        if i_rep_type & 1 != 0:
            l_rep_type.append(self.REP_TYPE_L1)
        if i_rep_type & 2 != 0:
            l_rep_type.append(self.REP_TYPE_ALU)
        if i_rep_type & 4 != 0:
            l_rep_type.append(self.REP_TYPE_SVA)
        if i_rep_type & 8 != 0:
            l_rep_type.append(self.REP_TYPE_HERV)
        if i_rep_type & 16 != 0:
            l_rep_type.append(self.REP_TYPE_MIT)
        if i_rep_type & 32 != 0:
            l_rep_type.append(self.REP_TYPE_MSTA)
        return l_rep_type

    def construct_interval_tree(self, sf_rmsk, i_min_copy_len, b_with_chr):
        xannotation = XAnnotation(sf_rmsk)
        #b_with_chr = True
        xannotation.set_with_chr(b_with_chr)
        xannotation.load_rmsk_annotation_with_extnd_div_with_lenth_cutoff(self._boundary_extnd, i_min_copy_len)
        xannotation.index_rmsk_annotation_interval_tree()
        return xannotation

    #iflag: type of repeats will work on
    # def run_pos_filter(self, sf_xtea_rslt, i_rep_type, sf_new_out):#run post filtering for each type
    #     l_rep_type=self.parse_rep_type(i_rep_type)
    #     return
    def get_min_copy_len_by_type(self, i_rep_type):
        l_rep_type = self.parse_rep_type(i_rep_type)
        l_min_copy_len=[]
        for s_type in l_rep_type:
            if s_type is self.REP_TYPE_ALU:
                l_min_copy_len.append(self.REP_ALU_MIN_LEN)
            elif s_type is self.REP_TYPE_L1:
                l_min_copy_len.append(self.REP_LINE_MIN_LEN)
            elif s_type is self.REP_TYPE_SVA:
                l_min_copy_len.append(self.REP_SVA_MIN_LEN)
            elif s_type is self.REP_TYPE_HERV:
                l_min_copy_len.append(self.REP_HERV_MIN_LEN)
            elif s_type is self.REP_TYPE_MIT:
                l_min_copy_len.append(self.REP_MIT_MIN_LEN)
            elif s_type is self.REP_TYPE_MSTA:
                l_min_copy_len.append(self.REP_MSTA_MIN_LEN)
        return l_min_copy_len

####
class XTEARsltParser():
    def __init__(self):
        self._two_side = "two_side"
        self._one_side = "one_side"
        self._one_half_side = "one_half"
        self._other="other"

    def set_rep_support_type(self, s_two, s_one_half, s_one, s_other):
        self._two_side = s_two
        self._one_side = s_one
        self._one_half_side = s_one_half
        self._other=s_other


    def load_in_xTEA_rslt(self, sf_rslt):
        l_rcd=[]
        with open(sf_rslt) as fin_in:
            for line in fin_in:
                fields = line.split()
                tmp_rcd=[]
                for s_tmp in fields:
                    tmp_rcd.append(s_tmp)
                tmp_rcd.append(line)#add the last line into list
                l_rcd.append(tmp_rcd)
        return l_rcd
####

    ####
    def get_ins_sub_type(self, rcd):
        s_type = rcd[32]
        if ("two_side" in s_type) or ("both-side" in s_type) or ("one_side_and_half_transduction" is s_type):
            return self._two_side
        elif "one_side" in s_type:#
            return self._one_side
        elif ("one_half" in s_type) or ("one-half" in s_type):
            return self._one_half_side
        else:
            return self._other

    ####return both-end/both-side or return one-side
    def get_ins_sub_type_alu(self, rcd):
        n_ef_lclip = int(rcd[5])
        n_ef_rclip = int(rcd[6])
        s_type = rcd[32]
        if ("two_side" in s_type) or ("both-side" in s_type) or ("one_side_and_half_transduction" is s_type):
            return self._two_side
        elif "one_side" in s_type:#
            return self._one_side
        elif ("one_half" in s_type) or ("one-half" in s_type):
            if n_ef_lclip>0 and n_ef_rclip>0:
                return self._two_side
            return self._one_half_side
        else:
            return self._other

####
    def get_ins_sub_type_sva(self, rcd):
        s_type = rcd[32]
        if ("two_side" in s_type) or ("both-side" in s_type) or ("one_side_and_half_transduction" is s_type):
            return self._two_side
        elif "one_side" in s_type:#
            return self._one_side
        elif ("one_half" in s_type) or ("one-half" in s_type):
            return self._one_half_side
        else:
            return self._other

####
class XTPRTFilter():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1] != "/":
            self.swfolder += "/"
        self.n_jobs = n_jobs
########Hard code here!!!!!!!!!!!!!!!!!!!
        self.f_side_polyA_cutoff=global_values.ONE_SIDE_POLYA_CUTOFF
        ####

####n_clip is clip reads cutoff
    ####both side are polyA, even though are two-side-tprt-both
    ####Will NOT filter out those left and right clipped reads form different cluster
    def is_polyA_dominant_two_side(self, rcd):
        f_cutoff = self.f_side_polyA_cutoff
        n_lpolyA = int(rcd[9])
        n_rpolyA = int(rcd[10])
        n_ef_lclip=int(rcd[5])
        n_ef_rclip=int(rcd[6])
        s_type = rcd[32]
        #n_ef_clip = n_ef_lclip + n_ef_rclip
        if n_ef_lclip <= 0 or n_ef_rclip<=0:
            return False

        if self._is_two_clip_form_different_cluster(rcd)==True:
            return False

        if ("two_side_tprt_both" in s_type) and (self._is_two_disc_form_different_cluster(rcd)==True):
            return False

        b_lpolyA= ((float(n_lpolyA) / float(n_ef_lclip)) > f_cutoff)
        b_rpolyA = ((float(n_rpolyA) / float(n_ef_rclip)) > f_cutoff)
        if b_lpolyA and b_rpolyA:
            return True
        return False


    ####n_clip is clip reads cutoff
    ####both side are polyA, even though are two-side-tprt-both
    def is_polyA_dominant_two_side_sva(self, rcd, i_pos_head):
        f_cutoff = self.f_side_polyA_cutoff
        n_lpolyA = int(rcd[9])
        n_rpolyA = int(rcd[10])
        n_ef_lclip = int(rcd[5])
        n_ef_rclip = int(rcd[6])
        s_type = rcd[32]
        # n_ef_clip = n_ef_lclip + n_ef_rclip
        if n_ef_lclip <= 0 or n_ef_rclip <= 0:
            return False

        if ("two_side_tprt_both" in s_type) and (self._disc_cluster_hit_cns_head(rcd, i_pos_head)==True):
            return False

        b_lpolyA = ((float(n_lpolyA) / float(n_ef_lclip)) > f_cutoff)
        b_rpolyA = ((float(n_rpolyA) / float(n_ef_rclip)) > f_cutoff)
        if b_lpolyA and b_rpolyA:
            return True
        return False

    ####two side clip, and both hit the end of the tail (polyA started region)
    def is_two_side_clip_both_polyA_sva(self, rcd, i_cns_tail):
        s_lclip_cluster = rcd[19]
        s_rclip_cluster = rcd[20]
        if ("-1" in s_lclip_cluster) or ("-1" in s_rclip_cluster):
            return False
        ll_fields = s_lclip_cluster.split(":")
        i_lstart = int(ll_fields[0])
        i_lend = int(ll_fields[1])
        lr_fields = s_rclip_cluster.split(":")
        i_rstart = int(lr_fields[0])
        i_rend = int(lr_fields[1])
        if i_lstart>=i_cns_tail and i_rstart>=i_cns_tail:
            return True
        return False

    ####
    def is_two_side_clip_both_non_polyA_sva(self, rcd, i_cns_tail):
        s_lclip_cluster = rcd[19]
        s_rclip_cluster = rcd[20]
        if ("-1" in s_lclip_cluster) or ("-1" in s_rclip_cluster):
            return False
        ll_fields = s_lclip_cluster.split(":")
        i_lstart = int(ll_fields[0])
        i_lend = int(ll_fields[1])
        lr_fields = s_rclip_cluster.split(":")
        i_rstart = int(lr_fields[0])
        i_rend = int(lr_fields[1])
        if i_lstart<i_cns_tail and i_rstart<i_cns_tail:
            return True
        return False

    ####if both left and right clipped reads form effective clusters
    ####two clusters are quite far from each other, then return True
    ####If one side is -1:-1, also return True
    def _is_two_clip_form_different_cluster(self, rcd):
        s_lclip_cluster=rcd[19]
        s_rclip_cluster=rcd[20]
        if ("-1" in s_lclip_cluster) or ("-1" in s_rclip_cluster):
            return True
        ll_fields=s_lclip_cluster.split(":")
        i_lstart=int(ll_fields[0])
        i_lend=int(ll_fields[1])
        lr_fields=s_rclip_cluster.split(":")
        i_rstart=int(lr_fields[0])
        i_rend=int(lr_fields[1])
        if abs(i_lend-i_rstart)<global_values.TWO_CLIP_CLUSTER_DIFF_CUTOFF \
                or abs(i_rend-i_lstart)<global_values.TWO_CLIP_CLUSTER_DIFF_CUTOFF:
            return False
        return True

    def _is_two_disc_form_different_cluster(self, rcd):
        s_ldisc_cluster = rcd[21]
        s_rdisc_cluster = rcd[22]
        if ("-1" in s_ldisc_cluster) or ("-1" in s_rdisc_cluster):
            return True
        ll_fields = s_ldisc_cluster.split(":")
        i_lstart = int(ll_fields[0])
        i_lend = int(ll_fields[1])
        lr_fields = s_rdisc_cluster.split(":")
        i_rstart = int(lr_fields[0])
        i_rend = int(lr_fields[1])
        if abs(i_lend - i_rstart) < global_values.TWO_CLIP_CLUSTER_DIFF_CUTOFF \
                or abs(i_rend - i_lstart) < global_values.TWO_CLIP_CLUSTER_DIFF_CUTOFF:
            return False
        return True

    def _disc_cluster_hit_cns_head(self, rcd, i_pos_head):
        s_ldisc_cluster = rcd[21]
        s_rdisc_cluster = rcd[22]
        if ("-1" in s_ldisc_cluster) or ("-1" in s_rdisc_cluster):
            return False
        ll_fields = s_ldisc_cluster.split(":")
        i_lstart = int(ll_fields[0])
        i_lend = int(ll_fields[1])
        lr_fields = s_rdisc_cluster.split(":")
        i_rstart = int(lr_fields[0])
        i_rend = int(lr_fields[1])
        if i_lstart<i_pos_head or i_rstart<i_pos_head:
            return True
        return False

    def has_polyA_signal(self, rcd):
        n_lpolyA = int(rcd[9])
        n_rpolyA = int(rcd[10])
        if n_lpolyA<=0 and n_rpolyA<=0:
            return False
        return True
####
    ####if one side is purly polyA, and the other side has raw clipped reads, but none aligned
    ####then this is viewed as a polyA dominant one
    def is_polyA_dominant_one_side(self, rcd, nclip_half_cutoff):
        f_cutoff = self.f_side_polyA_cutoff
        n_lpolyA = int(rcd[9])
        n_rpolyA = int(rcd[10])
        n_ef_lclip = int(rcd[5])
        n_ef_rclip = int(rcd[6])
        n_lr_clip = int(rcd[35])  # all qualified clipped reads

        if n_ef_lclip<=0 and n_ef_rclip<=0:
            return True

        b_l_no_signal=False
        if n_ef_lclip<=0:#
            b_l_no_signal=True
        b_r_no_signal = False
        if n_ef_rclip <= 0:#
            b_r_no_signal = True
####Potential bug here: if one side has large number of clipped reads, but small number ef clipped reads,
####then, it's still possible the other side have no or little clipped reads

        if b_l_no_signal is True:#left no signal
            b_rpolyA = ((float(n_rpolyA) / float(n_ef_rclip)) > f_cutoff)
            if (b_rpolyA is True) and ((n_lr_clip-n_ef_rclip)>=nclip_half_cutoff):
                return True
        if b_r_no_signal is True:
            b_lpolyA = ((float(n_lpolyA) / float(n_ef_lclip)) > f_cutoff)
            if (b_lpolyA is True) and ((n_lr_clip-n_ef_lclip)>=nclip_half_cutoff):
                return True
        return False

    def is_two_side_tprt_and_with_polyA(self, rcd, i_cns_end):####
        s_type = rcd[32]
        if "two_side_tprt" in s_type:
            s_lclip_cluster = rcd[19]
            s_rclip_cluster = rcd[20]
            if ("-1" in s_lclip_cluster) or ("-1" in s_rclip_cluster):
                return False
            ll_fields = s_lclip_cluster.split(":")
            i_lstart = int(ll_fields[0])
            i_lend = int(ll_fields[1])
            lr_fields = s_rclip_cluster.split(":")
            i_rstart = int(lr_fields[0])
            i_rend = int(lr_fields[1])
            if i_lend>i_cns_end or i_rend>i_cns_end:
                return True
        return False

    ####
    def is_two_side_clip_and_both_hit_end(self, rcd, i_cns_end):
        s_lclip_cluster = rcd[19]
        s_rclip_cluster = rcd[20]
        if ("-1" in s_lclip_cluster) or ("-1" in s_rclip_cluster):
            return False
        ll_fields = s_lclip_cluster.split(":")
        i_lstart = int(ll_fields[0])
        i_lend = int(ll_fields[1])
        lr_fields = s_rclip_cluster.split(":")
        i_rstart = int(lr_fields[0])
        i_rend = int(lr_fields[1])
        if i_lstart>i_cns_end and i_rstart>i_cns_end:
            return True
        return False

    def hit_L1_tail(self, rcd, i_cns_end):
        s_lclip_cluster = rcd[19]
        s_rclip_cluster = rcd[20]
        if ("-1" in s_lclip_cluster) or ("-1" in s_rclip_cluster):
            return True
        ll_fields = s_lclip_cluster.split(":")
        i_lstart = int(ll_fields[0])
        i_lend = int(ll_fields[1])
        lr_fields = s_rclip_cluster.split(":")
        i_rstart = int(lr_fields[0])
        i_rend = int(lr_fields[1])
        if i_lend>=i_cns_end or i_rend>=i_cns_end:
            return True
        return False
####
    def hit_end_of_cns(self, rcd):
        s_hit_end=rcd[34]
        if global_values.HIT_END_OF_CNS == s_hit_end:
            return True
        else:
            return False
####
    def has_enough_indel_reads(self, rcd, ncutoff):
        n_indel_reads=int(rcd[41])
        if n_indel_reads > ncutoff:
            return True
        else:
            return False

    #coverage is abnormal
    def cov_is_abnormal(self, rcd, f_cutoff):
        f_lcov=float(rcd[11])
        f_rcov=float(rcd[12])
        if f_lcov>f_cutoff or f_rcov>f_cutoff:
            return True
        return False
####
class AFConflictFilter():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1] != "/":
            self.swfolder += "/"
        self.n_jobs = n_jobs

    ####
    def get_rep_type(self):
        l_types = []
        l_types.append(global_values.ONE_SIDE_FLANKING)
        l_types.append(global_values.TWO_SIDE)
        l_types.append(global_values.TWO_SIDE_TPRT_BOTH)
        l_types.append(global_values.TWO_SIDE_TPRT)
        l_types.append(global_values.ONE_HALF_SIDE)
        l_types.append(global_values.ONE_HALF_SIDE_TRPT_BOTH)
        l_types.append(global_values.ONE_HALF_SIDE_TRPT)
        l_types.append(global_values.ONE_HALF_SIDE_POLYA_DOMINANT)
        l_types.append(global_values.ONE_SIDE)
        l_types.append(global_values.ONE_SIDE_COVERAGE_CONFLICT)
        l_types.append(global_values.ONE_SIDE_TRSDCT)
        l_types.append(global_values.ONE_SIDE_WEAK)
        l_types.append(global_values.ONE_SIDE_OTHER)
        l_types.append(global_values.ONE_SIDE_SV)
        l_types.append(global_values.ONE_SIDE_POLYA_DOMINANT)
        l_types.append(global_values.TWO_SIDE_POLYA_DOMINANT)
        l_types.append(global_values.HIGH_COV_ISD)
        l_types.append(global_values.OTHER_TYPE)
        return l_types

####
    ####
    def get_cutoff_by_type(self, l_types):
        m_cutoff = {}
        for s_type in l_types:
            if ("two_side" in s_type) or ("both-side" in s_type) or ("one_side_and_half_transduction" is s_type):
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
            elif "one_side" in s_type:  #
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
            elif ("one_half" in s_type) or ("one-half" in s_type):
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
            else:
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
        return m_cutoff

    def is_qualified_rcd(self, s_line, m_cutoff):
        fields = s_line.split()
        n_lpolyA = int(fields[9])
        n_rpolyA = int(fields[10])

        n_ef_clip = int(fields[5]) + int(fields[6])
        n_ef_disc = int(fields[7]) + int(fields[8])
        n_clip = int(fields[35])
        n_full_map = int(fields[36])
        n_disc = int(fields[39])
        n_concod = int(fields[40])

        f_ef_clip = 0.0
        if n_clip != 0:
            f_ef_clip = float(n_ef_clip) / float(n_clip)
        f_ef_disc = 0.0
        if n_disc != 0:
            f_ef_disc = float(n_ef_disc) / float(n_disc)
        f_clip_full_map = 0.0
        if (n_clip + n_full_map) != 0:
            f_clip_full_map = float(n_clip) / float(n_clip + n_full_map)
        f_disc_concod = 0.0
        if (n_disc + n_concod) != 0:
            f_disc_concod = float(n_disc) / float(n_disc + n_concod)

        s_type_ins = fields[32]
        b_pass = self.is_ins_pass_cutoff(m_cutoff, s_type_ins, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod)
        return b_pass
####

    ####
    def is_ins_pass_cutoff(self, m_cutoff, s_type, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod):
        (f_ef_clip_cutoff, f_ef_disc_cutoff, f_clip_full_cutoff, f_disc_concod_cutoff) = m_cutoff[s_type]
        b_ef_clip = (f_ef_clip > f_ef_clip_cutoff)
        b_ef_disc = f_ef_disc > f_ef_disc_cutoff
        b_clip_full = f_clip_full_map > f_clip_full_cutoff
        b_disc_concd = f_disc_concod > f_disc_concod_cutoff
        b_pass = b_ef_clip and b_ef_disc and b_clip_full and b_disc_concd
        return b_pass

    ####
    def is_qualified_mosaic_rcd(self, s_line, m_cutoff):
        fields = s_line.split()
        n_lpolyA = int(fields[9])
        n_rpolyA = int(fields[10])

        n_ef_clip = int(fields[5]) + int(fields[6])
        n_ef_disc = int(fields[7]) + int(fields[8])
        n_clip = int(fields[35])
        n_full_map = int(fields[36])
        n_disc = int(fields[39])
        n_concod = int(fields[40])

        f_ef_clip = 0.0
        if n_clip != 0:
            f_ef_clip = float(n_ef_clip) / float(n_clip)
        f_ef_disc = 0.0
        if n_disc != 0:
            f_ef_disc = float(n_ef_disc) / float(n_disc)
        f_clip_full_map = 0.0
        if (n_clip + n_full_map) != 0:
            f_clip_full_map = float(n_clip) / float(n_clip + n_full_map)
        f_disc_concod = 0.0
        if (n_disc + n_concod) != 0:
            f_disc_concod = float(n_disc) / float(n_disc + n_concod)

        s_type_ins = fields[32]
        b_pass = self.is_ins_pass_mosaic_cutoff(m_cutoff, s_type_ins, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod)
        #print n_ef_clip, n_ef_disc, n_clip, n_full_map, n_disc, n_concod, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod, b_pass
        return b_pass
####
    ####
    def is_ins_pass_mosaic_cutoff(self, m_cutoff, s_type, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod):
        # (f_ef_clip_cutoff, f_ef_disc_cutoff, f_clip_full_cutoff, f_disc_concod_cutoff) = m_cutoff[s_type]
        # b_ef_clip = (f_ef_clip > f_ef_clip_cutoff)
        # b_ef_disc = f_ef_disc > f_ef_disc_cutoff
        # b_clip_full = f_clip_full_map > f_clip_full_cutoff
        # b_disc_concd = f_disc_concod > f_disc_concod_cutoff
        #b_pass = b_ef_clip and b_ef_disc and b_clip_full and b_disc_concd
        (f_upper_af, f_lower_af) = m_cutoff[s_type]
        b_clip_af_qualified=False
        if f_clip_full_map>=f_lower_af and f_clip_full_map<=f_upper_af:
            b_clip_af_qualified = True
        b_disc_concd = True
        if f_disc_concod>f_upper_af:
            b_disc_concd=False
        b_pass=(b_clip_af_qualified and b_disc_concd)
        return b_pass

    ####
    def filter_by_af_conflict(self, sf_hc, n_clip_cutoff, n_disc_cutoff):
        l_types = self.get_rep_type()
        m_cutoff = self.get_cutoff_by_type(l_types)
        self.calc_ratio(m_cutoff, n_clip_cutoff, n_disc_cutoff, sf_hc)
    ####
    def calc_ratio(self, m_cutoff, n_clip_cutoff, n_disc_cutoff, sf_in):
        sf_out = sf_in + ".after_filter"
        with open(sf_in) as fin_in, open(sf_out, "w") as fout_af_filter:
            n_total = 0
            n_hard_pass = 0
            n_pass = 0
            for line in fin_in:
                fields = line.split()
                n_lpolyA = int(fields[9])
                n_rpolyA = int(fields[10])

                n_ef_clip = int(fields[5]) + int(fields[6])
                n_ef_disc = int(fields[7]) + int(fields[8])
                n_clip = int(fields[35])
                n_full_map = int(fields[36])
                n_disc = int(fields[39])
                n_concod = int(fields[40])
                n_total += 1

                if n_ef_clip < n_clip_cutoff:
                    print "ef_clip", line
                    continue
                if n_ef_disc < n_disc_cutoff:
                    print "ef_disc", line
                    continue
                if n_lpolyA + n_rpolyA < 1:
                    print "no polyA", line
                    continue

                n_hard_pass += 1

                f_ef_clip = 0.0
                if n_clip != 0:
                    f_ef_clip = float(n_ef_clip) / float(n_clip)
                f_ef_disc = 0.0
                if n_disc != 0:
                    f_ef_disc = float(n_ef_disc) / float(n_disc)
                f_clip_full_map = 0.0
                if (n_clip + n_full_map) != 0:
                    f_clip_full_map = float(n_clip) / float(n_clip + n_full_map)
                f_disc_concod = 0.0
                if (n_disc + n_concod) != 0:
                    f_disc_concod = float(n_disc) / float(n_disc + n_concod)

                s_type_ins = fields[32]
                b_pass = self.is_ins_pass_cutoff(m_cutoff, s_type_ins, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod)

                if b_pass is True:
                    n_pass += 1
                    fout_af_filter.write(line)
                else:
                    s = 1
                    print line.rstrip()
            print n_hard_pass, n_pass, n_total
####