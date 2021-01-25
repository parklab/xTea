from l_local_alignment import *
import global_values
import pysam

####Todo-list: 1. the TSD position is the refined insertion position

class LTSD():
    def __init__(self, sf_ref):
        self.sf_ref=sf_ref
        self.f_fa=None
        self.chr_map={}
        self._set_chr_map()
        self.m_ref_chrms = {}
        self.b_with_chr = False
        self.la=Local_alignment()
        self.f_match_cutoff=global_values.LRD_MIN_TSD_MATCH_RATIO

    def open_ref(self):
        self.f_fa = pysam.FastaFile(self.sf_ref)
        for tmp_chrm in self.f_fa.references:
            self.m_ref_chrms[tmp_chrm] = 1
        if "chr1" in self.m_ref_chrms:
            self.b_with_chr = True

    def close_ref(self):
        self.f_fa.close()

    def get_site_short_flanks(self, ins_chrm, ins_pos):
        i_lstart = ins_pos - global_values.LRD_MAX_TSD_LEN*2
        i_lend = ins_pos
        i_rstart = ins_pos
        i_rend = ins_pos + global_values.LRD_MAX_TSD_LEN*2
        ####
        ref_chrm = self.process_chrm_name(ins_chrm, self.b_with_chr)
        # 1. get the flank regions around the breakpoint
        s_lregion = self.f_fa.fetch(ref_chrm, i_lstart, i_lend)
        s_rregion = self.f_fa.fetch(ref_chrm, i_rstart, i_rend)
        return s_lregion, s_rregion

    ####check the small left/right flank region against the clip sequence
    ####if match, then are TSD candidate
    def check_tsd(self, s_lregion, s_rregion, s_clip):
        s_TSD=""
        b_tsd=False
        i_refined_pos=0

        #2. align against with them
        if self.la.is_seqs_matched(s_lregion, s_clip, self.f_match_cutoff) == True:
            return True, s_TSD, i_refined_pos
        if self.la.is_seqs_matched(s_rregion, s_clip, self.f_match_cutoff) == True:
            return True, s_TSD, i_refined_pos
        #
        #3. if no hit, check the reverse complementary one
        s_clip_rc=self.gnrt_reverse_complementary(s_clip)
        if self.la.is_seqs_matched(s_lregion, s_clip_rc, self.f_match_cutoff) == True:
            return True, s_TSD, i_refined_pos
        if self.la.is_seqs_matched(s_rregion, s_clip_rc, self.f_match_cutoff) == True:
            return True, s_TSD, i_refined_pos
        return b_tsd, s_TSD, i_refined_pos

    def gnrt_reverse_complementary(self, s_seq):
        s_rc=""
        for s in s_seq[::-1]:
            if s not in self.chr_map:
                s_rc+="N"
            else:
                s_rc+=self.chr_map[s]
        return s_rc

####
    def _set_chr_map(self):
        self.chr_map['A']='T'
        self.chr_map['a'] = 'T'
        self.chr_map['T'] = 'A'
        self.chr_map['t'] = 'A'
        self.chr_map['C'] = 'G'
        self.chr_map['c'] = 'G'
        self.chr_map['G'] = 'C'
        self.chr_map['g'] = 'C'
        self.chr_map['N'] = 'N'
        self.chr_map['n'] = 'N'

####
####this code is repeated used in several places, need to optimize
    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "b_with_chr"
    def process_chrm_name(self, chrm, b_with_chr):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":##Here remove the "chr"
            b_chrm_with_chr = True

        # print chrm, self.b_with_chr, b_chrm_with_chr ###############################################################
        if b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm
####
####
    def check_TSD_for_site(self, algnmt, ins_chrm, ins_pos, i_lclip, i_rclip, b_lclip, b_rclip, m_tsd):
        # check TSD matching
        b_lclip1=b_lclip
        b_rclip1=b_rclip
        s_lregion, s_rregion = self.get_site_short_flanks(ins_chrm, int(ins_pos))

        b_tsd = False
        if i_lclip < global_values.LRD_MAX_TSD_LEN and i_lclip>=global_values.LRD_MIN_TSD_LEN:
            s_clip = algnmt.query_sequence[:i_lclip]
            s_lregion1=s_lregion[-2*i_lclip:]
            s_rregion1=s_rregion[:2*i_lclip]

            print(s_lregion1, s_rregion1, s_clip, "left")
            b_tsd, s_TSD, i_refined_pos = self.check_tsd(s_lregion1, s_rregion1, s_clip)
            if b_tsd == True:
                b_lclip1 = False
                s_tmp_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                m_tsd[s_tmp_id] = s_clip
        if b_tsd == False and i_rclip < global_values.LRD_MAX_TSD_LEN and i_rclip>=global_values.LRD_MIN_TSD_LEN:
            s_clip = algnmt.query_sequence[-1 * i_rclip:]
            s_lregion1 = s_lregion[-2 * i_rclip:]
            s_rregion1 = s_rregion[:2 * i_rclip]
            print(s_lregion1, s_rregion1, s_clip, "right")
            b_tsd, s_TSD, i_refined_pos = self.check_tsd(s_lregion1, s_rregion1, s_clip)
            if b_tsd == True:
                b_rclip1 = False
                s_tmp_id = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                m_tsd[s_tmp_id] = s_clip
        return b_lclip1, b_rclip1
####