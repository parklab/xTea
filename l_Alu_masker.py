from l_rep_masker import *

class AluMasker(LRepMasker):
    def __init__(self, swfolder, n_jobs):
        LRepMasker.__init__(self, swfolder, n_jobs)
        self.i_alu_end = 275 #reach the polyA
        self._Alu = "Alu"
        self._min_alu_end_pos=240#allow some clip
        #self._n_cns_full_map=265#the consensus is fully contained inside the contig

    #Here i_max_clip is maximum allowed TSD
    def parse_Alu_from_algnmt(self, sf_algnmt, sf_ref, i_max_clip, sf_out):
        s_open_fmt = "rb"
        l_succeed = []
        l_tbd=[]
        # try:
        self.set_TSD(sf_ref)
        samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
        for algnmt in samfile.fetch():  # check each alignment
            if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                continue
            if algnmt.is_unmapped == True:  ##unmapped
                continue
            if algnmt.is_duplicate == True:  ##duplciation
                continue

            l_cigar = algnmt.cigar
            i_map_pos = algnmt.reference_start
            i_seq_start=0
            if l_cigar[0][0]==4:
                i_seq_start=l_cigar[0][1]

            i_map_end, n_mapped, i_seq_end = self._get_map_interval(l_cigar, i_map_pos)
            if i_map_pos>self.i_alu_end:#almost purely polyA, skip them
                continue

            i_max_all_clip = 2 * i_max_clip
            b_fully_map, i_lclip, i_rclip = self.is_fully_map(algnmt, i_max_all_clip)
            b_lclip = True  ###not used
            b_rclip = True  ###not used
            ####check the small clip region agains the flanking region for TSD
            b_lclip, b_rclip = self.check_TSD_for_site(algnmt, i_lclip, i_rclip, b_lclip, b_rclip)

            tmp_rcd = (i_map_pos, i_map_end, "+", i_seq_start, i_seq_end)
            if algnmt.is_reverse == True:
                tmp_rcd = (i_map_pos, i_map_end, "-", i_seq_start, i_seq_end)
            rslt_rcd = self.get_ins_information_from_algnmt_ori_seq(algnmt, self._Alu, tmp_rcd, None)
            if b_fully_map == False:
                if self.is_qualfied_tbd_case(n_mapped, i_lclip, i_rclip, i_max_clip)==True:
                    l_tbd.append(rslt_rcd)
                continue
####
            ####TBD whether we add this filter: which require reach the end
            ####require polyA???
            # if i_map_end<self._min_alu_end_pos:
            #     continue
            l_succeed.append(rslt_rcd)
        samfile.close()
        self.clean_TSD()
        # except ValueError:
        #     print sf_algnmt, "is empty"

        self.save_non_transduct_ins_info_to_file(l_succeed, global_values.LRD_MIN_INTERNAL_DEL, sf_out)
        sf_tbd=sf_out+"_to_be_decided_cases.txt"
        self.save_non_transduct_ins_info_to_file(l_tbd, global_values.LRD_MIN_INTERNAL_DEL, sf_tbd)

    ####if reach the end, and also only one side is clipped, then view as qualified
    def is_qualfied_tbd_case(self, n_mapped, i_lclip, i_rclip, i_max_clip):
        if n_mapped>=self._min_alu_end_pos:
            if i_lclip<i_max_clip or i_rclip<i_max_clip:
                return True
        else:
            return False

####
#minimap2 -k7 -w5 --sr --frag=yes -A2 -B4 -O4,12 -E2,1 -r150 -p.5 -N5 -f10000,50000 -n1 -m20 -s30 -g200 -2K50m
        # --heap-sort=yes --secondary=yes --cs -a -t 8

####
#report the original sequence, not the consensus

#TSD
#polyA
#
####