from l_rep_masker import *

class HERVMasker(LRepMasker):
    def __init__(self, swfolder, n_jobs):
        LRepMasker.__init__(self, swfolder, n_jobs)
        self._herv = "HERV"

    #Here i_max_clip is maximum allowed TSD
    def parse_HERV_from_algnmt(self, sf_algnmt, sf_ref, i_max_clip, sf_out):
        s_open_fmt = "rb"
        l_succeed = []

        self.set_TSD(sf_ref)
        samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
        for algnmt in samfile.fetch():  # check each alignment
            if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                continue
            if algnmt.is_unmapped == True:  ##unmapped
                continue
            if algnmt.is_duplicate == True:  ##duplciation
                continue

            i_max_all_clip = 2 * i_max_clip
            b_fully_map, i_lclip, i_rclip = self.is_fully_map(algnmt, i_max_all_clip)
            b_lclip = True  ###not used
            b_rclip = True  ###not used
            b_lclip, b_rclip = self.check_TSD_for_site(algnmt, i_lclip, i_rclip, b_lclip, b_rclip)
            if b_fully_map==False:
                continue

            l_cigar = algnmt.cigar
            i_map_pos = algnmt.reference_start
            i_seq_start=0
            if l_cigar[0][0]==4:
                i_seq_start=l_cigar[0][1]

            i_map_end, n_mapped, i_seq_end = self._get_map_interval(l_cigar, i_map_pos)

            tmp_rcd = (i_map_pos, i_map_end, "+", i_seq_start, i_seq_end)
            if algnmt.is_reverse == True:
                tmp_rcd = (i_map_pos, i_map_end, "-", i_seq_start, i_seq_end)
            rslt_rcd = self.get_ins_information_from_algnmt_ori_seq(algnmt, self._herv, tmp_rcd, None)
            l_succeed.append(rslt_rcd)
        samfile.close()
        self.clean_TSD()

        self.save_non_transduct_ins_info_to_file(l_succeed, global_values.LRD_MIN_INTERNAL_DEL, sf_out)
####
####