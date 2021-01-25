##01/07/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
from x_contig import *
from l_local_alignment import *

class LRD_Contig(XTEContig):#
    def __init__(self, s_working_folder, n_jobs):
        XTEContig.__init__(self, s_working_folder, n_jobs)

    ####for long read assembled contigs
    ##call is_qualified_TE_insertion_strict_version_v2
    # for each record in format: (sf_asm, sf_flanks, sf_algnmt, ins_chrm, ins_pos)
    def call_MEI_from_long_read_contig_flank_algnmt(self, l_asm_algnmt, sf_mei):
        m_constructed = {}
        m_new_pos={}
        with open(sf_mei, "a") as fout_mei_seqs:
            for rcd in l_asm_algnmt:
                sf_contig = rcd[0]
                sf_algnmt = rcd[2]
                ins_chrm = rcd[3]
                ins_pos = rcd[4]
                i_flk_len = global_values.LRD_EXTND_LEN
                s_open_fmt = "rb"  # open in bam format

                if os.path.isfile(sf_algnmt) == False:
                    continue

                f_min_map_ratio = 0.5 #
                i_slack = 25 #if smaller than this value, then clip will be ignored
                #this is old version
                # m_seqs = self.get_TE_insertion_seq_lrd(sf_contig, sf_algnmt, i_flk_len, ins_pos, f_clip_cutoff,
                #                                        i_slack, s_open_fmt)
                m_seqs, m_refine_pos = self.get_TE_insertion_seq_lrd2(sf_contig, sf_algnmt, i_flk_len, ins_pos,
                                                                      f_min_map_ratio, i_slack, s_open_fmt)
                s_site_info = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                if len(m_seqs) > 0:
                    ###Note, here use global_values.SEPERATOR to seperate chrm and pos, in case chrm contains "_"
                    i_cnt = 0
                    for s_tmp in m_seqs:
                        s_seq = m_seqs[s_tmp]
                        if s_seq != "" and len(s_seq)>=global_values.LRD_MIN_FOCUS_INS_LEN:
                            i_refine_pos=m_refine_pos[s_tmp]
                            m_new_pos[s_site_info]=i_refine_pos
                            s_site = "{0}{1}{2}{3}{4}".format(ins_chrm, global_values.SEPERATOR, i_refine_pos,
                                                              global_values.SEPERATOR, i_cnt)
                            fout_mei_seqs.write(">" + s_site + "\n")
                            fout_mei_seqs.write(s_seq + "\n")
                            i_cnt += 1
                            m_constructed[s_site_info] = 1
                else:
                    print("{0}:{1} is not constructed!".format(ins_chrm, ins_pos))
        return m_constructed, m_new_pos

    ####for long read assembled contigs
    # for each record in format: (sf_asm, sf_flanks, sf_algnmt, ins_chrm, ins_pos)
    def call_ref_copy_from_contig_flank_algnmt(self, l_asm_algnmt, sf_ref_copy):
        m_constructed = {}
        m_new_pos = {}
        with open(sf_ref_copy, "a") as fout_mei_seqs:
            for rcd in l_asm_algnmt:
                sf_contig = rcd[0]
                sf_algnmt = rcd[2]
                ins_chrm = rcd[3]
                ins_pos = rcd[4]
                i_flk_len = global_values.LRD_EXTND_LEN
                s_open_fmt = "rb"  # open in bam format

                if os.path.isfile(sf_algnmt) == False:
                    continue

                f_min_map_ratio = 0.5  #
                i_slack = 25  # if smaller than this value, then clip will be ignored
                # this is old version
                # m_seqs = self.get_TE_insertion_seq_lrd(sf_contig, sf_algnmt, i_flk_len, ins_pos, f_clip_cutoff,
                #                                        i_slack, s_open_fmt)
                m_seqs, m_refine_pos = self.get_TE_insertion_seq_lrd2(sf_contig, sf_algnmt, i_flk_len, ins_pos,
                                                                      f_min_map_ratio, i_slack, s_open_fmt)
                s_site_info = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
                if len(m_seqs) > 0:
                    ###Note, here use global_values.SEPERATOR to seperate chrm and pos, in case chrm contains "_"
                    i_cnt = 0
                    for s_tmp in m_seqs:
                        s_seq = m_seqs[s_tmp]
                        if s_seq != "" and len(s_seq) >= global_values.LRD_MIN_FOCUS_INS_LEN:
                            i_refine_pos = m_refine_pos[s_tmp]
                            m_new_pos[s_site_info] = i_refine_pos
                            s_site = "{0}{1}{2}{3}{4}".format(ins_chrm, global_values.SEPERATOR, ins_pos,
                                                              global_values.SEPERATOR, i_cnt)
                            fout_mei_seqs.write(">" + s_site + "\n")
                            fout_mei_seqs.write(s_seq + "\n")
                            i_cnt += 1
                            m_constructed[s_site_info] = 1
                else:
                    print("{0}:{1} is not constructed!".format(ins_chrm, ins_pos))
        return m_constructed, m_new_pos
####

    ####
    # Override from the parent class (XTEContig)
    # This is a strict version, that require the left and right aligned flanks should be:
    # most region are aligned, aligned to the same contig, of same orientation
    def get_TE_insertion_seq_lrd(self, sf_contig, sf_algnmt, flank_length, offset_ref, f_clip_cutoff, i_slack,
                                 s_open_fmt="r"):
        m_groups = {}  # each contig one group, if the group contain both "left" and "right" contig, then it's a candidates
        try:
            samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
            for algnmt in samfile.fetch():  # check each alignment, and find "left" and "right" flank
                if algnmt.is_secondary:  # filter out secondary, but keep the supplimentary
                    continue
                if algnmt.is_unmapped == True:  ##unmapped
                    continue
                if algnmt.is_duplicate == True:  ##duplciate
                    continue
                # as  the clipped long reads are from the alignment, so already corrected from rc.
                # Note, we are not sure about the assembler whether it generate a reverse ones
                if algnmt.is_reverse == True:
                    continue

                r_ref = algnmt.reference_name
                if r_ref not in m_groups:
                    m_groups[r_ref] = {}

                if (algnmt.query_name in m_groups[r_ref]) and (algnmt.is_supplementary == True):
                    continue
                m_groups[r_ref][algnmt.query_name] = algnmt
        except ValueError:
            print(sf_algnmt, "is empty")

        f_fa = self._open_contigs_file(sf_contig)
        if f_fa is None:
            return {}
        # now check those contigs have both left and right contig aligned
        m_tei = {}  # save the constructed insertion
        for s_contig in m_groups:
            s_left = global_values.LEFT_FLANK
            s_right = global_values.RIGHT_FLANK
            # both left and right contig are aligned to same contig
            if (s_left in m_groups[s_contig]) and (s_right in m_groups[s_contig]):
                # check whether they are of the same orientation
                l_algnmt = m_groups[s_contig][s_left]
                r_algnmt = m_groups[s_contig][s_right]

                # here we don't check whether they are clipped or not
                # we only check whether there are enough mapped region
                b_qualfied, istart, iend = self._get_MEI_seq_pos(l_algnmt, r_algnmt, flank_length, f_clip_cutoff,
                                                                 i_slack)
                if b_qualfied == False:
                    continue

                ##get the inserted part, return it to save into a file, then align all of them togethor.
                s_seq = self._get_inserted_part_lrd(f_fa, s_contig, istart, iend)
                m_tei[s_contig] = s_seq
        self._close_fa(f_fa)
        return m_tei

    ####
    ####
    # seqs for the long reads are from the alignment to the reference, so there is no reverse complementary for the flank
    def _get_MEI_seq_pos(self, l_algnmt, r_algnmt, flank_length, f_clip_cutoff, i_slack):
        b_left = True
        b_l_qlfd, i_l_map = self._is_qualified_algnmt(l_algnmt, b_left, flank_length, f_clip_cutoff)
        if b_l_qlfd == False:
            return False, None, None
        b_left = False
        b_r_qlfd, i_r_map = self._is_qualified_algnmt(r_algnmt, b_left, flank_length, f_clip_cutoff)
        if b_r_qlfd == False:
            return False, None, None

        istart = l_algnmt.reference_start + i_l_map
        iend = r_algnmt.reference_start
        if istart > iend:
            print("Error: start position {0} is larger than end position {1}".format(istart, iend))
            return False, None, None

        if (iend - istart) < i_slack:  ##left and right are concatenate
            return False, None, None

        return True, istart, iend

####
    ####
    # override from the parent class (XTEContig)
    # This is a strict version, that require the left and right aligned flanks should be:
    # most region are aligned, aligned to the same contig, of same orientation
    # i_slack: if the clip part is smaller than this value, then view it as full map
    #
    def get_TE_insertion_seq_lrd2(self, sf_contig, sf_algnmt, flank_length, raw_ins_pos, f_match_cutoff, i_slack,
                                 s_open_fmt="r"):
        m_groups = {}# each contig one group, if the group contain both "left" and "right" contig, then it's a candidates
        try:
            samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
            for algnmt in samfile.fetch():  # check each alignment, and find "left" and "right" flank
                if algnmt.is_secondary or algnmt.is_supplementary:# filter out secondary
                    continue
                if algnmt.is_unmapped == True:##unmapped
                    continue
                if algnmt.is_duplicate == True:##duplciate
                    continue
                # as  the clipped long reads are from the alignment, so already corrected from rc.
                # Note, we are not sure about the assembler whether it generate a reverse ones
                if algnmt.is_reverse == True:
                    continue

                r_ref = algnmt.reference_name
                if r_ref not in m_groups:
                    m_groups[r_ref] = {}

                m_groups[r_ref][algnmt.query_name] = algnmt
            samfile.close()
        except ValueError:
            print(sf_algnmt, "is empty")

        f_fa = self._open_contigs_file(sf_contig)
        if f_fa is None:
            return {}, {}
        # now check those contigs have both left and right contig aligned
        m_tei = {}  # save the constructed insertion
        la=Local_alignment()
        #i_slack=int(flank_length*f_max_clip_ratio)
        m_refined_pos={}
        for s_contig in m_groups:
            s_left = global_values.LEFT_FLANK
            s_right = global_values.RIGHT_FLANK
            # both left and right contig are aligned to same contig
            if (s_left in m_groups[s_contig]) and (s_right in m_groups[s_contig]):
                # check whether they are of the same orientation
                l_algnmt = m_groups[s_contig][s_left]
                r_algnmt = m_groups[s_contig][s_right]

                if len(l_algnmt.cigar)==0 or len(r_algnmt.cigar)==0:
                    continue
                # here we don't check whether they are clipped or not
                # we only check whether there are enough mapped region
                rslt_rcd = self._call_from_two_flank_algnmt_on_one_contig(l_algnmt, r_algnmt, raw_ins_pos, i_slack)
                b_qlfd, refined_ins_pos, istart, iend, s_flank_ck_seq, s_contig_seq_pos=rslt_rcd

                if b_qlfd == False:
                    continue

                if s_contig_seq_pos is not None:#
                    #check the flank_ck_seq and contig_seq match or not
                    i_ctg_start, i_cnt_end=s_contig_seq_pos
                    s_contig_seq=self._get_inserted_part_lrd(f_fa, s_contig, i_ctg_start, i_cnt_end)
                    #if s_contig_seq and s_flank_ck_seq are matched, then use the refined position, otherwise, use the old
                    if len(s_contig_seq)==0 or len(s_flank_ck_seq)==0:
                        refined_ins_pos = raw_ins_pos
                    elif la.is_seqs_matched(s_contig_seq, s_flank_ck_seq, f_match_cutoff)==False:
                        refined_ins_pos=raw_ins_pos
####
                ##get the inserted part, return it to save into a file, then align all of them together
                s_seq = self._get_inserted_part_lrd(f_fa, s_contig, istart, iend)
                m_tei[s_contig] = s_seq
                m_refined_pos[s_contig]=refined_ins_pos
        self._close_fa(f_fa)
        return m_tei, m_refined_pos
####
####
    ####this is the case for two flank regions aligned to the same contig
    ###for left/right flank alignment, we require at most 1 flank has large clip-seq
    ###the other may has < i_clip_slack clip
    ###we assume left/right flank is aligned at the left/right
    ###Note: raw_ins_pos is the original called out insertion position, which maybe not the exact position
    ###Also, in many case, the whole insertion will be called as a "deletion" within the flank-2-contig alignment
    def _call_from_two_flank_algnmt_on_one_contig(self, l_algnmt, r_algnmt, raw_ins_pos, i_clip_slack):
        refined_ins_pos = raw_ins_pos
        istart = -1  # insertion start position on the contig
        iend = -1  # insertion end position on the contig
        t_contig_target_pos = None
        s_flank_ck_seq = ""

        l_lcigar = l_algnmt.cigar
        l_rcigar = r_algnmt.cigar
        #flank length
        i_flk_len=len(l_algnmt.query_sequence)
        i_max_clip=int(i_flk_len/2) #at most half is clipped out

        i_lclip_total_len=0#left flank total clip length
        if l_lcigar[0][0]==5 or l_lcigar[0][0]==4:
            i_lclip_total_len+=l_lcigar[0][1]
        if l_lcigar[-1][0]==5 or l_lcigar[-1][0]==4:
            i_lclip_total_len+=l_lcigar[-1][1]
        i_rclip_total_len=0#right flank total clip length
        if l_rcigar[0][0]==5 or l_rcigar[0][0]==4:
            i_rclip_total_len+=l_rcigar[0][1]
        if l_rcigar[-1][0]==5 or l_rcigar[-1][0]==4:
            i_rclip_total_len+=l_rcigar[-1][1]
####
        b_left_flank_clip = False
        if l_lcigar[-1][0] == 4 and l_lcigar[-1][1] > i_clip_slack:
            b_left_flank_clip = True

        b_right_flank_clip = False
        if l_rcigar[0][0] == 4 and l_rcigar[0][1] > i_clip_slack:
            b_right_flank_clip = True

        b_qlfd = True
        # if both side clipped, then this is not qualified
        if (b_left_flank_clip is True) and (b_right_flank_clip is True):
            b_qlfd = False
        #if one side clipped length is quite large (more than half), then this is also viewed as not qualified
        elif (i_lclip_total_len>i_max_clip) or (i_rclip_total_len>i_max_clip):
            b_qlfd=False
        ####

        elif (b_left_flank_clip is False) and (b_right_flank_clip is False):  # both flank regions are well aligned
            # no need to refine the original insertion position
            l_flk_pos = l_algnmt.reference_start
            istart = self._get_right_side_pos(l_flk_pos, l_lcigar)  #
            iend = r_algnmt.reference_start

        elif (b_left_flank_clip is True) and (b_right_flank_clip is False):
            # check the clipped part aligned to the right side of the contig
            l_flk_pos = l_algnmt.reference_start
            i_lflk_clip_len = l_lcigar[-1][1]
            istart = self._get_right_side_pos(l_flk_pos, l_lcigar)  #
            #print raw_ins_pos, r_algnmt.reference_start, i_lflk_clip_len ###############################################
            iend = r_algnmt.reference_start - i_lflk_clip_len####
            refined_ins_pos -= i_lflk_clip_len
            # here need to check the sequence alignment
            s_lflk = l_algnmt.query_sequence
            s_lflk_clip_seq = s_lflk[-1 * i_lflk_clip_len:]
            s_flank_ck_seq = s_lflk_clip_seq
            t_contig_target_pos = (iend, iend + i_lflk_clip_len)

        elif (b_left_flank_clip is False) and (b_right_flank_clip is True):
            r_flk_clip_len = l_rcigar[0][1]
            refined_ins_pos += r_flk_clip_len
            iend = r_algnmt.reference_start
            s_rflk = r_algnmt.query_sequence
            r_flk_clip_seq = s_rflk[:r_flk_clip_len]
            s_flank_ck_seq = r_flk_clip_seq

            l_flk_pos = l_algnmt.reference_start
            i_lflk_rside = self._get_right_side_pos(l_flk_pos, l_lcigar)
            t_contig_target_pos = (i_lflk_rside, i_lflk_rside + r_flk_clip_len)
            istart = i_lflk_rside + r_flk_clip_len

        if istart>iend or istart<0 or iend<0:
            b_qlfd=False

        return b_qlfd, refined_ins_pos, istart, iend, s_flank_ck_seq, t_contig_target_pos
####
####
    # the two flank regions aligned to two different contigs
    # here have more strict checking:
    # 1. clipped part check (because of inexact insertion position
    # 2. TSD check
    # 3. polyA check
    def _call_from_two_flank_algnmt_on_two_contigs(self, l_algnmt, r_algnmt, raw_ins_pos, i_clip_slack):
        return

    ####
    def _get_right_side_pos(self, map_pos, l_cigar):
        ref_pos = map_pos
        n_mapped = 0
        for (opn, lth) in l_cigar:
            if opn == 0:  # alignment match (M)
                ref_pos += lth
                n_mapped += lth  # accumulate the mapped length
            elif opn == 2 or opn == 3:  # deletion (D) or skipped region (N)
                ref_pos += lth
            elif opn == 7 or opn == 8:  # sequence match (=) or sequence mismatch (X)
                n_mapped += lth
                ref_pos += lth
        return ref_pos

    ####
    def _is_qualified_algnmt(self, algnmt, b_left, flank_length, f_clip_cutoff):
        b_qulified = True
        l_cigar = algnmt.cigar
        f_map_ratio, cnt_map = self.cal_map_ratio(l_cigar, flank_length)  # mapped ratio over all sequences

        clip_len = 0
        if b_left == True:  # left flank
            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:
                clip_len = l_cigar[-1][1]
        else:
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:
                clip_len = l_cigar[0][1]
        if float(clip_len) / float(flank_length) > f_clip_cutoff:
            b_qulified = False
        return b_qulified, cnt_map

    def _open_contigs_file(self, sf_contig):
        f_fa = None
        try:
            pysam.faidx(sf_contig)  # index the fasta file
            f_fa = pysam.FastaFile(sf_contig)
        except pysam.SamtoolsError:
            print('{0} cannot be indexed by samtools faidx\n'.format(sf_contig))
        except ValueError:
            print('Cannot open {0}\n'.format(sf_contig))
        return f_fa

    # parse out the inserted part, from the flank region alignments
    # left_pos and right_pos are the map position on the contigs
    def _get_inserted_part_lrd(self, f_fa, s_contig, istart, iend):
        s_seq1 = ""
        try:
            s_seq1 = f_fa.fetch(s_contig, istart, iend)
        except IndexError:
            print('{0}:{1}-{2} the coordinates are out of range\n'.format(s_contig, istart, iend))
        except ValueError:
            print('{0}:{1}-{2} is invalid region\n'.format(s_contig, istart, iend))
        except:
            print('{0} cannot be opened or retrived correctly\n'.format(s_contig, istart, iend))####
        return s_seq1

    def _close_fa(self, f_fa):
        if f_fa is not None:
            f_fa.close()

    #load in contig (in fasta) to dict
    def _load_in_fasta(self, sf_contig):
        m_contig={}
        with open(sf_contig) as fin_contig:
            s_head=""
            s_seq=""
            for line in fin_contig:
                if len(line)<=0:
                    continue
                if line[0]==">":
                    if "" != s_head:
                        m_contig[s_head]=s_seq
                        s_seq=""
                    s_tmp_fields=line.rstrip()[1:].split()
                    s_head=s_tmp_fields[0]
                else:
                    s_seq+=(line.rstrip())
            #save the last one
            if "" != s_head:
                m_contig[s_head] = s_seq
        return m_contig


    def load_fasta_to_dict(self, sf_contig):
        m_contig=self._load_in_fasta(sf_contig)
        return m_contig

####
    ####this is not used
    def _load_in_fasta2(self, sf_contig):
        m_contig={}
        with pysam.FastxFile(sf_contig) as fh:
            for entry in fh:
                s_id=str(entry.name)
                s_seq=str(entry.sequence)
                m_contig[s_id]=s_seq
        return m_contig
####
####