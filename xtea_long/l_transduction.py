##11/17/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu, chong.simon.chu@gmail.com


from x_contig import *
from x_polyA import *
from x_reference import *
from cmd_runner import *
from l_output_fmt_parser import *
from l_TSD import *
from l_rep_masker import *

####Todo-list: 1. if no hit, but most part is aligned, then still save as aligned
####2. source is not saved
####3. TSD part
####4. only works for two-stage transduction only, how about more stages?

####
class LTransduction():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        self.n_jobs = n_jobs
        self.cmd_runner = CMD_RUNNER()
        self.lbinfo = LRegionNameBasicInfo()
        self._3mer = self.lbinfo.get_s_3mer()
        self._5mer = self.lbinfo.get_s_5mer()
        self._both_transdct = self.lbinfo.get_s_both_side()
        self.sf_ref="" #this is set at "call_transduction_from_candidates" function
        self._min_flank_len=45
####
    ####
    ####each record in l_transduct in format:
    # (t_segmt_sup, t_segmt_pri, algnmt, s_tmp_type, i_5mer_chk_len, i_3mer_chk_len)
    # where s_tmp_type is: "3mer", "5mer", or "both"
    # Two ways to search for transduction events:
    # 1. clearly >2 polyA signal
    # 2. the segments aligned to flanking regions (both reference and polymorphic full length copies)
    ####Here, require the full copy has length > i_min_copy_len
    ####i_chk_win is the maximum length for find the polyA signal
    def call_transduction_from_candidates(self, l_transduct, l_succeed, i_min_copy_len, i_chk_win, i_flank_len,
                                          sf_ref, sf_ref_fl_flank):
        # 1. check polyA for each candidates
        l_contain_polyA = self.call_transduction_from_polyA_signal(l_transduct, i_chk_win)

        # 2. find polymorphic full length copies
        ##2.1 find the non-transduction full length copies
        l_slct_fl_non_transdct = self.define_full_length_polymorphic_non_transduct_copy(l_succeed, i_min_copy_len)

        ##2.2 find the transduction full length copies (the transduction region is not included in the flanking)
        # 2.2.1 also collect the flanking regions for transduction polymorphic full length cases
        #collect the clipped 5' and 3' sgmts
        l_slct_fl_transdct, m_tsdc_polym_sgmt = self.define_full_length_polym_transduct_copy(l_transduct, i_min_copy_len)

        # 3. collect the flanking regions for non-transduction polymorphic cases
        xref = XReference()
        #the reference flanking regions (including the non-polym and polym (transduct) ones)
        l_ref_polymorphic_fl = l_slct_fl_non_transdct + l_slct_fl_transdct
        sf_polym_ref_fl_flank = self.swfolder + "polymorphic_ref_full_length_copy_flanks.fa"
        m_polym_fl_ref_flank=xref.gnrt_flank_regions_of_polymerphic_insertions(
            l_ref_polymorphic_fl, i_flank_len, sf_ref, sf_polym_ref_fl_flank)
####
        # 4. align the transduction regions to the flanking regions
        ##4.1 first merget the flank regions
        sf_merged_flank = self.swfolder + "merged_flanks.fa"
        sf_ghost_flank = self.swfolder + "ghost_rep_copy_flanks.fa"
        if os.path.isfile(sf_ghost_flank)==True:
            cmd = "cat {0} {1} ".format(sf_ref_fl_flank, sf_ghost_flank)
            self.cmd_runner.run_cmd_to_file(cmd, sf_merged_flank)
        else:
            sf_merged_flank=sf_ref_fl_flank
            #cmd = "cat {0} {1}".format(sf_ref_fl_flank, sf_polym_ref_fl_flank)
            #self.cmd_runner.run_cmd_to_file(cmd, sf_merged_flank)
####
        # if os.path.isfile(sf_ghost_flank)==True:
        #     cmd = "cat {0} {1} {2}".format(sf_ref_fl_flank, sf_polym_ref_fl_flank, sf_ghost_flank)
        #     self.cmd_runner.run_cmd_to_file(cmd, sf_merged_flank)
        # else:
        #     cmd = "cat {0} {1}".format(sf_ref_fl_flank, sf_polym_ref_fl_flank)
        #     self.cmd_runner.run_cmd_to_file(cmd, sf_merged_flank)

        # 5. algn and call the transduction event
        # each key in m_candidates with format:{0}_{1}_{2}".format(ins_chrm, ins_pos, self._5mer)
        # each key in m_with_polya with format: "{0}~{1}".format(ins_chrm, ins_pos)
        self.sf_ref=sf_ref
        #each record of m_candidates in format: m[s_id]=(s_source, s_seq)
        m_candidates, m_with_polya, m_transduct_tsd, m_tbd = self.call_transduction_for_sites(
            l_transduct, l_contain_polyA, m_polym_fl_ref_flank, sf_merged_flank, m_tsdc_polym_sgmt)
        # save the transduction information to the file
        #self._slct_save_transduction_to_file(l_transduct, m_candidates, m_with_polya, sf_rslts)
        return m_candidates, m_with_polya, m_transduct_tsd, m_tbd
####
####
    ####Given a list of candidate sites saved in l_transduct,
    ########each record in format:
    ####Algorithms: 1. first collect the clip regions of the candidate transduction seqs
    ####                1.1. also write the polymorphic full length clipped sgmts to a file for alignment
    ####                    all clipped sgmts are saved in m_tsdc_polym_sgmt
    ####            2. algn the clip-seqs to flank regions
    ####            3. call out the candidates
    def call_transduction_for_sites(self, l_transduct, l_contain_polyA, m_polym_fl_ref_flank, sf_flank, m_tsdc_polym_sgmt):
        sf_clip_seq = self.swfolder + "clip_seqs_for_transduction.fa"
        m_with_polya = {}
        m_slct_canddt={}#for temporary use to save the transduction full lenth copy offsprings
        m_slct_tsd={}#for temporary usage too
        with open(sf_clip_seq, "w") as fout_clip:
            for rcd, b_contain_polyA in zip(l_transduct, l_contain_polyA):
                algnmt = rcd[2]
                seq = algnmt.query_sequence  ###this is already converted back to consensus
                #s_transduction_type = rcd[3]
                i_5mer_chk_len = int(rcd[4])
                i_3mer_chk_len = int(rcd[5])

                s_5mer_clip = ""
                s_3mer_clip = ""

                if i_5mer_chk_len>self._min_flank_len:
                    s_5mer_clip = seq[:i_5mer_chk_len]
                if i_3mer_chk_len>self._min_flank_len:
                    s_3mer_clip = seq[-1 * i_3mer_chk_len:]
                (ins_chrm, ins_pos) = self.lbinfo.get_ins_info_from_qname(algnmt.query_name)
                s_id = "{0}~{1}".format(ins_chrm, ins_pos)
                m_with_polya[s_id] = b_contain_polyA

                #check the transduction full length copy first
                m_tmp_canddt, m_tmp_tsd=self.check_transduction_with_other_trsdct_fl_copy_merge_flank(
                    ins_chrm, ins_pos, s_5mer_clip, s_3mer_clip, m_tsdc_polym_sgmt, m_polym_fl_ref_flank)
                for s_tmp_id in m_tmp_tsd:
                    m_slct_tsd[s_tmp_id]=m_tmp_tsd[s_tmp_id]
                if len(m_tmp_canddt)>0:#find a hit from the transduction full length copy flanking region
                    for s_tmp_id in m_tmp_canddt:
                        m_slct_canddt[s_tmp_id]=m_tmp_canddt[s_tmp_id]
                    continue
####
                if s_5mer_clip != "":
                    sinfo = ">{0}~{1}~{2}".format(ins_chrm, ins_pos, self._5mer)
                    fout_clip.write(sinfo + "\n")
                    fout_clip.write(s_5mer_clip + "\n")
                if s_3mer_clip != "":
                    sinfo = ">{0}~{1}~{2}".format(ins_chrm, ins_pos, self._3mer)
                    fout_clip.write(sinfo + "\n")
                    fout_clip.write(s_3mer_clip + "\n")

        # algn the clip sequence to the flank regions
        xtea_contig = XTEContig(self.swfolder, self.n_jobs)
        sf_algnmt = self.swfolder + "transduction_clip_2_flank.sorted.bam"
        xtea_contig.align_short_contigs_2_cns_minimap2(sf_flank, sf_clip_seq, self.n_jobs, sf_algnmt)

        # check transduction from the algnmts
        m_candidates, m_transduct_tsd, m_tbd = self.call_transduction_from_clip_2_flank_algnmt(sf_algnmt)
        #merget the sites from transduction full length to the current call set
        for s_tmp_id in m_slct_canddt:
            m_candidates[s_tmp_id]=m_slct_canddt[s_tmp_id]
        for s_tmp_id in m_slct_tsd:
            m_transduct_tsd[s_tmp_id]=m_slct_tsd[s_tmp_id]
        return m_candidates, m_with_polya, m_transduct_tsd, m_tbd

####
    ####Align the potential transduction part (clipped part) to the reference genome
    ####to call out: 1) potential source of reference orphan transduction
    ####             2) Or, a truncated copy with transduction (same source with the orphan transduction)
    def call_orphan_transduction_for_sites(self, sf_ref, sf_clip_seq, m_fl_ins, sf_output):
        #write the clipped seqs to a file
        #align to ref
        #find the unique mapped ones
        #if the clipped part from a full length L1, then itself is the "source" (or secondary source)
        #otherwise, it is also a child of some unknown source


        return
####
    def check_transduction_with_other_trsdct_fl_copy(self, ins_chrm, ins_pos, s_5mer_clip, s_3mer_clip,
                                                     m_tsdc_polym_sgmt, m_polym_fl_ref_flank):
        # temp file to save the clipped sgmts for the current insertion
        sf_tmp_sgmt = self.swfolder + "tmp_site_clip_sgmts.fa"
        with open(sf_tmp_sgmt, "w") as fout_sgmt:
            if s_5mer_clip != "" and len(s_5mer_clip)>self._min_flank_len:
                sinfo = ">{0}~{1}~{2}".format(ins_chrm, ins_pos, self._5mer)
                fout_sgmt.write(sinfo + "\n")
                fout_sgmt.write(s_5mer_clip + "\n")
            if s_3mer_clip != "" and len(s_3mer_clip)>self._min_flank_len:
                sinfo = ">{0}~{1}~{2}".format(ins_chrm, ins_pos, self._3mer)
                fout_sgmt.write(sinfo + "\n")
                fout_sgmt.write(s_3mer_clip + "\n")

        ## temp file to save all the clipped sgnmts(exclude the current one)
        sf_tmp_all_sgmt = self.swfolder + "tmp_site_all_clip_sgmts.fa"
        with open(sf_tmp_all_sgmt, "w") as fout_all_clip:
            s_cur_id = "{0}{1}{2}".format(ins_chrm, global_values.S_DELIM, ins_pos)
            for s_tmp_id in m_polym_fl_ref_flank:
                if s_tmp_id == s_cur_id:
                    continue
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][0]+"\n")
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][1] + "\n")
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][2] + "\n")
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][3] + "\n")

            for s_tmp_id in m_tsdc_polym_sgmt:
                # s_flank_id=s_tmp_id+"{0}polymorphic{1}"
                if s_tmp_id == s_cur_id:
                    continue
                tmp_fields = s_tmp_id.split(global_values.S_DELIM)
                tmp_chrm = tmp_fields[0]
                tmp_pos = tmp_fields[1]
                sub_family = "polymorphic"
                bi_rc = 0
                s_left_head = ">{0}{1}{2}{3}{4}{5}{6}{7}{8}L".format(tmp_chrm, global_values.S_DELIM, tmp_pos,
                                                                     global_values.S_DELIM, tmp_pos,
                                                                     global_values.S_DELIM,
                                                                     sub_family, global_values.S_DELIM, bi_rc)

                s_lsgmt = m_tsdc_polym_sgmt[s_tmp_id][0]
                if len(s_lsgmt)>self._min_flank_len:
                    fout_all_clip.write(s_left_head + "\n")
                    fout_all_clip.write(s_lsgmt + "\n")
                s_right_head = ">{0}{1}{2}{3}{4}{5}{6}{7}{8}R".format(tmp_chrm, global_values.S_DELIM, tmp_pos,
                                                                      global_values.S_DELIM, tmp_pos,
                                                                      global_values.S_DELIM,
                                                                      sub_family, global_values.S_DELIM, bi_rc)
                s_rsgmt = m_tsdc_polym_sgmt[s_tmp_id][1]
                if len(s_rsgmt)>self._min_flank_len:
                    fout_all_clip.write(s_right_head + "\n")
                    fout_all_clip.write(s_rsgmt + "\n")
        # algn the current clipped to all other "transduction" full length ones
        xtea_contig = XTEContig(self.swfolder, self.n_jobs)
        sf_tmp_algmt = self.swfolder + "transduction_clip_2_all_other_clip.sorted.bam"
        xtea_contig.align_short_contigs_2_cns_minimap2(sf_tmp_all_sgmt, sf_tmp_sgmt, self.n_jobs, sf_tmp_algmt)
        m_candidates, m_transduct_tsd, m_tbd = self.call_transduction_from_clip_2_flank_algnmt(sf_tmp_algmt)
        return m_candidates, m_transduct_tsd
####
    ####
    #in this version, for full length transduction copy, we merge the clipped sgmt with the flanking region
    def check_transduction_with_other_trsdct_fl_copy_merge_flank(self, ins_chrm, ins_pos, s_5mer_clip, s_3mer_clip,
                                                     m_tsdc_polym_sgmt, m_polym_fl_ref_flank):
        # temp file to save the clipped sgmts for the current insertion
        sf_tmp_sgmt = self.swfolder + "tmp_site_clip_sgmts.fa"
        with open(sf_tmp_sgmt, "w") as fout_sgmt:
            if s_5mer_clip != "" and len(s_5mer_clip)>self._min_flank_len:
                sinfo = ">{0}~{1}~{2}".format(ins_chrm, ins_pos, self._5mer)
                fout_sgmt.write(sinfo + "\n")
                fout_sgmt.write(s_5mer_clip + "\n")
            if s_3mer_clip != "" and len(s_3mer_clip)>self._min_flank_len:
                sinfo = ">{0}~{1}~{2}".format(ins_chrm, ins_pos, self._3mer)
                fout_sgmt.write(sinfo + "\n")
                fout_sgmt.write(s_3mer_clip + "\n")

        ## temp file to save all the clipped sgnmts(exclude the current one)
        sf_tmp_all_sgmt = self.swfolder + "tmp_site_all_clip_sgmts.fa"
        with open(sf_tmp_all_sgmt, "w") as fout_all_clip:
            s_cur_id = "{0}{1}{2}".format(ins_chrm, global_values.S_DELIM, ins_pos)
            m_tmp_saved={}
            for s_tmp_id in m_tsdc_polym_sgmt:#
                # s_flank_id=s_tmp_id+"{0}polymorphic{1}"
                if s_tmp_id == s_cur_id:
                    continue

                tmp_fields = s_tmp_id.split(global_values.S_DELIM)
                tmp_chrm = tmp_fields[0]
                tmp_pos = tmp_fields[1]
                sub_family = "polymorphic"
                bi_rc = 0
                s_left_head = ">{0}{1}{2}{3}{4}{5}{6}{7}{8}L".format(tmp_chrm, global_values.S_DELIM, tmp_pos,
                                                                     global_values.S_DELIM, tmp_pos,
                                                                     global_values.S_DELIM,
                                                                     sub_family, global_values.S_DELIM, bi_rc)

                s_lsgmt = m_tsdc_polym_sgmt[s_tmp_id][0]
                if s_tmp_id in m_polym_fl_ref_flank:#concatenate the left flank with the "transduction" region
                    s_lsgmt+=m_polym_fl_ref_flank[s_tmp_id][1]
                if len(s_lsgmt)>self._min_flank_len:
                    fout_all_clip.write(s_left_head + "\n")
                    fout_all_clip.write(s_lsgmt + "\n")
                    m_tmp_saved[s_tmp_id]=1
                s_right_head = ">{0}{1}{2}{3}{4}{5}{6}{7}{8}R".format(tmp_chrm, global_values.S_DELIM, tmp_pos,
                                                                      global_values.S_DELIM, tmp_pos,
                                                                      global_values.S_DELIM,
                                                                      sub_family, global_values.S_DELIM, bi_rc)
                s_rsgmt = m_tsdc_polym_sgmt[s_tmp_id][1]
                if s_tmp_id in m_polym_fl_ref_flank:
                    s_rsgmt=m_polym_fl_ref_flank[s_tmp_id][3]+s_rsgmt
                if len(s_rsgmt)>self._min_flank_len:
                    fout_all_clip.write(s_right_head + "\n")
                    fout_all_clip.write(s_rsgmt + "\n")
                    m_tmp_saved[s_tmp_id] = 1

            ############################################################################################################
            for s_tmp_id in m_polym_fl_ref_flank:
                if s_tmp_id == s_cur_id:
                    continue
                if s_tmp_id in m_tmp_saved:
                    continue
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][0] + "\n")
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][1] + "\n")
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][2] + "\n")
                fout_all_clip.write(m_polym_fl_ref_flank[s_tmp_id][3] + "\n")
        # algn the current clipped to all other "transduction" full length ones
        xtea_contig = XTEContig(self.swfolder, self.n_jobs)
        sf_tmp_algmt = self.swfolder + "transduction_clip_2_all_other_clip.sorted.bam"
        xtea_contig.align_short_contigs_2_cns_minimap2(sf_tmp_all_sgmt, sf_tmp_sgmt, self.n_jobs, sf_tmp_algmt)
        m_candidates, m_transduct_tsd, m_tbd = self.call_transduction_from_clip_2_flank_algnmt(sf_tmp_algmt)

        return m_candidates, m_transduct_tsd
    ####
####
    def _prep_source_id(self, s_flank_id, i_hit_pos, s_rc, n_mapped):
        s_source_id = s_flank_id + global_values.SEPERATOR + str(i_hit_pos) + global_values.SEPERATOR + s_rc
        s_source_id += "{0}{1}".format(global_values.SEPERATOR, n_mapped)
        return s_source_id
####
    ####Check whether 3' and 5' are matched between the source and insertion.
    ####For each candidate, save the following information:
    ####1. source-position: with hit location (check whether there is a second stage hit);
    ####2. transduction sequence, 5' or 3';
    def call_transduction_from_clip_2_flank_algnmt(self, sf_algnmt):
        m_slct_transduction = {}
        l_tsd=LTSD(self.sf_ref)
        l_tsd.open_ref()
        m_transduct_tsd = {}  # for transduction tsd
        m_tbd={}
        #m_unprocessed={} #save for further check, like align to reference to check orphan transduction
        try:####
            s_open_fmt="rb"
            samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
            for algnmt in samfile.fetch():  # check each alignment
                s_id = algnmt.query_name  # in format: "{0}~{1}~{2}".format(ins_chrm, ins_pos, self._5mer/3mer)
                if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                    continue
                if algnmt.is_unmapped == True:  ##unmapped
                    #m_unprocessed[s_id]=algnmt.query_sequence
                    continue
                if algnmt.is_duplicate == True:  ##duplciation
                    continue
                mapq = algnmt.mapping_quality
                if mapq < global_values.MINIMAL_TRANSDUCT_MAPQ:
                    #m_unprocessed[s_id] = algnmt.query_sequence
                    continue
                l_cigar = algnmt.cigar
                if l_cigar[0][0]==5 or l_cigar[-1][0]==5:#skip hard clip
                    #m_unprocessed[s_id] = algnmt.query_sequence
                    continue

                s_id_fields = s_id.split(global_values.SEPERATOR)
                ins_chrm = s_id_fields[0]
                ins_pos = s_id_fields[1]
                b_ins_3prime=True
                if s_id_fields[-1] is self._5mer:
                    b_ins_3prime=False

                s_flank_id = algnmt.reference_name  #in format: chrY~3443551~3449565~L1HS~0L, where 0 means not rc
                b_source_3prime=self._get_source_flank_3_5_prime(s_flank_id)
####comment out temporarily
                # if b_ins_3prime is not b_source_3prime:#source and insertion flanking region direction doesn't match
                #     print s_id, s_flank_id, "source and insertion flank not match!", l_cigar
                #     m_unprocessed[s_id] = algnmt.query_sequence
                #     continue

                s_rc="1"
                if algnmt.is_reverse==False:
                    s_rc="0"
                i_hit_pos=algnmt.reference_start #
                #s_source_id=s_flank_id+global_values.SEPERATOR+str(i_hit_pos)+global_values.SEPERATOR+s_rc

                #1.check whether it is fully mapped
                i_seq_len = len(algnmt.query_sequence)
                if i_seq_len == 0:
                    continue

                i_lclip, i_rclip=self._get_clip_lenth(l_cigar)
                n_mapped=i_seq_len-i_lclip-i_rclip
                #s_source_id+="{0}{1}".format(global_values.SEPERATOR, n_mapped)
                s_source_id=self._prep_source_id(s_flank_id, i_hit_pos, s_rc, n_mapped)
                i_seq_end=i_seq_len-i_rclip
                # if i_rclip==0:
                #     i_seq_end=i_seq_len
                s_mapped_seq=algnmt.query_sequence[i_lclip:i_seq_end]
                #check TSD
                b_lclip=True
                b_rclip=True
                l_tsd.check_TSD_for_site(algnmt, ins_chrm, ins_pos, i_lclip, i_rclip, b_lclip, b_rclip, m_transduct_tsd)
                if b_lclip==False:
                    i_lclip=0
                if b_rclip==False:
                    i_rclip=0

################
                #print s_id, s_mapped_seq, s_source_id

                ####find the primary hit interval
                b_qualified = self._is_qualified_algnmt(i_lclip+i_rclip, i_seq_len)
                if b_qualified == True:
                    if s_id in m_slct_transduction:
                        print("Find more than one transduction source of {0}! Check further!".format(s_id))
                    m_slct_transduction[s_id] = (s_source_id, s_mapped_seq)
####
                else:
                    if algnmt.has_tag("SA") == False:  # no supplementary alignment
                    ####check whether most part is aligned
                        m_tbd[s_id] = (s_source_id, s_mapped_seq) ####need further checking, maybe source flank not included
                    else:  # check the SA field for a second hit
                        rep_msk=LRepMasker(self.swfolder, self.n_jobs)#
                        #each one in format: (s_first_source, map_pos, i_map_end, s_rc, i_seq_start, i_seq_end)
                        t_pri, t_sup=rep_msk.parse_sa_for_second_level_transduct(algnmt)
                        #first check whether there is large overlap
                        pri_start=t_pri[-2]
                        pri_end=t_pri[-1]
                        sa_start=t_sup[-2]
                        sa_end=t_sup[-1]
                        if self._two_segmts_overlap(pri_start, pri_end, sa_start, sa_end)==True:
                            m_tbd[s_id] = (s_source_id, s_mapped_seq)
                            #m_unprocessed[s_id] = algnmt.query_sequence
                            continue
                        ####
                        if float((pri_end-pri_start)+(sa_end-sa_start))/float(i_seq_len) \
                            < global_values.LRD_TRSDCT_MIN_MAP_RATION:
                            m_tbd[s_id] = (s_source_id, s_mapped_seq)
                            #m_unprocessed[s_id] = algnmt.query_sequence
                            continue
                        else:
                            #save the two level transduction
                            s_new_source=s_source_id+":"+t_sup[0]
                            s_new_seq=s_mapped_seq+":"+algnmt.query_sequence[sa_start:sa_end]
                            m_slct_transduction[s_id] = (s_new_source, s_new_seq)
####
                #2. if not, then check TSD
                #3. Then check SA
                #4. in the end check whether most part is aligned
                # make sure most of the part is aligned
            samfile.close()
        except ValueError:
            print(sf_algnmt, "is empty")

        l_tsd.close_ref()
        return m_slct_transduction, m_transduct_tsd, m_tbd

#####check the supplementary alignemnt for a second hit searching
    # def check_second_hit(self, algnmt):
    #     s_tag = algnmt.get_tag("SA")
    #     rcd_tag = self._parse_SA_fields(s_tag)
    #     i_sa_pos = rcd_tag[0]
    #     b_sa_rc = rcd_tag[1]
    #     s_sa_cigar = rcd_tag[2]
    #     i_sa_mapq = rcd_tag[3]  # this will not be checked
    #     l_sa_cigar = self._cvt_scigar_lcigar(s_sa_cigar)
    #     s_sa_rc = "+"
    #     if b_sa_rc == True:
    #         s_sa_rc = "-"
    #     i_sa_lclip, i_sa_rclip=self._get_clip_lenth(l_sa_cigar)

    def _two_segmts_overlap(self, i_seq_start, i_seq_end, i_sa_start, i_sa_end):
        if i_seq_start<i_sa_start and i_seq_end> i_sa_end:
            #print "debug 3"
            return False
        if i_sa_start < i_seq_start and i_sa_end>i_seq_end:
            #print "debug 4"
            return False
        #if overlap, then overlap region is short:
        if i_seq_end>i_sa_start and i_seq_end<i_sa_end:
            if (i_seq_end-i_sa_start)>global_values.LRD_PRI_SUP_MAX_OVRLAP:
                #print "debug 5"
                return False
        if i_sa_end>i_seq_start and i_sa_end<i_seq_end:
            if (i_sa_end-i_seq_start)>global_values.LRD_PRI_SUP_MAX_OVRLAP:
                #print "debug 6"
                return False
        return True
    ####
    def _get_source_flank_3_5_prime(self, s_name):
        fields=s_name.split(global_values.SEPERATOR)
        rc_lr=fields[-1]
        b_rc=False
        b_left=False
        if rc_lr[0]=="1":
            b_rc=True
        if rc_lr[1]=="L":
            b_left=True

        b_3prime=True
        if b_left==True and b_rc==False:
            b_3prime=False
        elif b_left==False and b_rc==True:
            b_3prime=False

        return b_3prime

    ####get clip length
    def _get_clip_lenth(self, l_cigar):
        i_lclip = 0
        i_rclip = 0
        if l_cigar[0][0] == 4:
            i_lclip = l_cigar[0][1]
        if l_cigar[-1][0] == 4:
            i_rclip = l_cigar[-1][1]
        return i_lclip, i_rclip
    ####
    def _is_qualified_algnmt(self, i_clip, i_total):
        if i_total==0:
            return False
        if (float(i_clip) / float(i_total)) > (1 - global_values.LRD_TRSDCT_MIN_MAP_RATION):
            return False
        else:
            return True

####
####This code duplicates the same function located within l_rep_masker.py
    # def _cvt_scigar_lcigar(self, s_cigar):
    #     print s_cigar
    #     m_op = {'S': 4, 'H': 5, 'M': 0, 'D': 2, 'I': 1, 'P': 6, 'X': 8, '=': 7, 'N': 3, 'B': 9}
    #     s_tmp = ""
    #     l_op = []
    #     for ch in s_cigar:
    #         if ch in m_op:
    #             i_op = m_op[ch]
    #             i_len = int(s_tmp)
    #             l_op.append((i_op, i_len))
    #             s_tmp = ''
    #         else:
    #             s_tmp += ch
    #     return l_op

####
####This code duplicates the same function located within l_rep_masker.py
    ####return value in format: "LINE1,6067,+,1382S23M1942S,0,0;"
    ####rname ,pos ,strand ,CIGAR ,mapQ ,NM;
    # def _parse_SA_fields(self, s_tag):
    #     s_tag_fields = s_tag.split(",")
    #     i_sa_pos = int(s_tag_fields[1])
    #     b_rc = True
    #     if s_tag_fields[2] == "+":
    #         b_rc = False
    #     s_sa_cigar = s_tag_fields[3]
    #     i_sa_mapq = int(s_tag_fields[4])
    #
    #     rcd = (i_sa_pos, b_rc, s_sa_cigar, i_sa_mapq)
    #     return rcd
#
####This code duplicates the same function located within l_rep_masker.py
    # get the mapped interval
    # def _get_map_interval(self, l_cigar, i_map_pos):
        # ref_pos = i_map_pos
        # seq_pos = 0
        # n_mapped = 0
        # for (opn, lth) in l_cigar:
        #     if opn == 0:  # alignment match (M)
        #         ref_pos += lth
        #         seq_pos += lth
        #         n_mapped += lth  # accumulate the mapped length
        #     elif opn == 1:  # insertion
        #         seq_pos += lth
        #     elif opn == 2 or opn == 3:  # deletion (D) or skipped region (N)
        #         ref_pos += lth
        #     elif opn == 4:  # soft-clip (S)
        #         seq_pos += lth
        #     elif opn == 5 or opn == 6:  # hard-clip (H) or padding (P)
        #         seq_pos += lth
        #     elif opn == 7 or opn == 8:  # sequence match (=) or sequence mismatch (X)
        #         ref_pos += lth
        #         seq_pos += lth
        #         n_mapped += lth
        # if len(l_cigar)>1 and (l_cigar[-1][0]==4 or l_cigar[-1][0]==5):
        #     seq_pos-=l_cigar[-1][1]
        # # map position
        # return ref_pos, n_mapped, seq_pos

    ####Call out the full length polymorphic insertion, return a lits of candidates
    ####Each record in format:
    # (ins_chrm, ins_pos, rep_type, structure_type, s_region1, s_region2, TSD, transduct_src, s_seq, transduct_seq)
    # s_region in format: "{0}:{1}:{2}".format(i_start, i_end, t_seg1[2])
    def define_full_length_polymorphic_non_transduct_copy(self, l_succeed, i_min_len):
        l_slcted = []
        for rcd in l_succeed:
            s_region1 = rcd[4]
            s_region2 = rcd[5]
            i_reg1_len = 0
            i_reg2_len = 0
            if s_region1 is not "None":
                i_reg1_len = self.lbinfo.parse_region_fields(s_region1)
            if s_region2 is not "None":
                i_reg2_len = self.lbinfo.parse_region_fields(s_region2)
            i_ins_len = i_reg1_len + i_reg2_len
            if i_ins_len >= i_min_len:
                l_slcted.append((rcd[0], int(rcd[1])))
        return l_slcted
    ####

    # Here, if a candidate transduction event has long enough mapped region on the consensus, then
    # will also be collected as potential full length source copy
    ###each record in format: t_segmt_pri, t_segmt_sup, algnmt, s_transduction_type, i_5mer_chk_len, i_3mer_chk_len
    ########t_segmt_pri in format (start, end, direction[+/-])
    def define_full_length_polym_transduct_copy(self, l_transduct, i_min_len):
        l_slcted = []
        m_tsdc_polym_sgmt={}
        for rcd in l_transduct:
            t_segmt_pri = rcd[0]
            t_segmt_sup = rcd[1]
            i_mapped = 0
            if t_segmt_pri is not None:
                i_mapped += (int(t_segmt_pri[1]) - int(t_segmt_pri[0]) + 1)
            if t_segmt_sup is not None:
                i_mapped += (int(t_segmt_sup[1]) - int(t_segmt_sup[0]) + 1)
            ####
            algnmt = rcd[2]
            (ins_chrm, ins_pos) = self.lbinfo.get_ins_info_from_qname(algnmt.query_name)
            if i_mapped >= i_min_len:
                l_slcted.append((ins_chrm, int(ins_pos)))
                s_seq = algnmt.query_sequence
                i_5mer_chk_len = rcd[4]
                i_3mer_chk_len = rcd[5]
                s_5mer_clip = s_seq[:i_5mer_chk_len]
                s_3mer_clip = s_seq[-1 * i_3mer_chk_len:]
                s_id="{0}{1}{2}".format(ins_chrm, global_values.S_DELIM, ins_pos)
                m_tsdc_polym_sgmt[s_id]=(s_5mer_clip, s_3mer_clip)
        return l_slcted, m_tsdc_polym_sgmt

####
    # # Here, if a candidate transduction event has long enough mapped region on the consensus, then
    # # will also be collected as potential full length source copy
    # ###each record in format: t_segmt_pri, t_segmt_sup, algnmt, s_transduction_type, i_5mer_chk_len, i_3mer_chk_len
    # ########t_segmt_pri in format (start, end, direction[+/-])
    # ###also collect the clipped part as partial flanking regions
    # def define_full_length_polymorphic_transduct_with_clip_seq(self, l_transduct, i_min_len):
    #     l_slcted = []
    #     for rcd in l_transduct:
    #         t_segmt_pri = rcd[0]
    #         t_segmt_sup = rcd[1]
    #         i_mapped = 0
    #         if t_segmt_pri is not None:
    #             i_mapped += (int(t_segmt_pri[1]) - int(t_segmt_pri[0]) + 1)
    #         if t_segmt_sup is not None:
    #             i_mapped += (int(t_segmt_sup[1]) - int(t_segmt_sup[0]) + 1)
    #         ####
    #         algnmt = rcd[2]
    #         (ins_chrm, ins_pos) = self.lbinfo.get_ins_info_from_qname(algnmt.query_name)
    #         if i_mapped >= i_min_len:
    #             l_slcted.append((ins_chrm, int(ins_pos)))
    #
    #         ####Here need to define the clipped part
    #         s_seq=algnmt.query_sequence
    #         i_5mer_chk_len=rcd[4]
    #         i_3mer_chk_len=rcd[5]
    #         s_5mer_clip=s_seq[:i_5mer_chk_len]
    #         s_3mer_clip=s_seq[-1*i_3mer_chk_len:]
    #     return l_slcted
####

    ####this is to check the sequence contain polyA signal (here we only select a small portion close to the tail)
    ###each record in format: t_segmt_pri, t_segmt_sup, algnmt, s_transduction_type, i_5mer_chk_len, i_3mer_chk_len
    ########t_segmt_pri in format (start, end, direction[+/-])
    # i_chk_win: is the windows size where polyA will be check within this window
    def call_transduction_from_polyA_signal(self, l_transduct, i_chk_win):
        l_contain_polyA = []
        polya = PolyA()
        for rcd in l_transduct:
            algnmt = rcd[2]
            seq = algnmt.query_sequence  ###this is already converted back to consensus

            s_transduction_type = rcd[3]
            i_5mer_chk_len = int(rcd[4])
            i_3mer_chk_len = int(rcd[5])
            s_5mer_polyA = ""
            s_3mer_polyA = ""
            if s_transduction_type == self._both_transdct:
                if i_5mer_chk_len > 0:
                    s_5mer_polyA = seq[:i_chk_win]
                if i_3mer_chk_len > 0:
                    s_3mer_polyA = seq[-1 * i_chk_win:]
            elif s_transduction_type == self._5mer:
                if i_5mer_chk_len > 0:
                    s_5mer_polyA = seq[:i_chk_win]
            elif s_transduction_type == self._3mer:
                if i_3mer_chk_len > 0:
                    s_3mer_polyA = seq[-1 * i_chk_win:]
            ####
            b_contain = False
            b_5mer=False
            if polya.contain_polyA_T(s_5mer_polyA, False):  # b_rc is alway False, as consensus is non-rc
                b_contain = True
                b_5mer=True
            if polya.contain_polyA_T(s_3mer_polyA, False):
                b_contain = True
            l_contain_polyA.append(b_contain)
        return l_contain_polyA
####