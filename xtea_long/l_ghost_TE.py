##03/21/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu


####this module will search for those ghost polymorphic full length copies, but doesn't show as an insertion
####Potential sources: 1. repeat copies come from the telomere region or gap region
###################### 2. population specific repeat copies (not in the reference genome - European biased)
##Major observation: 1) reads contain the copies will still map to some copy in the reference with two-side clip
##                   2) different copies can be distinguished using the flanking region

#####Criteria for two reads are originated from the same location:
####1. both flank sides are overlapped, and overlap length >1k, or
####2. One end has large overlap (>5kb), while the other side is short (not qualified, or unmapped)

####How about unmapped reads?
####How about duplications???!!!!

####
####
import os
import pysam
from x_reference import *
from rmsk_parser import *
from union_find_set import *
from l_asm import *
from x_contig import *
####
def unwrap_self_collect_polymorphic_copy(arg, **kwarg):
    return LNonRefPolymorphic.collect_polymorphic_by_chrm(*arg, **kwarg)

class LNonRefPolymorphic():
    def __init__(self, n_jobs, s_wfolder):
        self.n_jobs=n_jobs
        self.s_wfolder=s_wfolder

####
    def collect_ghost_polymorphic_rep_reads(self, sf_ref_copy_rmsk, i_min_rep_len, sf_bam_list, sf_ref, sf_out):
        '''
        :param sf_ref_copy_bed: This is the reference copy in bed format;
        :param sf_bam_list: The bam list;
        :param sf_out: output file
        :return: report the polymorphic copies
        '''
        sf_algnmt=self._get_algnmt_from_bam_list(sf_bam_list)
        if sf_algnmt is None or os.path.isfile(sf_algnmt)==False:
            print("Error: File {0} doesn't exist!!!\n".format(sf_algnmt))
            return

        samfile = pysam.AlignmentFile(sf_algnmt, "rb", reference_filename=sf_ref)
        references = samfile.references

        xchrom = XChromosome()
        l_chrm_records = []
        for chrm in references:
            if xchrom.is_decoy_contig_chrms(chrm) == True:  ###decoy sequnces and contigs are not considered
                continue
            sf_tmp_out = self.s_wfolder + chrm + "_tmp_candidate_polymorphic.fa"
            tmp_rcd=(chrm, sf_algnmt, sf_ref, sf_ref_copy_rmsk, i_min_rep_len, sf_tmp_out)
            l_chrm_records.append(tmp_rcd)
            #self.collect_polymorphic_by_chrm(tmp_rcd)
        samfile.close()
####
        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_polymorphic_copy, list(zip([self] * len(l_chrm_records), l_chrm_records)), 1)
        pool.close()
        pool.join()

        ####merge the candidate fasta files, and then classify based on the flanking regions
        sf_out_masked=sf_out+".masked.fa"
        sf_out_flank=sf_out+".flanking.fa"
        sf_out_lr_flanks = sf_out + ".separate_flanking.fa"##
        with open(sf_out, "w") as fout_merged, open(sf_out_masked, "w") as fout_mask, \
                open(sf_out_flank, "w") as fout_flank, open(sf_out_lr_flanks, "w") as fout_lr_flank:
            for chrm in references:
                sf_tmp_out = self.s_wfolder + chrm + "_tmp_candidate_polymorphic.fa"
                if os.path.isfile(sf_tmp_out):
                    with open(sf_tmp_out) as fin_tmp:
                        for line in fin_tmp:
                            fout_merged.write(line)
                sf_tmp_mask=sf_tmp_out + ".masked.fa"
                if os.path.isfile(sf_tmp_mask):
                    with open(sf_tmp_mask) as fin_mask:
                        for line in fin_mask:
                            fout_mask.write(line)
                sf_tmp_flank=sf_tmp_out + ".flanking.fa"
                if os.path.isfile(sf_tmp_flank):
                    with open(sf_tmp_flank) as fin_flank:
                        for line in fin_flank:
                            fout_flank.write(line)
                sf_tmp_flank2 = sf_tmp_out + ".separate_flanking.fa"
                if os.path.isfile(sf_tmp_flank2):
                    with open(sf_tmp_flank2) as fin_sprt_flank:
                        for line in fin_sprt_flank:
                            fout_lr_flank.write(line)
        ####
####
    ####Input: sf_flank: for left and right flank regions, they are saved as separate flank region
    ####If clip part is longer than "i_max_clip", then viewed as clip
    def cluster_reads_by_flank_region(self, sf_algnmt, sf_ori_fa, sf_flank, b_pacbio,
                                      i_max_clip, i_min_overlap,  iset_cutoff, sf_out_folder):
        if os.path.isdir(sf_out_folder)==False:
            os.mkdir(sf_out_folder)

        # 1. first align the flanking regions to them self
        if b_pacbio:
            self._algn_pacbio_reads_to_themselves(sf_flank, sf_flank, sf_algnmt)
        else:  # nanopore
            self._algn_nanopore_reads_to_themselves(sf_flank, sf_flank, sf_algnmt)
        # 2. cluster the reads based on the overalp between the flanking regions
        m_info, m_reads, l_reads = self._parse_self_aligned_reads(sf_algnmt, i_max_clip)
        m_cmpont = self._cluster_reads(m_info, m_reads, i_min_overlap)

        if sf_out_folder[-1] == "/":
            sf_out_folder += "/"
        sf_out_cluster = sf_out_folder + "tmp_component_ids.txt"
        sf_out_cluster2 = sf_out_folder + "tmp_component_ids2.txt"
        with open(sf_out_cluster, "w") as fout_compnt, open(sf_out_cluster2, "w") as fout_compnt2:
            for l_set in m_cmpont:  #
                # m_cluster[i_tmp]={}
                #l_set = m_cmpont[i_tmp]
                for i_id in l_set:
                    fout_compnt2.write(str(i_id) + "\t")
                fout_compnt2.write("\n")

                for i_id in l_set:
                    s_id = l_reads[i_id]
                    fout_compnt.write(s_id + "\t")
                fout_compnt.write("\n")
        # 3. get the sequence for each id and put each cluster in a separate file
        l_cluster_fa=self._export_cluster_seqs(sf_ori_fa, sf_out_cluster, iset_cutoff, self.s_wfolder)
        return l_cluster_fa
####

    ####for each cluster, assemble the collected reads
    def asm_collect_cluster_reads(self, l_cluster_list, sf_merged_copies):
        with open(sf_merged_copies, "w") as fout_merged:
            i_cnt=1
            for sf_fa in l_cluster_list:
                if os.path.isfile(sf_fa) == False:
                    continue
                lasm = L_Local_ASM()
                lasm.construct_short_cns_wtdbg2_given_file(sf_fa, self.n_jobs, self.s_wfolder)  # by default use 4 cores
                sf_lay_fa = "{0}_wtdbg2.ctg.lay.fa".format(sf_fa)
                if os.path.isfile(sf_lay_fa)==False:
                    continue
                s_slct_seq=self._load_in_fa_slct_longest(sf_lay_fa)
                if s_slct_seq!="" and len(s_slct_seq)>1:
                    s_id=">ghost_rep_copy{0}".format(i_cnt)
                    i_cnt+=1
                    fout_merged.write(s_id+"\n")
                    fout_merged.write(s_slct_seq+"\n")
####

    ####align the contigs to consensus to separate the copy and also the flanks
    def seprt_cns_flank_of_asm_contig(self, s_sample_id, sf_contig, sf_cns, sf_ref, sf_out_prefix):
        if len(sf_out_prefix)==0:
            print("Prefix is empty!!!")
            return
        if sf_out_prefix[-1]!="/":
            sf_out_prefix+="/"
        xtea_contig = XTEContig(self.s_wfolder, self.n_jobs)
        sf_algnmt=sf_out_prefix+"cns_2_contig.sorted.bam"
        xtea_contig.align_short_contigs_minimap2_v2(sf_cns, sf_contig, self.n_jobs, sf_algnmt)
        #parse the bam and call out the copies and flanking regions in separate files
        sf_rep=sf_out_prefix+"ghost_rep_copies.fa"#
        sf_flank=sf_out_prefix+"ghost_rep_copy_flanks.fa"
        self._parse_bam_to_separate_rep_flank(s_sample_id, sf_algnmt, sf_rep, sf_flank)
        #align the flank regions to reference genome
        sf_flank_algnmt=sf_flank+".algn_2_ref.sorted.bam"
        xtea_contig.align_contigs_2_reference_genome(sf_ref, sf_flank, self.n_jobs, sf_flank_algnmt)##

####
    def _parse_bam_to_separate_rep_flank(self, s_sample_id, sf_algnmt, sf_out_rep, sf_out_flanks):
        samfile = pysam.AlignmentFile(sf_algnmt, "rb")
        i_cnt=1
        with open(sf_out_rep, "w") as fout_rep, open(sf_out_flanks, "w") as fout_flanks:
            for algnmt in samfile.fetch():
                if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                    continue
                if algnmt.is_unmapped == True:  ##unmapped
                    continue
                if algnmt.is_duplicate == True:  ##duplciation
                    continue

                l_cigar = algnmt.cigar
                if len(l_cigar) < 2:  # wrong alignment
                    continue

                s_lflank=""
                i_seq_len = len(algnmt.query_sequence)
                i_map_start=0
                i_map_end=i_seq_len
                if l_cigar[0][0]==4:
                    s_lflank=algnmt.query_sequence[:l_cigar[0][1]]
                    i_map_start+=l_cigar[0][1]
                s_lflank_id = s_sample_id + "_" + str(i_cnt) + "_L"
                fout_flanks.write(">" + s_lflank_id + "\n")
                fout_flanks.write(s_lflank + "\n")

                s_rflank=""
                if l_cigar[-1][0]==4:
                    s_rflank=algnmt.query_sequence[:l_cigar[-1][1]]
                    i_map_end-=l_cigar[-1][1]
                s_rflank_id=s_sample_id+"_"+str(i_cnt)+"_R"
                fout_flanks.write(">"+s_rflank_id+"\n")
                fout_flanks.write(s_rflank + "\n")

                s_seq_copy=algnmt.query_sequence[i_map_start:i_map_end]
                s_copy_id=s_sample_id+"_"+str(i_cnt)+"_"+str(i_map_end-i_map_start)
                fout_rep.write(">"+s_copy_id+"\n")
                fout_rep.write(s_seq_copy+"\n")
                i_cnt+=1
        samfile.close()
####

####
    ####
    def _load_in_fa_slct_longest(self, sf_fa):
        n_max_len=0
        s_seq_slct=""
        with pysam.FastxFile(sf_fa) as freads:
            for entry in freads:
                s_tmp=entry.sequence
                i_tmp_len=len(s_tmp)
                if i_tmp_len > n_max_len:
                    n_max_len=i_tmp_len
                    s_seq_slct=s_tmp
        return s_seq_slct
    ####
    def _algn_pacbio_reads_to_themselves(self, sf_ref, sf_flank, sf_algnmt):
        cmd = "{0} -x ava-pb -c -a -t {1} {2} {3} | samtools view -hSb - | " \
              "samtools sort -o {4} -".format(global_values.MINIMAP2, self.n_jobs, sf_ref, sf_flank, sf_algnmt)
        cmd_runner = CMD_RUNNER()
        cmd_runner.run_cmd_small_output(cmd)

    def _algn_nanopore_reads_to_themselves(self, sf_ref, sf_flank, sf_algnmt):
        cmd = "{0} -x ava-ont -c -a -t {1} {2} {3} | samtools view -hSb - | " \
              "samtools sort -o {4} -".format(global_values.MINIMAP2, self.n_jobs, sf_ref, sf_flank, sf_algnmt)
        cmd_runner = CMD_RUNNER()
        cmd_runner.run_cmd_small_output(cmd)

    ####
    def collect_polymorphic_by_chrm(self, rcd):
        chrm=rcd[0]
        sf_bam=rcd[1]
        sf_ref=rcd[2]
        sf_rmsk=rcd[3]
        min_copy_len=rcd[4]
        sf_out=rcd[5]

        sf_out_masked=sf_out+".masked.fa"
        sf_out_flanking = sf_out + ".flanking.fa"
        sf_out_flanks=sf_out + ".separate_flanking.fa"

        b_tmplt_with_chr=False
        if len(chrm)>3 and chrm[:3]=="chr":
            b_tmplt_with_chr=True

        #m_rep_ref_copies=self._load_in_bed(sf_bed, b_tmplt_with_chr)
        rmsk_parser=RMSK_Parser(sf_rmsk)
        m_rep_ref_copies = rmsk_parser.slct_copies_with_min_len(min_copy_len, b_tmplt_with_chr)

        if chrm not in m_rep_ref_copies:
            return
        with open(sf_out, "w") as fout_slct_reads, open(sf_out_masked, "w") as fout_masked, \
                open(sf_out_flanking, "w") as fout_flank, open(sf_out_flanks,"w") as fout_lr_flanks:

            samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=sf_ref)
            l_regions=[]
            for istart in m_rep_ref_copies[chrm]:
                iend=m_rep_ref_copies[chrm][istart][0]
                b_rc=m_rep_ref_copies[chrm][istart][1]
                l_regions.append((istart, iend, b_rc))
            ####
            for (istart, iend, b_rc) in l_regions:
                rgn_start = istart - global_values.LRD_CLIP_SEARCH_WIN
                if rgn_start < 0:
                    rgn_start = 0
                rgn_end = iend + global_values.LRD_CLIP_SEARCH_WIN
                for algnmt in samfile.fetch(chrm, rgn_start, rgn_end):##fetch reads mapped to "chrm"
                    ##here need to skip the secondary and supplementary alignments?
                    # if algnmt.is_secondary or algnmt.is_supplementary:
                    if algnmt.is_supplementary or algnmt.is_secondary:  # skip supplementary, but keep the secondary alignment
                        continue
                    if algnmt.is_duplicate == True:  ##duplciate
                        continue
                    if algnmt.is_unmapped == True:  ##unmapped
                        continue
                    l_cigar = algnmt.cigar
                    if len(l_cigar) < 2:  # wrong alignment
                        continue
                    # if algnmt.mapping_quality < global_values.LRD_MIN_MAPQ:  #by default this is set to 20
                    #     continue

                    map_pos = algnmt.reference_start
                    query_name = algnmt.query_name ####note, this one may have space in the middle
                    query_seq = algnmt.query_sequence
                    if query_seq is None:
                        continue
                    # query_quality = algnmt.query_qualities  ##this is different from the one saved in the fastq/sam, no offset 33 to subtract

                    #find the left and right clip positions
                    i_lside=map_pos
                    i_rside, n_mapped=self._get_right_side_pos_from_cigar(map_pos, l_cigar)

                    #if both side clipped and both clip position are close to the repeat copy breakpoints
                    if l_cigar[0][0] == 4 and l_cigar[-1][0] == 4:  # both side soft-clipped
                        i_lclip_len = l_cigar[0][1]
                        i_rclip_len = l_cigar[-1][1]
                        if i_lclip_len < global_values.LRD_POLYMORPHIC_MIN_CLIP_LEN and \
                                        i_rclip_len < global_values.LRD_POLYMORPHIC_MIN_CLIP_LEN:
                            continue
                        #check the clip positions against the copy breakpoints
                        b_lclose=self._is_clip_pos_close_copy_end(i_lside, istart, global_values.LRD_POLYM_BRK_CHK_WIN)
                        b_rclose=self._is_clip_pos_close_copy_end(i_rside, iend, global_values.LRD_POLYM_BRK_CHK_WIN)
                        if b_lclose and b_rclose:
                            #save the read
                            l_tmp_fields=query_name.split()
                            s_new_name="~".join(l_tmp_fields)#read name

                            s_rc="+"
                            if b_rc==True:
                                s_rc="C"
                            sinfo=">"+s_new_name+"_{0}_{1}_{2}_{3}_{4}_{5}_{6}".format(chrm, istart, iend, s_rc,
                                                                                       i_lclip_len, i_rclip_len, n_mapped)
                            #sinfo2=sinfo+" "+chrm+":"+str(map_pos)+" "+algnmt.cigarstring+"\n"
                            sinfo2 = sinfo + " " + chrm + ":" + str(map_pos) + "\n"
                            fout_slct_reads.write(sinfo2)
                            fout_slct_reads.write(str(query_seq)+"\n")

                            i_mask_len=len(query_seq)-i_lclip_len-i_rclip_len
                            s_masked_read=query_seq[:i_lclip_len]+ ("N"*i_mask_len) +query_seq[-1*i_rclip_len:]
                            fout_masked.write(sinfo2)
                            fout_masked.write(s_masked_read+"\n")
                            fout_flank.write(sinfo2)
                            fout_flank.write(query_seq[:i_lclip_len]+query_seq[-1*i_rclip_len:]+"\n")

                            if i_lclip_len >= global_values.LRD_POLYMORPHIC_MIN_CLIP_LEN:
                                fout_lr_flanks.write(sinfo+"_L"+ " " + chrm + ":" + str(map_pos) + "\n")
                                fout_lr_flanks.write(query_seq[:i_lclip_len]+"\n")

                            if i_rclip_len >= global_values.LRD_POLYMORPHIC_MIN_CLIP_LEN:
                                fout_lr_flanks.write(sinfo + "_R" + " " + chrm + ":" + str(map_pos) + "\n")
                                fout_lr_flanks.write(query_seq[-1*i_rclip_len:] + "\n")
####
####
                    # # Or one side clipped, the the right side of the read also close to the repeat copy breakpoints
                    # if l_cigar[-1][0] == 4:  # right clipped
                    #     rclip_len = l_cigar[-1][1]
                    #     if rclip_len > global_values.LRD_MIN_CLIP_LTH:
                    #         clip_pos = rclip_pos  # this is calculated from cigar
                    #         if rclip_pos != -1:
                    #             self._save_pos_to_dict(chrm, clip_pos, m_breakpoint_pos, 1)
            samfile.close()
####
####
    #Given reads id by cluster, output the sequences of each cluster
    def _export_cluster_seqs(self, sf_ori_reads, sf_comp_id, iset_cutoff, sf_folder):
        l_cluster_fa=[]
        m_reads={}
        with pysam.FastxFile(sf_ori_reads) as freads:
            for entry in freads:
                tmp_fields=entry.name.split("_")
                read_ori_name="_".join(tmp_fields[:-7]) #note there is no direction, so -7 not -8
                m_reads[read_ori_name]=entry.sequence
        icnt=0
        with open(sf_comp_id) as fin_comp:
            for line in fin_comp:
                fields = line.split()
                if len(fields) < iset_cutoff:
                    continue
                sf_out_tmp=sf_folder+"ghost_cluster{0}.fa".format(icnt)
                icnt+=1
                l_cluster_fa.append(sf_out_tmp)
                with open(sf_out_tmp, "w") as fout_tmp:
                    for s_tmp in fields:
                        if s_tmp not in m_reads:
                            print("Error ", s_tmp, "not found in reads!")
                            continue
                        fout_tmp.write(">"+s_tmp+"\n")
                        fout_tmp.write(m_reads[s_tmp]+"\n")
        return l_cluster_fa
####
####
    ####parse the alignment to cluster reads
    def _parse_self_aligned_reads(self, sf_bam, i_max_clip):
        m_cluster={} #this is to save those have two sides flanking overlap
        m_reads={}
        l_reads=[]
        i_cnt=0
        samfile = pysam.AlignmentFile(sf_bam, "rb")
        for algnmt in samfile.fetch(until_eof=True):  ##fetch reads
            if algnmt.is_supplementary:  # skip supplementary, but keep the secondary alignment
                continue
            if algnmt.is_unmapped == True:  ##unmapped
                continue

            map_pos = algnmt.reference_start #
            query_name = algnmt.query_name
            l_cigar = algnmt.cigar
            ref_name=algnmt.reference_name
            if query_name == ref_name: #skip those aligned to itself
                continue

            b_rc = algnmt.is_reverse
            #if one contig is contained within another, then they are grouped
            if self._is_contained(l_cigar, i_max_clip)==True:
                qname_ori, qname_flank_dir, s_rc_qrep = self._parse_info_from_name(query_name)
                ref_ori, ref_flank_dir, s_rc_rrep = self._parse_info_from_name(ref_name)
                if b_rc == True and (s_rc_qrep == s_rc_rrep):  # conflict: rc, but the repeats are not rc
                    continue
                if b_rc == False and (s_rc_qrep != s_rc_rrep):  # conflict: not-rc, but the repeats are rc
                    continue
                if b_rc == True and (qname_flank_dir == ref_flank_dir):  # conflict: rc, but RR/LL
                    continue
                if b_rc == False and (qname_flank_dir != ref_flank_dir):  # conflict: not-rc, but RL/LR
                    continue

                s_id_1 = "{0}~{1}".format(ref_ori, qname_ori)
                s_id_2 = "{0}{1}".format(ref_flank_dir, qname_flank_dir)
                if qname_ori>ref_ori:
                    s_id_1="{0}~{1}".format(qname_ori, ref_ori)
                    s_id_2 = "{0}{1}".format(qname_flank_dir, ref_flank_dir)

                ##if the mapped region is long, then save to m_cluster_long_overlap
                i_rside, i_overlap = self._get_right_side_pos_from_cigar(map_pos, l_cigar)
                if s_id_1 in m_cluster:
                    s_exist_id_2=m_cluster[s_id_1][0]
                    s_chk=self._get_rc_direction(s_id_2)
                    if s_chk != s_exist_id_2:
                        continue
                else:
                    m_cluster[s_id_1]=[]
                m_cluster[s_id_1].append((s_id_2, qname_ori, ref_ori, i_overlap))

                if qname_ori not in m_reads:
                    m_reads[qname_ori]=i_cnt
                    l_reads.append(qname_ori)
                    i_cnt+=1
                if ref_ori not in m_reads:
                    m_reads[ref_ori]=i_cnt
                    l_reads.append(ref_ori)
                    i_cnt+=1

        samfile.close()
        return m_cluster, m_reads, l_reads

####
    def _cluster_reads(self, m_info, m_reads, i_min_overlap):
        n_reads=len(m_reads)
        ufs=UnionFindSet(n_reads)
        ufs.setIdSz()
        for s_id in m_info:
            #here require both flank are overlapped
            if len(m_info[s_id])==1:#if only have one side overlap, then require the overlpa length is long
                i_overlap_len=m_info[s_id][0][3]
                if i_overlap_len<i_min_overlap:
                    continue
            elif len(m_info[s_id])>2:
                print("Wrong alignment: has more than 2 overlaps between two reads!!!\n")
                continue
            qname=m_info[s_id][0][1]
            rname=m_info[s_id][0][2]
            i_qid=m_reads[qname]
            i_rid=m_reads[rname]
            ufs.union(i_qid, i_rid)
        m_cmpont=ufs.outputComponents()
        return m_cmpont

####
    def _get_rc_direction(self, s_dir):
        if s_dir=="RR":
            return "LL"
        elif s_dir=="LL":
            return "RR"
        elif s_dir=="LR":
            return "RL"
        else:
            return "LR"
    #
    def _parse_info_from_name(self, s_name):
        qname_fields = s_name.split("_")
        # chrm, istart, iend, s_rc, i_lclip_len, i_rclip_len, n_mapped, L/R
        s_rc_rep = qname_fields[-5]
        qname_ori = "_".join(qname_fields[:-8])
        qname_flank_ori = qname_fields[-1]
        return qname_ori, qname_flank_ori, s_rc_rep

    ####one read is contained in another
    def _is_contained(self, l_cigar, i_max_clip):
        b_left_clip=False
        b_right_clip=False
        if len(l_cigar)==1:
            return True

        if l_cigar[0][0] == 4:  # left-clip
            lclip_len = l_cigar[0][1]
            if lclip_len>i_max_clip:
                b_left_clip=True
        if l_cigar[-1][0] == 4:  # right-clip
            lclip_len = l_cigar[-1][1]
            if lclip_len > i_max_clip:
                b_right_clip = True
        if b_left_clip or b_right_clip:
            return False
        return True

####
    #clip positon is close to the copy end
    def _is_clip_pos_close_copy_end(self, i_clip_pos, i_copy_end, i_slack):
        if abs(i_clip_pos-i_copy_end) <= i_slack:
            return True
        else:
            return False

    #load in the repeat copy information saved in bed format
    def _load_in_bed(self, sf_bed, b_tmplt_with_chr):
        m_copies={}
        with open(sf_bed) as fin_bed:
            for line in fin_bed:
                fields=line.split("\t")
                if len(fields)<3:
                    print("wrong bed format")
                    continue
                if line[0]=="#":
                    continue
                ori_chrm=fields[0]
                chrm=self._process_chrm_name(b_tmplt_with_chr, ori_chrm)
                istart=int(fields[1])
                iend=int(fields[2])

                if chrm not in m_copies:
                    m_copies[chrm]={}
                m_copies[chrm][istart]=iend
        return m_copies

####
    ####
    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, b_tmplt_with_chr, chrm):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        if b_tmplt_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_tmplt_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_tmplt_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    ####get the right side mapping position
    def _get_right_side_pos_from_cigar(self, map_pos, l_cigar):
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
            elif opn == 4 or opn == 5 or opn == 6 or opn == 1:  # hard-clip (H) or padding (P) or soft-clip (S)
                continue
        return ref_pos, n_mapped
####
####
    def _get_algnmt_from_bam_list(self, sf_bam_list):
        if os.path.isfile(sf_bam_list)==False:
            print("Error: File {0} doesn't exist!!!\n".format(sf_bam_list))
            return None
        sf_algnmt=None
        with open(sf_bam_list) as fin_bam_list:
            for line in fin_bam_list:
                fields=line.split()
                sf_algnmt=fields[-1]
        return sf_algnmt

####
####read in the rmsk output file, and select the potential full length L1s that within the centromere regions
    #Here only check "Alpha", "Beta" and "HSATII" repeats in the flanking regions
    def slct_centromere_with_flank_rmsk(self, sf_rmsk, sf_out):
        rmsk_parser=RMSK_Parser(sf_rmsk)
        #m_rcds_sites in format: [ins_chrm][ins_pos].append(hit_rcd)
        #m_rcds_ctg in format: [ins_chrm][ins_pos][s_contig].append(hit_rcd)
        #hit_rcd in format: (i_contig_start, i_contig_end, b_rc, s_sub_family, s_rep_family, i_cns_start, i_cns_end)
        m_rcds_sites, m_rcds_ctg=rmsk_parser.parse_rmsk()
        m_slcted={}
        #for chrm 
        pass