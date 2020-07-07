##11/17/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu, chong.simon.chu@gmail.com

#To-do-list:
# 1. separate HERV insertions from others (only HERV with duplication cases)
# 2. save to a separate fasta file (candidate dimorphic HERV), need to run l_dimorphic_HERV step after this.

import os
import pysam
from l_contig import *
from l_rep_classification import *
from l_vcf import *

class Complex_TE_SV(LRD_Contig):
    def __init__(self, s_working_folder, n_jobs):
        LRD_Contig.__init__(self, s_working_folder, n_jobs)
        self.TBD_TE_SEQ = "to_be_decide_TE_seq.fa"

    ####
    def parse_trim_seq_from_flank_contig_algnmt(self, l_rcds, f_map_ratio, sf_ltrimmed, sf_rtrimmed):
        m_lcontig_len = {}
        m_rcontig_len = {}

        with open(sf_ltrimmed, "w") as fout_lseq, open(sf_rtrimmed, "w") as fout_rseq:
            for (sf_asm, sf_flank_fa, sf_algnmt, chrm, pos) in l_rcds:
                if os.path.isfile(sf_algnmt) == False:
                    continue
                b_left, s_contig, i_seq_start, i_seq_end = self.parse_trim_one_algnmt(sf_algnmt, f_map_ratio)
                if i_seq_start == 0 and i_seq_end == 0:
                    continue

                if os.path.isfile(sf_asm)==False:
                    print("File %s doesn't exist!" % sf_asm)
                    continue
                m_contig=self._load_in_fasta(sf_asm)
                if s_contig not in m_contig:
                    print("Couldn't find contig %s!" % s_contig)
                    continue
                s_seq1 = m_contig[s_contig]
                s_segmt = s_seq1[i_seq_start:i_seq_end]
                if len(s_segmt)<50:
                    print "Debug", len(s_seq1), i_seq_start, i_seq_end
                if len(s_segmt)==0:
                    continue

                s_id = "{0}{1}{2}".format(chrm, global_values.SEPERATOR, pos)
                if b_left == True:
                    fout_lseq.write(">" + s_id + "\n" + s_segmt + "\n")
                    m_lcontig_len[s_id] = len(s_segmt)
                else:
                    fout_rseq.write(">" + s_id + "\n" + s_segmt + "\n")
                    m_rcontig_len[s_id] = len(s_segmt)
        return m_lcontig_len, m_rcontig_len
####

    ####
    # if b_change_ori is set, then it means we are aligning left (right) flank to right(left) trimmed contig
    def parse_trim_one_algnmt(self, sf_algnmt, f_map_ratio, b_change_ori=False):
        i_seq_start = 0
        i_seq_end = 0
        s_contig = ""
        b_left = True
        b_rc = False
        if os.path.isfile(sf_algnmt)==False:
            return b_left, s_contig, i_seq_start, i_seq_end
        try:
            samfile = pysam.AlignmentFile(sf_algnmt)  # read in the sam file
            for algnmt in samfile.fetch():  # check each alignment, and find "left" and "right" flank
                if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary
                    continue
                if algnmt.is_unmapped == True:  ##unmapped
                    continue
                if algnmt.is_duplicate == True:  ##duplciate
                    continue
                # as  the clipped long reads are from the alignment, so already corrected from rc.
                # Note, we are not sure about the assembler whether it generate a reverse ones
                b_rc = algnmt.is_reverse
                if algnmt.is_reverse == True and b_change_ori == False:
                    # if this is align the flank to get the trimmed seq, then require it non reverse-complementary
                    continue

                if algnmt.mapping_quality < global_values.MINIMUM_DISC_MAPQ:  ###############anchor mapping quality
                    continue
                # make sure most region is aligned
                l_cigar = algnmt.cigar

                # flank length
                i_flk_len = len(algnmt.query_sequence)
                i_max_clip = int(i_flk_len * (1 - f_map_ratio))  # at most half is clipped out

                i_lclip_len = 0  # left-clip length
                i_rclip_len = 0
                if l_cigar[0][0] == 5 or l_cigar[0][0] == 4:
                    i_lclip_len = l_cigar[0][1]
                if l_cigar[-1][0] == 5 or l_cigar[-1][0] == 4:
                    i_rclip_len = l_cigar[-1][1]

                i_clip_total_len = i_lclip_len + i_rclip_len  # left flank total clip length
                if i_clip_total_len > i_max_clip:
                    continue

                map_pos = algnmt.reference_start
                s_qname = algnmt.query_name
                s_contig = algnmt.reference_name
                s_tmp_ref = global_values.LEFT_FLANK
                if b_rc == True:
                    s_tmp_ref = global_values.RIGHT_FLANK
                if s_qname == s_tmp_ref:#
                    map_pos += (i_flk_len - i_lclip_len - i_rclip_len)
                    i_seq_start = map_pos
                    i_seq_end = -1
                else:
                    i_seq_start = 0
                    i_seq_end = map_pos
                    b_left = False
            samfile.close()
        except ValueError:
            print sf_algnmt, "is empty"

        if b_change_ori == True:
            b_left = b_rc
        return b_left, s_contig, i_seq_start, i_seq_end
        ####

    def parse_flk_2_trim_algnmt_for_internal_seq(self, l_rcd, f_map_ratio):
        m_seqs = {}
        b_change_ori = True
        for (sf_r_trimmed, sf_flank_fa, sf_algnmt, chrm, pos) in l_rcd:
            b_rc, s_contig, i_seq_start, i_seq_end = self.parse_trim_one_algnmt(sf_algnmt, f_map_ratio,
                                                                                b_change_ori)
            if i_seq_start == 0 and i_seq_end == 0:
                continue
            s_id = "%s%s%d" % (chrm, global_values.SEPERATOR, pos)
            m_seqs[s_id] = (s_contig, i_seq_start, i_seq_end)
        return m_seqs
        ####

    #classify the complex SV to detailed events
    #each record in format: "s_id:(s_contig_id, i_seq_start, i_seq_end)"
    #m_l_seq_rcds: record of left-flank to right-trimmed contigs
    def classify_complex_SV(self, m_l_seq_rcds, m_r_seq_rcds, sf_l_trimmed, sf_r_trimmed, ratio_cutoff):
        m_lcontig = self._load_in_fasta(sf_l_trimmed)#left-flank already trimmed contigs
        m_rcontig=self._load_in_fasta(sf_r_trimmed)#right flank already trimmed contigs
        m_simple_sv={}#save the simple deletion/translocation
        m_candidates={}#candidate sites
        for s_l_id in m_l_seq_rcds:
            # 1. the two breakpoints should match each other, and the parsed seq length should be similar
            (s_r_id, i_seq_start_r, i_seq_end_r)=m_l_seq_rcds[s_l_id]
            if s_r_id not in m_r_seq_rcds:
                print s_r_id, "not found in right rcds!"
                continue
            (s_l_id_tmp, i_seq_start_l, i_seq_end_l) = m_r_seq_rcds[s_r_id]
            if s_l_id_tmp != s_l_id:
                print s_l_id_tmp, s_l_id, "are not equal!"
                continue

            #check whether exist
            if s_l_id not in m_lcontig:
                print "%s is not found in dict!" % (s_l_id)
                continue
            if s_r_id not in m_rcontig:
                print "%s is not found in dict!" % (s_r_id)
                continue
            s_l_seq=m_lcontig[s_l_id][i_seq_start_l:i_seq_end_l]
            s_r_seq=m_rcontig[s_r_id][i_seq_start_r:i_seq_end_r]#
            i_lseq_len=len(s_l_seq)
            i_rseq_len=len(s_r_seq)

            s_new_id = s_l_id + global_values.SEPERATOR + s_r_id
            #checkc length consistency
            #2. the length should be larger then minimal length
            if (i_lseq_len<global_values.LRD_MIN_INS_LTH) or (i_rseq_len<global_values.LRD_MIN_INS_LTH):
                m_simple_sv[s_new_id]=(s_l_id, s_r_id)#either simple deletion or translocation
                continue

            #separate chrm and pos fields
            l_tmp1=s_l_id.split(global_values.SEPERATOR)
            s_lpos=l_tmp1[-1]
            s_lchrm=global_values.SEPERATOR.join(l_tmp1[:-1])

            l_tmp2=s_r_id.split(global_values.SEPERATOR)
            s_rpos=l_tmp2[-1]
            s_rchrm=global_values.SEPERATOR.join(l_tmp2[:-1])
            m_candidates[s_new_id]=(s_lchrm, s_lpos, s_rchrm, s_rpos, s_l_seq, s_r_seq)
            #print i_seq_start_l, i_seq_end_l, i_seq_start_r, i_seq_end_r, s_new_id, "seq_pos"
            #print s_new_id, s_l_seq, s_r_seq
        return m_candidates, m_simple_sv
####

####
    def classify_events_by_seq_algnmt_mask(self, m_candidates, sf_ref, f_mratio, f_cover_ratio, sf_l1_rslt, sf_sva_rslt,
                                           flk_lenth, sf_wfolder, sf_rep_folder, b_hg19, m_simple_sv, sf_out_prefix):
        l_inter, l_intra=self.define_inter_intra_chrm_events(m_candidates)
        #for intra events
        m_for_te_mask_intra, m_te_del, m_te_inv, m_te_dup, m_non_te_sv, m_algnmts=self.classify_intra_events(
            sf_ref, l_intra, f_mratio, f_cover_ratio, sf_wfolder)
        l_rslts=self.mask_seq_for_TE(m_for_te_mask_intra, sf_l1_rslt, sf_sva_rslt, sf_ref, flk_lenth, sf_wfolder,
                        sf_rep_folder, b_hg19, sf_out_prefix)
        #select the qualified intra candidates
        sf_te_promoted_rslt = sf_out_prefix + "_complex_sv.txt"
        m_te_sv=self.slct_te_complex_sv_from_mask(l_rslts, m_te_del, m_te_inv, m_te_dup, m_algnmts, sf_te_promoted_rslt)
        #for inter events
        #no candidates???

        #save simple sv events and other complex events
        sf_non_te_related = sf_out_prefix + "_no_TE_related_SV.txt"
        self.save_non_te_related_SV(m_for_te_mask_intra, m_non_te_sv, m_algnmts, m_te_sv, m_te_dup, m_simple_sv, sf_non_te_related)

####
    def save_non_te_related_SV(self, m_all_complex, m_non_te_sv, m_algnmts, m_te_sv, m_te_dup, m_simple_sv, sf_out):
        with open(sf_out,"w") as fout_rslt, open(sf_out+".fa","w") as fout_fa:
            #first output the other complex events
            for s_id in m_all_complex:
                if s_id in m_te_sv:
                    continue
                t_rcd = m_algnmts[s_id]  # in format: (sf_algnmt, sf_contig, s_contig_seq)
                s_contig_seq = t_rcd[2]
                #s_type=m_te_sv[s_id]
                fout_rslt.write(s_id + "\tnon_TE\t"+ "Complex_with_insertion" + "\t" + s_contig_seq + "\n")
                if s_id in m_te_dup:#save for further HERV checking
                    fout_fa.write(">"+s_id+"\n"+s_contig_seq+"\n")
            #second, output the simple non-TE-SV events
            for s_id in m_non_te_sv:
                t_rcd = m_algnmts[s_id]  # in format: (sf_algnmt, sf_contig, s_contig_seq)
                s_contig_seq = t_rcd[2]
                s_type = m_non_te_sv[s_id]
                fout_rslt.write(s_id + "\tnon_TE\t" + s_type + "\t" + s_contig_seq + "\n")

            #third, output the simple SV (deletion or translocation)
            for s_id in m_simple_sv:
                fout_rslt.write(s_id + "\tsimple\t" + "Deletion_or_Translocation" + "\t" + "." + "\n")
####
####
    #for intra events, separet by whether it's duplication or (deletion, inversion)
    def classify_intra_events(self, sf_ref, l_intra, f_map_ratio, f_full_cover_ratio, sf_wfolder):
        if sf_wfolder[-1]!="/":
            sf_wfolder+="/"

        m_te_del={}
        m_te_dup={}
        m_te_inv={}
        m_simple_sv={}#

        xref = XReference()
        #first get the reference sequence between the two breakpoints
        l_ref_seqs = xref.get_ref_seqs_of_sites(sf_ref, l_intra)
        #save to file for each seq, and return the path
        m_ref_seqs, m_contigs, m_contig_seq = self.dump_ref_seq_to_file(l_intra, l_ref_seqs, sf_wfolder)
        #prepare the alignment record, and align the sequences
        l_algn_rcd=[]
        m_algnmt={}
        for s_id in m_ref_seqs:
            sf_ref_sgmt=m_ref_seqs[s_id]
            sf_contig=m_contigs[s_id]
            sf_algnmt=sf_wfolder+s_id+".ref_2_contig.bam"
            l_tmp1 = s_id.split(global_values.SEPERATOR)
            s_lpos = l_tmp1[-1]
            s_chrm = global_values.SEPERATOR.join(l_tmp1[:-1])
            l_algn_rcd.append((sf_contig, sf_ref_sgmt, sf_algnmt, s_chrm, int(s_lpos)))
            s_contig_seq=m_contig_seq[s_id]
            m_algnmt[s_id]=(sf_algnmt, sf_contig, s_contig_seq)
        self.align_flanks_to_contig_2(l_algn_rcd)

####first, collect all the sequences need to be masked
####second, run the mask module (same module as the classification module)
####third, classify based on the mask results

        m_for_te_mask={}
        m_sv_type={}
        i_idx=0
        for (s_chrm, s_lpos, s_rpos, s_l_seq, s_r_seq) in l_intra:
            # s_l_id = "%s%s%s" % (s_chrm, global_values.SEPERATOR, s_lpos)
            # s_r_id = "%s%s%s" % (s_chrm, global_values.SEPERATOR, s_rpos)
            s_ref_sgmt=l_ref_seqs[i_idx]
            i_idx+=1
            s_id2 = s_chrm + global_values.SEPERATOR + s_lpos + global_values.SEPERATOR + s_rpos
####        #Here duplicate breakpoints means: right (left) clip formed breakpoint at right (left)

            if int(s_lpos)>int(s_rpos):  # this is for duplication (including tandem duplication)
                #2. find a hit
                sf_algnmt_tmp = m_algnmt[s_id2][0]  # alignment is generated by algning ref-segmt to contig
                s_contig = m_algnmt[s_id2][2]  #
                b_dup=True#if this is duplication, then no mapq cutoff
                t_rcd = self.parse_segmt_contig_algnmt(sf_algnmt_tmp, s_contig, f_map_ratio, f_full_cover_ratio, b_dup)
                #3. mask the rest of the sequence
                if t_rcd is None:
                    if len(s_contig)<global_values.LRD_MIN_INS_LTH:
                        m_simple_sv[s_id2]="tandem_duplication"
                        continue
                    #this means a duplication
                    m_te_dup[s_id2]=t_rcd
                    m_for_te_mask[s_id2] = s_contig
                    m_sv_type[s_id2]="Insertion_with_duplication"

                #3.1 if no TE seq, then this is just a tandem duplication
                #3.2, also pay attention to full-length HERV, if the reference genome is LTR
            else:#this is for deletion or inversion
                #align the ref-segmt to contig, and if aligned with reverse-complementary, then inversion candidates
                sf_algnmt_tmp=m_algnmt[s_id2][0]#alignment is generated by algning ref-segmt to contig
                s_contig=m_algnmt[s_id2][2]#
                t_rcd=self.parse_segmt_contig_algnmt(sf_algnmt_tmp, s_contig, f_map_ratio, f_full_cover_ratio)
                #if not aligned, then further check whether the segmnt is a TE, if yes, then that's a promoted del
                if t_rcd is None:#
                    if len(s_contig)<global_values.LRD_MIN_INS_LTH:
                        m_simple_sv[s_id2] = "deletion"
                        continue
                    m_for_te_mask[s_id2]=s_contig#save the contig for further TE masking
                    m_sv_type[s_id2] = "Insertion_with_deletion"
                    m_te_del[s_id2]=t_rcd#save the id
                else:
                    #rcd in format:(b_rc, i_map_start, i_map_end, len(s_contig))
                    (b_rc, i_map_start, i_map_end, i_ctg_len)=t_rcd
                    #get the unaligned sequence for masking
                    s_seq=s_contig[i_map_end:]
                    if i_map_start>(i_ctg_len-i_map_end):
                        s_seq=s_contig[:i_map_start]

                    if b_rc==True:#check for potential TE promoted inversion
                        if len(s_seq) < global_values.LRD_MIN_INS_LTH:
                            m_simple_sv[s_id2] = "inversion"
                            continue
                        m_te_inv[s_id2]=t_rcd
                        m_for_te_mask[s_id2]=s_seq
                        m_sv_type[s_id2] = "Insertion_with_inversion"
####
        #save the tbd sequence to a file
        sf_tbd_te_seq=sf_wfolder+self.TBD_TE_SEQ
        with open(sf_tbd_te_seq,"w") as fout_tbd_seq:
            for s_tmp_id in m_for_te_mask:
                s_tmp_seq=m_for_te_mask[s_tmp_id]
                s_type=m_sv_type[s_tmp_id]
                fout_tbd_seq.write(">"+s_tmp_id+"\t"+s_type+"\n"+s_tmp_seq+"\n")

        print m_te_del, "TE_deletion"
        print m_te_inv, "TE_inversion"
        print m_te_dup, "TE_duplication"
        return m_for_te_mask, m_te_del, m_te_inv, m_te_dup, m_simple_sv, m_algnmt#sequence for masking
####
        #if aligned with reverse complementary:
        #1. the whole contig is covered, then it's purely an inversion
        #2. the remaining part is TE seq, then it's a promoted inversion
        #print t_rcd, s_chrm, s_lpos, s_rpos
####

####
####mask the sequence to specific TE type
####Here to also mask the transduction regions, we integrate all the already masked L1 and SVA
    def mask_seq_for_TE(self, m_for_te_mask, sf_l1_rslt, sf_sva_rslt, sf_ref, flk_lenth, swfolder,
                        sf_rep_folder, b_hg19, sf_out_prefix):
        # merge the fasta files
        lrr = L_Raw_Rslt()
        m_l1_seqs = lrr.load_in_results(sf_l1_rslt)
        m_sva_seqs=lrr.load_in_results(sf_sva_rslt)
        sf_rep_ins=swfolder+"TE_promoted_complex_SV_to_be_masked.fa"
        with open(sf_rep_ins, "w") as fout_tbd_seq:
            for s_tmp_id in m_for_te_mask:
                s_tmp_seq=m_for_te_mask[s_tmp_id]
                fout_tbd_seq.write(">"+s_tmp_id+"\n"+s_tmp_seq+"\n")
            for s_id_tmp in m_l1_seqs:#
                s_seq_tmp=m_l1_seqs[s_id_tmp]
                fout_tbd_seq.write(">"+s_id_tmp+"\n"+s_seq_tmp+"\n")
            for s_id_tmp in m_sva_seqs:#
                s_seq_tmp=m_sva_seqs[s_id_tmp]
                fout_tbd_seq.write(">"+s_id_tmp+"\n"+s_seq_tmp+"\n")

        #mask the sequences
        lrc = LRepClassification(swfolder, self.n_jobs)
        i_type = lrc.get_mask_complex_sv_rep_type()
        lrc.set_rep_configuration(i_type, sf_rep_folder, b_hg19)
        lrc.classify_ins_seqs(sf_rep_ins, sf_ref, flk_lenth, sf_out_prefix)
        l_rslt_list=lrc.get_merged_results(sf_out_prefix, i_type)
        return l_rslt_list

### #First, select those sites with seq that masked as TE
### #Second, separate by type. Pay attention to HERV for "duplication" events
    def slct_te_complex_sv_from_mask(self, l_rslts, m_te_del, m_te_inv, m_te_dup, m_algnmts, sf_rslt):
        m_te_del2=self._process_id_in_dict(m_te_del)#only take the first two fields in the id
        m_te_inv2 = self._process_id_in_dict(m_te_inv)
        m_te_dup2 = self._process_id_in_dict(m_te_dup)
        m_slcted={}
        lrr = L_Raw_Rslt()
        with open(sf_rslt, "w") as fout_rslt:
            for (s_type, sf_one_rep) in l_rslts:
                if os.path.isfile(sf_one_rep)==False:
                    continue
                m_rslt_seqs = lrr.load_in_results(sf_one_rep)
                for s_tmp_id in m_rslt_seqs:
                    if s_tmp_id in m_te_del2:
                        s_id=m_te_del2[s_tmp_id][0]
                        m_slcted[s_id]=1
                        t_rcd=m_algnmts[s_id]#in format: (sf_algnmt, sf_contig, s_contig_seq)
                        s_contig_seq=t_rcd[2]
                        fout_rslt.write(s_id+"\t"+s_type+"\tTE_promoted_deletion\t"+s_contig_seq+"\n")
                    elif s_tmp_id in m_te_inv2:
                        s_id = m_te_inv2[s_tmp_id][0]
                        m_slcted[s_id] = 1
                        t_rcd = m_algnmts[s_id]  # in format: (sf_algnmt, sf_contig, s_contig_seq)
                        s_contig_seq = t_rcd[2]
                        fout_rslt.write(s_id + "\t" + s_type + "\tTE_promoted_inversion\t"+s_contig_seq+"\n")
                    elif s_tmp_id in m_te_dup2:
                        s_id = m_te_dup2[s_tmp_id][0]
                        m_slcted[s_id] = 1
                        t_rcd = m_algnmts[s_id]  # in format: (sf_algnmt, sf_contig, s_contig_seq)
                        s_contig_seq = t_rcd[2]
                        fout_rslt.write(s_id + "\t" + s_type + "\tTE_promoted_duplication\t"+s_contig_seq+"\n")
        return m_slcted
####

    def _process_id_in_dict(self, m_sv):
        m_new_sv={}
        for s_id in m_sv:
            s_fields=s_id.split(global_values.SEPERATOR)
            s_new_id=global_values.SEPERATOR.join(s_fields[:-1])
            m_new_sv[s_new_id]=(s_id, m_sv[s_id])
        return m_new_sv

    ####
    def classify_inter_events(self, l_inter):#for inter-translocation
        pass

####
    #f_map_ratio: if
    #f_full_cover_ratio: if larger than this portion of sgmt is covered, then fully covered
    def parse_segmt_contig_algnmt(self, sf_algnmt, s_contig, f_map_ratio, f_full_cover_ratio, b_dup=False):
        t_rcd=None
        i_contig_len=len(s_contig)
        try:
            samfile = pysam.AlignmentFile(sf_algnmt)  # read in the sam file
            for algnmt in samfile.fetch():  # check each alignment, and find "left" and "right" flank
                if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                    continue
                if algnmt.is_unmapped == True:##unmapped
                    continue
                if algnmt.is_duplicate == True:##duplciate records
                    continue
                # as  the clipped long reads are from the alignment, so already corrected from rc.
                # Note, we are not sure about the assembler whether it generate a reverse ones
                b_rc = algnmt.is_reverse
                if (algnmt.mapping_quality < global_values.MINIMUM_DISC_MAPQ) and b_dup==False:
                    print "Low mapping quality"
                    continue
                # make sure most region is aligned
                l_cigar = algnmt.cigar

                i_lclip_len = 0  # left-clip length
                i_rclip_len = 0
                if l_cigar[0][0] == 5 or l_cigar[0][0] == 4:
                    i_lclip_len = l_cigar[0][1]
                if l_cigar[-1][0] == 5 or l_cigar[-1][0] == 4:
                    i_rclip_len = l_cigar[-1][1]

                i_clip_total_len = i_lclip_len + i_rclip_len  # left flank total clip length
                i_seq_len = len(algnmt.query_sequence)
                i_max_clip = int(i_seq_len * (1 - f_map_ratio))  # at most half is clipped out
                if i_clip_total_len>i_max_clip:
                    continue
                i_map_start = algnmt.reference_start
                i_map_end=i_map_start+(i_seq_len-i_clip_total_len)

                # if float(i_seq_len-i_clip_total_len)/float(i_contig_len)<f_full_cover_ratio:
                #     continue
                t_rcd=(b_rc, i_map_start, i_map_end, len(s_contig))
            samfile.close()
        except ValueError:
            print sf_algnmt, "is empty"
        return t_rcd
####

    # parse the alignment for duplication: tandem duplication or TE promoted duplication
    # f_full_cover_ratio: if larger than this portion of sgmt is covered, then fully covered
    def parse_segmt_contig_algnmt_for_dup(self, sf_algnmt, s_contig, f_map_ratio, f_full_cover_ratio):
        t_rcd = None
        i_contig_len = len(s_contig)
        try:
            samfile = pysam.AlignmentFile(sf_algnmt)  # read in the sam file
            for algnmt in samfile.fetch():  # check each alignment, and find "left" and "right" flank
                if algnmt.is_secondary or algnmt.is_supplementary:  # filter out secondary and supplementary
                    continue
                if algnmt.is_unmapped == True:  ##unmapped
                    continue
                if algnmt.is_duplicate == True:  ##duplciate records
                    continue
                # as  the clipped long reads are from the alignment, so already corrected from rc.
                # Note, we are not sure about the assembler whether it generate a reverse ones
                b_rc = algnmt.is_reverse
                if algnmt.mapping_quality < global_values.MINIMUM_DISC_MAPQ:
                    print "Low mapping quality"
                    continue
                # make sure most region is aligned
                l_cigar = algnmt.cigar

                i_lclip_len = 0  # left-clip length
                i_rclip_len = 0
                if l_cigar[0][0] == 5 or l_cigar[0][0] == 4:
                    i_lclip_len = l_cigar[0][1]
                if l_cigar[-1][0] == 5 or l_cigar[-1][0] == 4:
                    i_rclip_len = l_cigar[-1][1]

                i_clip_total_len = i_lclip_len + i_rclip_len  # left flank total clip length
                i_seq_len = len(algnmt.query_sequence)
                i_max_clip = int(i_seq_len * (1 - f_map_ratio))  # at most half is clipped out
                if i_clip_total_len > i_max_clip:
                    continue
                i_map_start = algnmt.reference_start
                i_map_end = i_map_start + (i_seq_len - i_clip_total_len)

                # if float(i_seq_len-i_clip_total_len)/float(i_contig_len)<f_full_cover_ratio:
                #     continue
                t_rcd = (b_rc, i_map_start, i_map_end, len(s_contig))
            samfile.close()
        except ValueError:
            print sf_algnmt, "is empty"
        return t_rcd
####
    ####
    def dump_ref_seq_to_file(self, l_intra, l_seq, sf_wfolder):
        m_ref_seq={}
        m_contig = {}
        m_contig_seq={}
        i_idx=0#
        for (s_chrm, s_lpos, s_rpos, s_l_seq, s_r_seq) in l_intra:
            s_id=s_chrm+global_values.SEPERATOR+s_lpos+global_values.SEPERATOR+s_rpos
            s_seq=l_seq[i_idx]
            i_idx+=1
            sf_fa=sf_wfolder+s_id+".ref_seq.fa"
            with open(sf_fa,"w") as fout_fa:
                fout_fa.write(">"+s_id+"_ref_sgmt\n")
                fout_fa.write(s_seq+"\n")
            m_ref_seq[s_id]=sf_fa

            sf_fa2 = sf_wfolder + s_id + ".contig.fa"
            s_tmp_seq=s_l_seq#
            with open(sf_fa2, "w") as fout_fa2:
                fout_fa2.write(">" + s_id + "\n")
                if len(s_l_seq) > len(s_r_seq):
                    fout_fa2.write(s_l_seq + "\n")
                else:
                    fout_fa2.write(s_r_seq + "\n")
                    s_tmp_seq=s_r_seq
            m_contig[s_id] = sf_fa2
            m_contig_seq[s_id] = s_tmp_seq
        return m_ref_seq, m_contig, m_contig_seq


    ####
    def define_inter_intra_chrm_events(self, m_candidates):
        l_inter_chrm=[]
        l_intra_chrm=[]
        for s_id in m_candidates:
            (s_lchrm, s_lpos, s_rchrm, s_rpos, s_l_seq, s_r_seq)=m_candidates[s_id]
            if s_lchrm==s_rchrm:
                l_intra_chrm.append((s_lchrm, s_lpos, s_rpos, s_l_seq, s_r_seq))
            else:
                l_inter_chrm.append((s_lchrm, s_lpos, s_rchrm, s_rpos, s_l_seq, s_r_seq))
        return l_inter_chrm, l_intra_chrm

####

####
####
    ####Get the duplication-like breakpoints, which:
    ####1. right clipped reads have breakpoints on the right
    ####2. left clipped reads have breakpoints on the left
    # def _define_dup_brkpnts(self, sf_lbrkpnts, sf_rbrkpnts):
    #     m_dup_lbrkpnts={}
    #     m_dup_rbrkpnts={}
    #     with open(sf_lbrkpnts) as fin_lbrk:
    #         for line in fin_lbrk:
    #             fields=line.split()
    #             n_lclip=int(fields[3])
    #             n_rclip=int(fields[4])
    #
    #             if n_lclip>n_rclip:
    #                 s_id="%s%s%s" % (fields[0], global_values.SEPERATOR, fields[1])
    #                 m_dup_lbrkpnts[s_id]=int(fields[2])
    #
    #     with open(sf_rbrkpnts) as fin_rbrk:
    #         for line in fin_rbrk:
    #             fields=line.split()
    #             n_lclip = int(fields[3])
    #             n_rclip = int(fields[4])
    #
    #             if n_rclip>n_lclip:
    #                 s_id = "%s%s%s" % (fields[0], global_values.SEPERATOR, fields[1])
    #                 m_dup_rbrkpnts[s_id] = int(fields[2])
    #     return m_dup_lbrkpnts, m_dup_rbrkpnts
####
####