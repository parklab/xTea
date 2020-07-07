##3/08/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

####More filters:
#1. mapq< 20
#2. aligned region < 1k
#3. standard derivation of the cluster breakpoints (what's the cutoff)?

#TD list:
#simple deletions that saved in cigar are saved, but not called out. At line ~185.
####

import pysam
from multiprocessing import Pool
from x_alignments import *
from x_intermediate_sites import *

def unwrap_self_collect_clip_pos_lrd(arg, **kwarg):
    return LRDBreakpoints.collect_breakpoints_by_chrm(*arg, **kwarg)

####
class LRDBreakpoints():
    def __init__(self, n_jobs, sf_bam="", sf_ref="", s_working_folder=""):
        self.n_jobs = n_jobs
        self.sf_bam = sf_bam
        self.sf_reference = sf_ref
        self.working_folder = s_working_folder
        self._LEFT_BRKPNT=".left_breakpoint"
        self._RIGHT_BRKPNT=".right_breakpoint"

    #for each bam we have a separate output
    def collect_breakpoints(self, i_peak_win, i_clip_cutoff, i_max_cutoff, std_dev_cutoff, sf_brkpnts,
                            sf_lbrkpnts, sf_rbrkpnts):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        references = samfile.references

        xchrom = XChromosome()
        l_chrm_records = []
        for chrm in references:
            if xchrom.is_decoy_contig_chrms(chrm) == True:  ###decoy sequnces and contigs are not considered
                continue
            sf_tmp_out=self.working_folder + "_" + chrm + "_tmp_candidate_breakpoints.txt"
            l_chrm_records.append((chrm, self.sf_bam, i_peak_win, sf_tmp_out))
            #for test only
            # tmp_rcd=(chrm, self.sf_bam, i_peak_win, sf_tmp_out)
            # self.collect_breakpoints_by_chrm(tmp_rcd)
        samfile.close()
####
        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_clip_pos_lrd, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        self._save_brkpnts_info(sf_brkpnts, l_chrm_records, i_clip_cutoff, i_max_cutoff, std_dev_cutoff)
        sf_lbrk_tmp=sf_lbrkpnts+".tmp"
        self._save_brkpnts_info(sf_lbrk_tmp, l_chrm_records, i_clip_cutoff, i_max_cutoff, std_dev_cutoff,
                                self._LEFT_BRKPNT, "a")
        sf_rbrk_tmp=sf_rbrkpnts+".tmp"
        self._save_brkpnts_info(sf_rbrk_tmp, l_chrm_records, i_clip_cutoff, i_max_cutoff, std_dev_cutoff,
                                self._RIGHT_BRKPNT, "a")

        #filter out left and right sites that are too close to each other
        self.filter_left_right_brkpnts_by_distance(sf_lbrk_tmp, sf_rbrk_tmp, global_values.LRD_MIN_INTERNAL_DEL,
                                                   sf_lbrkpnts, sf_rbrkpnts)

####
    ####
    def _save_brkpnts_info(self, sf_brkpnts, l_chrm_records, i_clip_cutoff, i_max_cutoff, std_dev_cutoff,
                           s_suffix="", s_w_flag="w"):
        #now merge the chromosomes and also filter with cutoff
        with open(sf_brkpnts, "w") as fout_brkpnts, open(sf_brkpnts+".1", s_w_flag) as fout_brkpnts1:
            for rcd in l_chrm_records:
                chrm =rcd[0]
                sf_tmp_out = self.working_folder + "_"+ chrm + "_tmp_candidate_breakpoints.txt"+s_suffix
                if os.path.isfile(sf_tmp_out)==False:
                    continue
                x_intermediate_sites = XIntemediateSites()
                m_chrm_sites=x_intermediate_sites.load_in_candidate_list_str_version(sf_tmp_out)
                for tmp_chrm in m_chrm_sites:
                    for tmp_pos in m_chrm_sites[tmp_chrm]:
                        n_clip=int(m_chrm_sites[tmp_chrm][tmp_pos][0])
                        n_lclip=int(m_chrm_sites[tmp_chrm][tmp_pos][1])
                        n_rclip=int(m_chrm_sites[tmp_chrm][tmp_pos][2])
                        n_contained=int(m_chrm_sites[tmp_chrm][tmp_pos][3])
                        n_focal_clip_contain=int(m_chrm_sites[tmp_chrm][tmp_pos][4])
                        f_std_dev=float(m_chrm_sites[tmp_chrm][tmp_pos][5])
                        b_save, b_lclip, b_rclip = self.filter_candidate_by_cutoff(m_chrm_sites, tmp_chrm, tmp_pos,
                                                                                   i_clip_cutoff, i_max_cutoff,
                                                                                   std_dev_cutoff)####
                        s_info="{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(n_clip, n_lclip, n_rclip, n_contained,
                                                                     n_focal_clip_contain, f_std_dev)
                        if s_suffix==self._LEFT_BRKPNT and b_rclip==True:
                            fout_brkpnts.write(tmp_chrm + "\t" + str(tmp_pos) + "\t" + s_info + "\n")
                        elif s_suffix==self._RIGHT_BRKPNT and b_lclip==True:
                            fout_brkpnts.write(tmp_chrm + "\t" + str(tmp_pos) + "\t" + s_info + "\n")
                        elif s_suffix=="" and b_save == True:
                            fout_brkpnts.write(tmp_chrm+"\t"+str(tmp_pos)+"\t"+s_info+"\n")
                        fout_brkpnts1.write(tmp_chrm+"\t"+str(tmp_pos)+"\t"+s_info+"\n") ####save all the sites

####
    ####Filter 1: If the n_singal > cutoff, and std_derivation < cutoff2, then pass
    ####Or pass filter 2: the focal region > cutoff singals * 85%, then pass
    def filter_candidate_by_cutoff(self, m_chrm_sites, tmp_chrm, tmp_pos, i_clip_cutoff, i_max_cutoff, std_dev_cutoff):
        b_save=False
        b_lclip=False
        b_rclip=False
        n_total_clip = int(m_chrm_sites[tmp_chrm][tmp_pos][0])
        n_focal_clip=int(m_chrm_sites[tmp_chrm][tmp_pos][4])
        f_std_dev=float(m_chrm_sites[tmp_chrm][tmp_pos][5])

        #focal region has enough support
        if n_focal_clip >= int(i_clip_cutoff*global_values.LRD_BRKPNT_FOCAL_CLIP_RATIO):
            b_save=True

        #total have enough support, and standard derivation is small enough
        if n_total_clip >= i_clip_cutoff and f_std_dev<=std_dev_cutoff:
            b_save = True

        if n_total_clip > i_max_cutoff: #too many clipped reads
            b_save=False
        if int(f_std_dev)>global_values.LRD_BRKPNT_MAXIMUM_STD:#over the allowed maximum standard derivation
            b_save=False

        n_lclip = int(m_chrm_sites[tmp_chrm][tmp_pos][1])
        n_rclip = int(m_chrm_sites[tmp_chrm][tmp_pos][2])
        n_contained = int(m_chrm_sites[tmp_chrm][tmp_pos][3])
        if n_lclip == 0 and n_contained == 0:  # if only have right clipped, then skip
            if b_save==True:
                b_rclip = True
            b_save = False

        if n_rclip == 0 and n_contained == 0:  # if only have left clipped, then skip
            if b_save==True:
                b_lclip = True
            b_save = False
        return b_save, b_lclip, b_rclip
####
####
    #collect the breakpoints (clip breakpoints and contained insertion, deletion breakpoints) for each chrom
    def collect_breakpoints_by_chrm(self, record):
        chrm = record[0]
        sf_bam = record[1]
        i_peak_win=int(record[2])
        sf_out = record[3]

        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        m_breakpoint_pos = {}
        m_lbreakpoint_pos = {}#left breakpoints
        m_rbreakpoint_pos = {}  # right breakpoints
        for algnmt in samfile.fetch(chrm):  ##fetch reads mapped to "chrm"
            ##here need to skip the secondary and supplementary alignments?
            #if algnmt.is_secondary or algnmt.is_supplementary:
            if algnmt.is_secondary:#skip supplementary, but keep the secondary alignment
                continue
            if algnmt.is_duplicate == True:  ##duplciate
                continue
            if algnmt.is_unmapped == True:  ##unmapped
                continue
            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue
            if algnmt.mapping_quality < global_values.LRD_MIN_MAPQ:#by default this is set to 20
                continue

            #query_name = algnmt.query_name ####note, this one may have space in the middle
            #query_seq = algnmt.query_sequence
            #query_quality = algnmt.query_qualities  ##this is different from the one saved in the fastq/sam, \
            # no offset 33 to subtract
            map_pos = algnmt.reference_start
            # if query_seq == None:
            #     continue
            #seq_lenth = len(query_seq)
            lclip_len = 0
            rclip_len = 0
####
            # now check the contained ones:
            l_contain_ins_pos, l_ldel_pos, l_rdel_pos, rclip_pos, n_mapped = \
                self._get_contained_ins_del_rclip_pos_mapped_bases_from_cigar(map_pos, l_cigar,
                                                                              global_values.LRD_MIN_INTERNAL_DEL)
            if n_mapped<global_values.LRD_MIN_MAP_LEN:#require at least the mapped region is 1k (by default)
                continue
####
            for pos in l_contain_ins_pos:
                self._save_pos_to_dict(chrm, pos, m_breakpoint_pos)

###################################
#this is the simple deletions that saved in the cigar!!!! need special processing!!!!!
# for pos in l_ldel_pos:#for left breakpoint (from contained deletion), view as right clip
#     self._save_pos_to_dict(chrm, pos, m_lbreakpoint_pos, 1)
# for pos in l_rdel_pos:#for right breakpoint (from contained deletion), view as left clip
#     self._save_pos_to_dict(chrm, pos, m_rbreakpoint_pos, 0)#
###################################

            if l_cigar[0][0] == 4:#left clipped, right breakpoints
                lclip_len = l_cigar[0][1]
                clip_pos = map_pos
                if lclip_len > global_values.LRD_MIN_CLIP_LTH:
                    self._save_pos_to_dict(chrm, clip_pos, m_breakpoint_pos, 0)
                    self._save_pos_to_dict(chrm, clip_pos, m_rbreakpoint_pos, 0)

            if l_cigar[-1][0] == 4:  #right clipped, left breakpoints
                rclip_len = l_cigar[-1][1]
                if rclip_len > global_values.LRD_MIN_CLIP_LTH:
                    clip_pos = rclip_pos #this is calculated from cigar
                    if rclip_pos!=-1:
                        self._save_pos_to_dict(chrm, clip_pos, m_breakpoint_pos, 1)
                        self._save_pos_to_dict(chrm, clip_pos, m_lbreakpoint_pos, 1)
        samfile.close()
####
        x_intermediate_sites = XIntemediateSites()

#Here add a filtering step to remove those "random" clip reads????
#as these reads will bridge the clusters (especially for raw pacbio/nanopore reads)
#because the alignment is not accurate at the breakpoints!!!!
        b_save=True
        m_peak, m_tmp_brkpnts_cluster=x_intermediate_sites.call_peak_candidate_sites_lrd(m_breakpoint_pos,
                                                                                         i_peak_win, b_save)
        #for left-breakpoints only
        m_lpeak, m_tmp_lbrkpnts_cluster = x_intermediate_sites.call_peak_candidate_sites_lrd(m_lbreakpoint_pos,
                                                                                           i_peak_win, b_save)
        # for right-breakpoints only
        m_rpeak, m_tmp_rbrkpnts_cluster = x_intermediate_sites.call_peak_candidate_sites_lrd(m_rbreakpoint_pos,
                                                                                             i_peak_win, b_save)

####for temporary usage only:#######################
        sf_tmp_out=sf_out+".tmp_cluster"
        sf_tmp_out_std=sf_out+".tmp_cluster_std_dev"
        with open(sf_tmp_out,"w") as fout_tmp, open(sf_tmp_out_std, "w") as fout_std_dev:
            for tmp_chrm in m_tmp_brkpnts_cluster:
                for tmp_pos in m_tmp_brkpnts_cluster[tmp_chrm]:
                    l_tmp_pos=m_tmp_brkpnts_cluster[tmp_chrm][tmp_pos]
                    f_tmp_std=x_intermediate_sites.calc_std_derivation(l_tmp_pos)
                    s_tmp_pos=""
                    for i_pos in l_tmp_pos:
                        s_tmp_pos+="{0}\t".format(i_pos)
                    fout_tmp.write(tmp_chrm+"\t"+str(tmp_pos)+"\t"+s_tmp_pos+"\n")
                    fout_std_dev.write(tmp_chrm + "\t" + str(tmp_pos) + "\t" + s_tmp_pos + "\t" + str(f_tmp_std) + "\n")

        sf_tmp_out = sf_out + ".left_tmp_cluster"
        sf_tmp_out_std = sf_out + ".left_tmp_cluster_std_dev"
        with open(sf_tmp_out, "w") as fout_tmp, open(sf_tmp_out_std, "w") as fout_std_dev:
            for tmp_chrm in m_tmp_lbrkpnts_cluster:
                for tmp_pos in m_tmp_lbrkpnts_cluster[tmp_chrm]:
                    l_tmp_pos = m_tmp_lbrkpnts_cluster[tmp_chrm][tmp_pos]
                    f_tmp_std = x_intermediate_sites.calc_std_derivation(l_tmp_pos)
                    s_tmp_pos = ""
                    for i_pos in l_tmp_pos:
                        s_tmp_pos += "{0}\t".format(i_pos)
                    fout_tmp.write(tmp_chrm + "\t" + str(tmp_pos) + "\t" + s_tmp_pos + "\n")
                    fout_std_dev.write(tmp_chrm + "\t" + str(tmp_pos) + "\t" + s_tmp_pos + "\t" + str(f_tmp_std) + "\n")

        sf_tmp_out = sf_out + ".right_tmp_cluster"
        sf_tmp_out_std = sf_out + ".right_tmp_cluster_std_dev"
        with open(sf_tmp_out, "w") as fout_tmp, open(sf_tmp_out_std, "w") as fout_std_dev:
            for tmp_chrm in m_tmp_rbrkpnts_cluster:
                for tmp_pos in m_tmp_rbrkpnts_cluster[tmp_chrm]:
                    l_tmp_pos = m_tmp_rbrkpnts_cluster[tmp_chrm][tmp_pos]
                    f_tmp_std = x_intermediate_sites.calc_std_derivation(l_tmp_pos)
                    s_tmp_pos = ""
                    for i_pos in l_tmp_pos:
                        s_tmp_pos += "{0}\t".format(i_pos)
                    fout_tmp.write(tmp_chrm + "\t" + str(tmp_pos) + "\t" + s_tmp_pos + "\n")
                    fout_std_dev.write(tmp_chrm + "\t" + str(tmp_pos) + "\t" + s_tmp_pos + "\t" + str(f_tmp_std) + "\n")
################################################
####
        ####output the candidate sites
        x_intermediate_sites.output_candidate_sites(m_peak, sf_out)
        x_intermediate_sites.output_candidate_sites(m_lpeak, sf_out + self._LEFT_BRKPNT)
        x_intermediate_sites.output_candidate_sites(m_rpeak, sf_out + self._RIGHT_BRKPNT)
####
####
    #if two breakpoints from left and right dict specifically are close to each other, then they should be filtered out
    def filter_left_right_brkpnts_by_distance(self, sf_lbrk_tmp, sf_rbrk_tmp, i_dist, sf_lbrkpnts, sf_rbrkpnts):
        x_intermediate_sites=XIntemediateSites()
        m_lpeak=x_intermediate_sites.load_in_candidate_list_str_version(sf_lbrk_tmp)
        m_rpeak=x_intermediate_sites.load_in_candidate_list_str_version(sf_rbrk_tmp)
        m_new_lpeak={}
        m_new_rpeak={}
        m_l_need_rm={}
        m_r_need_rm={}
        for chrm in m_lpeak:
            if chrm not in m_rpeak:
                continue
            for pos in m_lpeak[chrm]:
                i_start=pos-i_dist
                i_end=pos+i_dist
                b_hit=False
                for i_tmp in range(i_start,i_end):
                    if i_tmp in m_rpeak[chrm]:
                        b_hit=True
                        if chrm not in m_r_need_rm:
                            m_r_need_rm[chrm]={}
                        m_r_need_rm[chrm][i_tmp]=1
                        continue
                if b_hit==True:
                    if chrm not in m_l_need_rm:
                        m_l_need_rm[chrm]={}
                    m_l_need_rm[chrm][pos]=1
        for chrm in m_lpeak:
            for pos in m_lpeak[chrm]:
                if (chrm in m_l_need_rm) and (pos in m_l_need_rm[chrm]):
                    continue
                if chrm not in m_new_lpeak:
                    m_new_lpeak[chrm]={}
                m_new_lpeak[chrm][pos]=m_lpeak[chrm][pos]
        for chrm in m_rpeak:
            for pos in m_rpeak[chrm]:
                if (chrm in m_r_need_rm) and (pos in m_r_need_rm[chrm]):
                    continue
                if chrm not in m_new_rpeak:
                    m_new_rpeak[chrm] = {}
                m_new_rpeak[chrm][pos] = m_rpeak[chrm][pos]
        x_intermediate_sites.output_candidate_sites(m_new_lpeak, sf_lbrkpnts)
        x_intermediate_sites.output_candidate_sites(m_new_rpeak, sf_rbrkpnts)

    ####save position to dictionary
    def _save_pos_to_dict(self, chrm, pos, m_pos, iflag=2):
        if chrm not in m_pos:
            m_pos[chrm] = {}
        if pos not in m_pos[chrm]:
            m_pos[chrm][pos] = []
            if iflag==0:#left clip case
                m_pos[chrm][pos].append([1, 0, 0])
            elif iflag==1:#right clip case
                m_pos[chrm][pos].append([0, 1, 0])
            else:#contained case
                m_pos[chrm][pos].append([0, 0, 1])
        else:
            if iflag == 0:
                m_pos[chrm][pos][0][0] += 1
            elif iflag==1:
                m_pos[chrm][pos][0][1] += 1
            else:
                m_pos[chrm][pos][0][2] += 1

    ####get contained and also rclip position from cigar
    def _get_contained_ins_del_rclip_pos_mapped_bases_from_cigar(self, map_pos, l_cigar, i_min_del_len):
        ref_pos = map_pos
        l_contain_pos = []
        l_ldel_pos=[]#left deletion position
        l_rdel_pos=[]#right deletion position
        i_rclip_pos = -1#right clip position
        n_mapped=0
        for (opn, lth) in l_cigar:
            if opn == 0:  # alignment match (M)
                ref_pos += lth
                n_mapped+=lth #accumulate the mapped length
            elif opn == 2 or opn == 3:  # deletion (D) or skipped region (N)
                b_qualified_del=False
                if lth>i_min_del_len:
                    b_qualified_del=True
                if b_qualified_del==True:
                    l_ldel_pos.append(ref_pos)#deletion has two breakpoints
                ref_pos += lth
                if b_qualified_del == True:
                    l_rdel_pos.append(ref_pos)#the other breakpoint
            elif opn == 5 or opn == 6:  # hard-clip (H) or padding (P)
                continue
            elif opn == 7 or opn == 8:  # sequence match (=) or sequence mismatch (X)
                n_mapped+=lth
                ref_pos += lth
            elif opn == 1:  # insertion
                if lth >= global_values.LRD_MIN_INS_LTH:
                    l_contain_pos.append(ref_pos)
            elif opn == 4 and len(l_cigar) > 0 and l_cigar[-1][0] == 4:#right clip
                i_rclip_pos = ref_pos
        return l_contain_pos, l_ldel_pos, l_rdel_pos, i_rclip_pos, n_mapped
####

    ####This is not used
    ####get contain insertion position from cigar
    def _get_contained_insertion_pos_from_cigar(self, map_pos, l_cigar):
        ref_pos = map_pos
        l_pos = []
        for (opn, lth) in l_cigar:
            if opn == 0:  # alignment match (M)
                ref_pos += lth
            elif opn == 2 or opn == 3:  # deletion (D) or skipped region (N)
                ref_pos += lth
            elif opn == 5 or opn == 6:  # hard-clip (H) or padding (P)
                continue
            elif opn == 7 or opn == 8:  # sequence match (=) or sequence mismatch (X)
                ref_pos += lth
            elif opn == 1:  # insertion
                if lth >= global_values.LRD_MIN_INS_LTH:
                    l_pos.append(ref_pos)
        return l_pos
####