import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
from x_alignments import *
from bwa_align import *
from x_reference import *
from x_consensus import *

BWA_PATH = "bwa"
BWA_REALIGN_CUTOFF = 9
SAMTOOLS_PATH = "samtools"
FLAG_LEFT_CLIP = "L"
FLAG_RIGHT_CLIP = "R"
SEPERATOR = '~'
NOT_TRANSDUCTION = "not_transduction"
MINIMAL_TRANSDUCT_MAPQ = 30
NEARBY_CLIP = 50
MAX_CLIP_CLIP_LEN = 6
MIN_POLYMORPHIC_SOURCE_DIST=500

class XTransduction():
    def __init__(self, working_folder, n_jobs, sf_reference):
        self.working_folder = working_folder
        self.n_jobs = n_jobs
        self.sf_reference=sf_reference

    ####call the transductions whose source from reference FL-L1
    def call_candidate_transductions(self, m_sites, sf_clip_algnmt, sf_disc_algnmt, sf_flank, i_flank_lenth, ndisc_cutoff):
        # extract the clipped and discordant reads of the selected sites
        # sf_clip_algnmt, sf_disc_algnmt, m_picked_sites, min_clip_lth,sf_disc_fa, sf_clip_fa
        # if clipped part is larger than this , then the read is saved
        # this is to filter out those reads already fully mapped to consensus
        min_clip_lth = 10
        sf_picked_clip = self.working_folder + "picked_candidate_sites_all_clip.fa"
        sf_picked_disc = self.working_folder + "picked_candidate_sites_all_disc.fa"
        self.collect_transduction_clip_disc_reads(sf_clip_algnmt, sf_disc_algnmt, m_sites, min_clip_lth,
                                                  sf_picked_clip, sf_picked_disc)

        # align the picked disc and clip parts to the FL-L1 flank regions
        bwa_align=BWAlign(BWA_PATH, BWA_REALIGN_CUTOFF, self.n_jobs)
        if os.path.isfile(sf_flank + ".sa") == False:
            cmd = "{0} index {1}".format(BWA_PATH, sf_flank)
            Popen(cmd, shell=True, stdout=PIPE).communicate()
        sf_picked_clip_sam = self.working_folder + "picked_clip.sam"
        bwa_align.realign_clipped_reads(sf_flank, sf_picked_clip, sf_picked_clip_sam)
        sf_picked_disc_sam = self.working_folder + "picked_disc.sam"
        bwa_align.realign_disc_reads(sf_flank, sf_picked_disc, sf_picked_disc_sam)

        # parse the disc and clip alignment to flank regions
        # m_lclip_transd format: m_ltransduct[ori_chrm][ori_insertion_pos][ref_mpos].append(s_transduct)
        #Here need to rule out those polymorphic full-length L1s, as their flank regions will be collected
        m_lclip_transd, m_rclip_transd, m_td_polyA = self.parse_transduction_clip_algnmt(sf_picked_clip_sam,
                                                                                         i_flank_lenth)
        m_disc_transd = self.parse_transduction_disc_algnmt(sf_picked_disc_sam, i_flank_lenth, ndisc_cutoff)

        return m_lclip_transd, m_rclip_transd, m_td_polyA, m_disc_transd


####
    def construct_novel_flanks(self, l_fl_polymorphic, sf_flank, i_flank_lenth, sf_new_flank):
        xref=XReference()
        sf_polymerphic_flanks=sf_new_flank+".polymerphic_only.fa"
        xref.gnrt_flank_regions_of_polymerphic_insertions(l_fl_polymorphic, i_flank_lenth, self.sf_reference,
                                                          sf_polymerphic_flanks)
        #merge the two flanks files, and index(bwa) the new file
        cmd="cat {0} {1} > {2}".format(sf_flank, sf_polymerphic_flanks, sf_new_flank)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    # if clipped, and the clipped part is larger than "min_clip_lth", then return True
    def _is_large_clipped(self, l_cigar, min_clip_lth):
        if len(l_cigar) < 1:  # wrong alignment
            return False
        if len(l_cigar) > 2:
            ####check the cigar
            ###if both clipped, and the clipped part is large, then skip
            b_left_clip = False
            i_left_clip_len = 0
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                b_left_clip = True
                i_left_clip_len = l_cigar[0][1]
            b_right_clip = False
            i_right_clip_len = 0
            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                b_right_clip = True
                i_right_clip_len = l_cigar[-1][1]

            ####
            if (b_left_clip == True and i_left_clip_len > min_clip_lth) or (
                            b_right_clip == True and i_right_clip_len > min_clip_lth):
                return True
        return False


    # for each site, collect the disc and clip reads not aligned to consensus
    def collect_transduction_clip_disc_reads(self, sf_clip_algnmt, sf_disc_algnmt, m_picked_sites, min_clip_lth,
                                             sf_clip_fa, sf_disc_fa):
        m_sites = {}
        for ins_chrm in m_picked_sites:
            for ins_pos in m_picked_sites[ins_chrm]:
                s_site = "{0}{1}{2}".format(ins_chrm, SEPERATOR, ins_pos)
                m_sites[s_site] = 1
        # collect the clipped reads
        # For clip reads, e.g.:4~11019868~R~0~11019545~98~0 [chrm, map-pos, clip_dir, rev-comp, ins_pos, cnt, sample]
        with open(sf_clip_fa, "w") as fout_clip:
            samfile = pysam.AlignmentFile(sf_clip_algnmt, "r", reference_filename=self.sf_reference)
            for algnmt in samfile.fetch():
                read_info = algnmt.query_name
                read_fields = read_info.split(SEPERATOR)
                rid_tmp = "{0}{1}{2}".format(read_fields[0], SEPERATOR, read_fields[-3])
                if rid_tmp not in m_sites:
                    continue
                read_seq = algnmt.query_sequence
                # if unmapped, then save the unmapped clipped reads
                b_save = False
                if algnmt.is_unmapped == True:
                    b_save = True
                else:
                    l_cigar = algnmt.cigar
                    b_save = self._is_large_clipped(l_cigar, min_clip_lth)

                if b_save == True:
                    fout_clip.write(">" + read_info + "\n")
                    fout_clip.write(read_seq + "\n")
            samfile.close()
        ####
        # collect the discordant reads
        # For disc reads, e.g.: read-id~0~181192671~3~181192845~0 [rid, is-first, pos, ins-chrm, ins-pos, sample]
        with open(sf_disc_fa, "w") as fout_disc:
            samfile = pysam.AlignmentFile(sf_disc_algnmt, "r", reference_filename=self.sf_reference)
            for algnmt in samfile.fetch():
                read_info = algnmt.query_name
                read_fields = read_info.split(SEPERATOR)
                rid_tmp = "{0}{1}{2}".format(read_fields[-3], SEPERATOR, read_fields[-2])
                if rid_tmp not in m_sites:
                    continue
                read_seq = algnmt.query_sequence
                # if unmapped, then save the unmapped clipped reads
                b_save = False
                if algnmt.is_unmapped == True:
                    b_save = True
                else:
                    l_cigar = algnmt.cigar
                    b_save = self._is_large_clipped(l_cigar, min_clip_lth)

                if b_save == True:
                    fout_disc.write(">" + read_info + "\n")
                    fout_disc.write(read_seq + "\n")
            samfile.close()

    # sf_clip_alignmt: the alignment on flank regions (reads not mapped (or clipped) to consensus)
    def parse_transduction_clip_algnmt(self, sf_clip_alignmt, flank_lth):
        samfile = pysam.AlignmentFile(sf_clip_alignmt, "r", reference_filename=self.sf_reference)
        # m_transduction = {}  # save the candidate transduction positions
        m_ltransduct = {}
        m_rtransduct = {}
        m_polyA = {}
        for algnmt in samfile.fetch():
            read_info = algnmt.query_name
            read_info_fields = read_info.split(SEPERATOR)

            ori_chrm = read_info_fields[0]
            ref_mpos = int(read_info_fields[1])  ##clip position on the reference
            ori_clip_flag = read_info_fields[2]
            ori_insertion_pos = int(read_info_fields[4])

            #####Hard code here#######################################################################################
            if abs(ref_mpos - ori_insertion_pos) > NEARBY_CLIP:  # only focus on the clipped reads nearby
                continue
            if algnmt.is_unmapped == True:  ####skip the unmapped reads
                continue

            b_left = True  # left clip
            if ori_clip_flag == FLAG_RIGHT_CLIP:
                b_left = False  # right clip

            bmapped_cutoff = 0.85
            l_cigar = algnmt.cigar
            b_clip_qualified_algned, n_map_bases = self.is_clipped_part_qualified_algnmt(l_cigar, bmapped_cutoff)
            if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
                continue

            clipped_seq = algnmt.query_sequence
            ##if it is mapped to the left flank regions
            # here make sure it is not poly-A
            # b_polya = self.contain_poly_A_T(clipped_seq, N_MIN_A_T)  # by default, at least 5A or 5T
            b_polya = self.is_consecutive_polyA_T(clipped_seq)
            if b_polya == True:
                if ori_chrm not in m_polyA:
                    m_polyA[ori_chrm] = {}
                if ori_insertion_pos not in m_polyA[ori_chrm]:
                    m_polyA[ori_chrm][ori_insertion_pos] = []
                    m_polyA[ori_chrm][ori_insertion_pos].append(0)
                    m_polyA[ori_chrm][ori_insertion_pos].append(0)
                if b_left == True:
                    m_polyA[ori_chrm][ori_insertion_pos][0] += 1
                else:
                    m_polyA[ori_chrm][ori_insertion_pos][1] += 1

            # calc the map position on the flank [convert to the reference genome]
            map_pos = algnmt.reference_start
            hit_flank_id = algnmt.reference_name  # e.g. chrY~3443451~3449565~L1HS~0L
            flank_id_fields = hit_flank_id.split(SEPERATOR)
            chrm_fl_L1 = flank_id_fields[0]
            source_start = int(flank_id_fields[1])
            source_end = int(flank_id_fields[2])
            sourc_rc = flank_id_fields[-1][0]

            if (ori_chrm==chrm_fl_L1) and (abs(ref_mpos-source_start)<MIN_POLYMORPHIC_SOURCE_DIST
                                           or abs(ref_mpos-source_end)<MIN_POLYMORPHIC_SOURCE_DIST):
                continue

            # here need to process the chrom, and keep it consistent with the final list
            if len(chrm_fl_L1) > 3 and chrm_fl_L1[:3] == "chr":  # contain
                if len(ori_chrm) <= 3 or ori_chrm[:3] != "chr":
                    chrm_fl_L1 = chrm_fl_L1[3:]  # trim the "chr"

            hit_ref_pos = int(flank_id_fields[2]) + map_pos
            if flank_id_fields[-1][-1] == "L":  # left-flank region
                hit_ref_pos = int(flank_id_fields[1]) - flank_lth + map_pos

            if b_left == True:  # for left-clipped reads
                if ori_chrm not in m_ltransduct:
                    m_ltransduct[ori_chrm] = {}
                if ori_insertion_pos not in m_ltransduct[ori_chrm]:
                    m_ltransduct[ori_chrm][ori_insertion_pos] = {}
                if ref_mpos not in m_ltransduct[ori_chrm][ori_insertion_pos]:
                    m_ltransduct[ori_chrm][ori_insertion_pos][ref_mpos] = []
                s_transduct = "{0}:{1}-{2}~{3}~{4}".format(chrm_fl_L1, source_start, source_end, sourc_rc, hit_ref_pos)
                m_ltransduct[ori_chrm][ori_insertion_pos][ref_mpos].append(s_transduct)
            else:
                if ori_chrm not in m_rtransduct:
                    m_rtransduct[ori_chrm] = {}
                if ori_insertion_pos not in m_rtransduct[ori_chrm]:
                    m_rtransduct[ori_chrm][ori_insertion_pos] = {}
                if ref_mpos not in m_rtransduct[ori_chrm][ori_insertion_pos]:
                    m_rtransduct[ori_chrm][ori_insertion_pos][ref_mpos] = []
                s_transduct = "{0}:{1}-{2}~{3}~{4}".format(chrm_fl_L1, source_start, source_end, sourc_rc, hit_ref_pos)
                m_rtransduct[ori_chrm][ori_insertion_pos][ref_mpos].append(s_transduct)
        samfile.close()
        return m_ltransduct, m_rtransduct, m_polyA

    ####Parse the alignment to the FL-L1 flank regions, and call out the transductions.
    def parse_transduction_disc_algnmt(self, sf_disc_alignmt, flank_lth, n_disc_cutoff):
        samfile = pysam.AlignmentFile(sf_disc_alignmt, "r", reference_filename=self.sf_reference)
        m_transduction = {}  # save the candidate transduction positions
        for algnmt in samfile.fetch():
            read_info = algnmt.query_name
            read_info_fields = read_info.split(SEPERATOR)
            anchor_map_pos = int(read_info_fields[-4])
            ins_chrm = read_info_fields[-3]
            ins_pos = int(read_info_fields[-2])
            sample_id = read_info_fields[-1]

            if algnmt.is_unmapped == True:  ####skip the unmapped reads
                continue
            if algnmt.mapping_quality <= MINIMAL_TRANSDUCT_MAPQ:
                continue
            l_cigar = algnmt.cigar
            bmapped_cutoff = 0.85
            b_clip_qualified_algned, n_map_bases = self.is_clipped_part_qualified_algnmt(l_cigar, bmapped_cutoff)
            if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
                continue

            map_pos = algnmt.reference_start
            hit_flank_id = algnmt.reference_name  # e.g. chrY~3443451~3449565~L1HSL
            flank_id_fields = hit_flank_id.split(SEPERATOR)
            chrm_fl_L1 = flank_id_fields[0]
            source_start = int(flank_id_fields[1])
            source_end = int(flank_id_fields[2])
            sourc_rc = flank_id_fields[-1][0]
            # here need to process the chrom, and keep it consistent with the final list
            if len(chrm_fl_L1) > 3 and chrm_fl_L1[:3] == "chr":  # contain
                if len(ins_chrm) <= 3 or ins_chrm[:3] != "chr":  # keep it same as the ins_chrm
                    chrm_fl_L1 = chrm_fl_L1[3:]  # trim the "chr"

            hit_ref_pos = int(flank_id_fields[2]) + map_pos
            if flank_id_fields[-1][-1] == "L":  # left-flank region
                hit_ref_pos = int(flank_id_fields[1]) - flank_lth + map_pos
            if ins_chrm not in m_transduction:
                m_transduction[ins_chrm] = {}
            if ins_pos not in m_transduction[ins_chrm]:
                m_transduction[ins_chrm][ins_pos] = {}

            s_source = "{0}:{1}-{2}~{3}".format(chrm_fl_L1, source_start, source_end, sourc_rc)

            ####this is to rule out those polymerphic Fl-L1 cases, which aligned to themselves' flank regions
            if (ins_chrm==chrm_fl_L1) and (abs(ins_pos-source_start)<200 or abs(ins_pos-source_end)<300):
                continue

            if s_source not in m_transduction[ins_chrm][ins_pos]:
                m_transduction[ins_chrm][ins_pos][s_source] = []
            m_transduction[ins_chrm][ins_pos][s_source].append((hit_ref_pos, anchor_map_pos))
        samfile.close()
        m_picked = {}
        # for each candidate, it may link to several sources, here select the one with largest # of disc reads support
        for ins_chrm in m_transduction:
            for ins_pos in m_transduction[ins_chrm]:
                max_cnt = 0
                max_source = ""
                for tmp_source in m_transduction[ins_chrm][ins_pos]:
                    tmp_cnt = len(m_transduction[ins_chrm][ins_pos][tmp_source])
                    if tmp_cnt > max_cnt:
                        max_cnt = tmp_cnt
                        max_source = tmp_source
                if max_cnt < n_disc_cutoff:  ####Here require at least "n_disc_cutoff" support
                    continue
                if ins_chrm not in m_picked:
                    m_picked[ins_chrm] = {}
                l_pos = []
                for tmp_record in m_transduction[ins_chrm][ins_pos][max_source]:
                    l_pos.append(tmp_record)  # each record is (hit_ref_pos, anchor_pos)
                m_picked[ins_chrm][ins_pos] = (max_source, l_pos)
        return m_picked

    # check the clipped part is qualified aligned or not
    def is_clipped_part_qualified_algnmt(self, l_cigar, ratio_cutoff):
        if len(l_cigar) < 1:  # wrong alignment
            return False, 0
        if len(l_cigar) > 2:
            ####check the cigar
            ###if both clipped, and the clipped part is large, then skip
            b_left_clip = False
            i_left_clip_len = 0
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                b_left_clip = True
                i_left_clip_len = l_cigar[0][1]
            b_right_clip = False
            i_right_clip_len = 0
            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                b_right_clip = True
                i_right_clip_len = l_cigar[-1][1]

            if b_left_clip == True and b_right_clip == True:
                if (i_left_clip_len > MAX_CLIP_CLIP_LEN) and (i_right_clip_len > MAX_CLIP_CLIP_LEN):
                    return False, 0

        ####for the alignment (of the clipped read), if the mapped part is smaller than the clipped part,
        ####then skip
        n_total = 0
        n_map = 0
        for (type, lenth) in l_cigar:
            if type == 0:
                n_map += lenth
            if type != 2:  # deletion is not added to the total length
                n_total += lenth

        if n_map < (n_total * ratio_cutoff):  ########################require at least 3/4 of the seq is mapped !!!!!!!!
            return False, 0
        return True, n_map

    def is_consecutive_polyA_T(self, seq):
        if ("AAAAA" in seq) or ("TTTTT" in seq) or ("AATAA" in seq) or ("TTATT" in seq):
            return True
        else:
            return False

    # call out the candidate transductions only use the disc info (m_disc_transdct)
    # here still consdier the left/right anchor, so count seperately for left and right disc reads
    def call_out_candidate_transduct(self, m_repsnt_pos, m_disc_transdct, icorcord, ratio):
        m_transduction = {}
        ck_cluster=ClusterChecker()
        for ins_chrm in m_repsnt_pos:
            for ins_pos in m_repsnt_pos[ins_chrm]:
                clip_pos = m_repsnt_pos[ins_chrm][ins_pos]
                b_l_trsdct = False
                b_r_trsdct = False
                trsdct_l = -1
                trsdct_r = -1
                n_disc_trsdct = 0
                if (ins_chrm in m_disc_transdct) and (ins_pos in m_disc_transdct[ins_chrm]):
                    trsdct_source, l_trsdct_pos = m_disc_transdct[ins_chrm][ins_pos]
                    l_l_trsdct_pos = []
                    l_r_trsdct_pos = []
                    for (hit_ref_pos, anchor_pos) in l_trsdct_pos:
                        if anchor_pos < clip_pos:
                            l_l_trsdct_pos.append(hit_ref_pos)
                        else:
                            l_r_trsdct_pos.append(hit_ref_pos)
                        if len(l_l_trsdct_pos) > len(l_r_trsdct_pos):
                            b_l_trsdct, trsdct_l, trsdct_r = ck_cluster._is_disc_cluster(l_l_trsdct_pos, icorcord, ratio)
                            n_disc_trsdct = len(l_l_trsdct_pos)
                        else:
                            b_r_trsdct, trsdct_l, trsdct_r = ck_cluster._is_disc_cluster(l_r_trsdct_pos, icorcord, ratio)
                            n_disc_trsdct = len(l_r_trsdct_pos)
                            ####
                    if b_l_trsdct == True or b_r_trsdct == True:
                        if ins_chrm not in m_transduction:
                            m_transduction[ins_chrm] = {}
                        trsdct_info = "{0}~{1}-{2}".format(trsdct_source, trsdct_l, trsdct_r)
                        m_transduction[ins_chrm][ins_pos] = (b_l_trsdct, b_r_trsdct, n_disc_trsdct, trsdct_info)
        return m_transduction
####
    # In this version, the transduction candidates are called out only by discordant reads
    # then, in this function, for each candidate, extract the related information
    def _extract_transduct_info(self, m_list, m_transduct, m_lclip_transduct, m_rclip_transduct, m_td_polyA, BIN_SIZE):
        m_info = {}
        ckcluster = ClusterChecker()
        for ins_chrm in m_list:
            for ins_pos in m_list[ins_chrm]:
                if (ins_chrm not in m_transduct) or (ins_pos not in m_transduct[ins_chrm]):
                    continue
                record = m_transduct[ins_chrm][ins_pos]
                b_l_trsdct = record[0]
                n_disc_trsdct = record[2]
                trsdct_info = record[3]

                # for the clip information
                n_clip_trsdct = 0
                peak_ref_clip_pos = -1
                if b_l_trsdct == True:  # left-transduction
                    if (ins_chrm in m_lclip_transduct) and (ins_pos in m_lclip_transduct[ins_chrm]):
                        # m_lclip_pos = m_lclip_transduct[ins_chrm][ins_pos]
                        l1_pos, l1_nbr, l1_acm, l2p, l2n, l2a = ckcluster.find_first_second_peak(m_lclip_transduct,
                                                                                            ins_chrm,
                                                                                            ins_pos, BIN_SIZE)
                        n_clip_trsdct = l1_acm
                        peak_ref_clip_pos = l1_pos
                else:
                    if (ins_chrm in m_rclip_transduct) and (ins_pos in m_rclip_transduct[ins_chrm]):
                        # m_rclip_pos = m_rclip_transduct[ins_chrm][ins_pos]
                        l1_pos, l1_nbr, l1_acm, l2p, l2n, l2a = ckcluster.find_first_second_peak(m_rclip_transduct,
                                                                                            ins_chrm,
                                                                                            ins_pos, BIN_SIZE)
                        n_clip_trsdct = l1_acm
                        peak_ref_clip_pos = l1_pos

                # for the poly-A
                # number of polyA reads
                n_polyA = 0
                if (ins_chrm in m_td_polyA) and (ins_pos in m_td_polyA[ins_chrm]):
                    if b_l_trsdct == True:
                        n_polyA = m_td_polyA[ins_chrm][ins_pos][0]
                    else:
                        n_polyA = m_td_polyA[ins_chrm][ins_pos][1]

                # save the information
                if ins_chrm not in m_info:
                    m_info[ins_chrm] = {}
                rcd_info = (b_l_trsdct, n_clip_trsdct, peak_ref_clip_pos, n_disc_trsdct, n_polyA)
                m_info[ins_chrm][ins_pos] = (trsdct_info, rcd_info)
        return m_info
