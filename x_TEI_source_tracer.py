##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import sys
import networkx as nx
import pysam
from x_annotation import *
from x_alignments import *

S_DELIM = "~"


class SourceTracer():
    def __init__(self, sf_vcf, sf_bam_coordinate, sf_bam_barcode, sf_annotation, n_jobs, s_working_folder, sf_ref):
        self.sf_vcf = sf_vcf
        self.sf_bam = sf_bam_coordinate
        self.sf_bam_barcode = sf_bam_barcode
        self.m_sites = {}
        self.n_jobs = n_jobs
        self.b_with_chr = False
        self.s_working_folder = s_working_folder
        if len(self.s_working_folder) > 0 and self.s_working_folder[-1] != "/":
            self.s_working_folder += "/"
        self.out_header = None  ##this is set in the "_is_chrm_contain_chr" function
        self.xannotation = XAnnotation(sf_annotation)
        self.sf_reference=sf_ref

    def load_files(self):
        self._is_chrm_contain_chr()
        self.xannotation.set_with_chr(self.b_with_chr)
        self.xannotation.load_rmsk_annotation_L1()
        self.xannotation.index_rmsk_annotation()
        self.load_sites_from_vcf()

    # load in the MEI sites from vcf file
    def load_sites_from_vcf(self):
        with open(self.sf_vcf) as fin_vcf:
            for line in fin_vcf:
                if line[0] == "#":
                    continue
                fields = line.split()
                if fields[-1] == "0|0":
                    continue
                tmp_chrm = fields[0]
                chrm = self._process_chrm_name(tmp_chrm)
                pos = int(fields[1])
                if chrm not in self.m_sites:
                    self.m_sites[chrm] = set()
                self.m_sites[chrm].add(pos)

    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, chrm):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        if self.b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif self.b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif self.b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    # if the chromosome name in format: chr1, then return true,
    # else if in format: 1, then return false
    def _is_chrm_contain_chr(self):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb", reference_filename=self.sf_reference)
        self.out_header = samfile.header  # it's a dictionary
        samfile.close()

        l_chrms = self.out_header['SQ']
        m_chrms = {}
        for record in l_chrms:
            chrm_id = record['SN']
            m_chrms[chrm_id] = 1
        if "chr1" in m_chrms:
            self.b_with_chr = True
        else:
            self.b_with_chr = False

    ##for this network, use the total number of connections as the weight for each edge
    # m_sites_algnmnts in structure of {barcode, [alignments]}
    def _construct_simple_network_for_site(self, site_chrm, site_pos, m_sites_algnmnts, sf_algnmt):
        m_network = {}
        s_algnmt_fa = sf_algnmt + ".fasta"
        sf_fa = open(s_algnmt_fa, "w")
        for barcode in m_sites_algnmnts:
            for alnmt in m_sites_algnmnts[barcode]:
                # f_sam.write(alnmt)  ###write the alignment to file
                # write to fastq file
                xalmt = XAlignment()
                s_fa = xalmt._cvt_algnmt_to_fa(alnmt)
                sf_fa.write(s_fa)

                ##get the map position and check whether mapped to a repeat region
                chrm = alnmt.get_tag('ch')
                pos = int(alnmt.get_tag('cp'))

                mate_chrm = ""
                mate_pos = 0
                if alnmt.mate_is_unmapped == False:
                    mate_chrm = alnmt.get_tag('nh')
                    mate_pos = int(alnmt.next_reference_start)
                ##check the mapping position, whether belong to L1 region
                ##collect the positions

                bhit, region_pos = self.xannotation.is_within_repeat_region(chrm, pos)
                bmate_hit, mate_region_pos = self.xannotation.is_within_repeat_region(mate_chrm, mate_pos)

                if bhit == True:
                    s_hit_id = site_chrm + S_DELIM + str(site_pos) + S_DELIM + chrm + S_DELIM + str(region_pos)
                    if s_hit_id in m_network:
                        m_network[s_hit_id] += 1
                    else:
                        m_network[s_hit_id] = 1

                if bmate_hit == True:
                    s_mate_hit_id = site_chrm + S_DELIM + str(site_pos) + S_DELIM + mate_chrm + S_DELIM + str(
                        mate_region_pos)
                    if s_mate_hit_id in m_network:
                        m_network[s_mate_hit_id] += 1
                    else:
                        m_network[s_mate_hit_id] = 1

                ##for a pair, the relationship is only saved for one time
                if bhit == True and bmate_hit == True and alnmt.is_read1:
                    s_connect_id = chrm + S_DELIM + str(region_pos) + S_DELIM + mate_chrm + S_DELIM + str(
                        mate_region_pos)
                    if s_connect_id in m_network:
                        m_network[s_connect_id] += 1
                    else:
                        m_network[s_connect_id] = 1
        sf_fa.close()

        # construct the graph
        G = nx.MultiDiGraph()
        for s_connect_id in m_network:
            iweight = m_network[s_connect_id]
            fields = s_connect_id.split(S_DELIM)
            left_node = fields[0] + ":" + fields[1]
            right_node = fields[2] + ":" + fields[3]
            G.add_edge(left_node, right_node, weight=iweight)
        return G

    # initialize the bins for a given TE region
    def _init_bins_for_region(self, region_lth, bin_size):
        l_hit_cnt = []
        n_bin = region_lth / bin_size
        if region_lth % bin_size != 0:
            n_bin += 1
        for i in range(n_bin):
            l_hit_cnt.append(0)
        return l_hit_cnt

    # increase the number of hit of the related bin the position falls in
    def _increase_hit_cnt_at_pos(self, pos, start_pos, bin_size, l_hit_cnt):
        lth = pos - start_pos
        if lth < 0:
            print "Error happen with given position smaller than region start position!!!!"
            return
        i_hit = lth / bin_size
        if i_hit >= len(l_hit_cnt):
            print "Error happen with bin index out of range!!!"
            #print pos, start_pos, len(l_hit_cnt) ########################################################################
            return
        l_hit_cnt[i_hit] += 1

    # calculate the weight of each TE region
    def _calc_weight_of_source_region(self, max_cell_cnt, max_cell_weight, l_hit_cnt):
        weight = 0.0
        for bin_cnt in l_hit_cnt:
            if bin_cnt >= max_cell_cnt:
                weight += float(max_cell_weight)
            else:
                weight += float(bin_cnt) * float(max_cell_weight) / float(max_cell_cnt)
        return weight


    # construct the weighted graph for each site
    # not only consider the total number of counts, but also consider the region information
    def _construct_weighted_network_for_site(self, site_chrm, site_pos, m_sites_algnmnts, sf_algnmt):
        m_weighted_network = {}
        bin_size = 100  #######################################hard code here###########################################
        max_cell_cnt = 5
        max_cell_weight = 1.0

        s_algnmt_fa = sf_algnmt + ".fasta"
        sf_fa = open(s_algnmt_fa, "w")
        for barcode in m_sites_algnmnts:
            for alnmt in m_sites_algnmnts[barcode]:
                # f_sam.write(alnmt)  ###write the alignment to file
                # write to fastq file
                xalmt = XAlignment()
                s_fa = xalmt._cvt_algnmt_to_fa(alnmt)
                sf_fa.write(s_fa)

                ##get the map position and check whether mapped to a repeat region
                chrm = alnmt.get_tag('ch')
                pos = int(alnmt.get_tag('cp'))

                mate_chrm = ""
                mate_pos = 0
                if alnmt.mate_is_unmapped == False:
                    mate_chrm = alnmt.get_tag('nh')
                    mate_pos = int(alnmt.next_reference_start)
                ##check the mapping position, whether belong to L1 region
                ##collect the positions

                bhit, region_pos = self.xannotation.is_within_repeat_region(chrm, pos)
                bmate_hit, mate_region_pos = self.xannotation.is_within_repeat_region(mate_chrm, mate_pos)
                if bhit == True:
                    s_hit_id = site_chrm + S_DELIM + str(site_pos) + S_DELIM + chrm + S_DELIM + str(region_pos)
                    if s_hit_id not in m_weighted_network:
                        region_lth = self.xannotation.get_annotation_TE_region_lenth(chrm, region_pos)
                        m_weighted_network[s_hit_id] = self._init_bins_for_region(region_lth, bin_size)
                    else:
                        self._increase_hit_cnt_at_pos(pos, region_pos, bin_size, m_weighted_network[s_hit_id])

                if bmate_hit == True:
                    s_mate_hit_id = site_chrm + S_DELIM + str(site_pos) + S_DELIM + mate_chrm + S_DELIM + str(
                        mate_region_pos)
                    if s_mate_hit_id not in m_weighted_network:
                        region_lth = self.xannotation.get_annotation_TE_region_lenth(mate_chrm, mate_region_pos)
                        m_weighted_network[s_mate_hit_id] = self._init_bins_for_region(region_lth, bin_size)
                    else:
                        self._increase_hit_cnt_at_pos(mate_pos, mate_region_pos, bin_size,
                                                      m_weighted_network[s_mate_hit_id])

                ##for a pair, the relationship is only saved for one time
                if bhit == True and bmate_hit == True and alnmt.is_read1:
                    s_connect_id = chrm + S_DELIM + str(region_pos) + S_DELIM + mate_chrm + S_DELIM + str(
                        mate_region_pos)
                    region_lth1 = self.xannotation.get_annotation_TE_region_lenth(chrm, region_pos)
                    region_lth2 = self.xannotation.get_annotation_TE_region_lenth(mate_chrm, mate_region_pos)

                    short_region_lth = region_lth1
                    short_pos = pos
                    short_region_pos = region_pos
                    if region_lth2 < region_lth1:
                        short_region_lth = region_lth2
                        short_pos = mate_pos
                        short_region_pos = mate_region_pos

                    if s_connect_id not in m_weighted_network:
                        m_weighted_network[s_connect_id] = self._init_bins_for_region(short_region_lth, bin_size)
                    else:
                        self._increase_hit_cnt_at_pos(short_pos, short_region_pos, bin_size,
                                                      m_weighted_network[s_connect_id])
        sf_fa.close()

        G = nx.MultiDiGraph() # construct the graph
        for s_connect_id in m_weighted_network:
            weight = self._calc_weight_of_source_region(max_cell_cnt, max_cell_weight, m_weighted_network[s_connect_id])
            fields = s_connect_id.split(S_DELIM)
            left_node = fields[0] + ":" + fields[1]
            right_node = fields[2] + ":" + fields[3]
            G.add_edge(left_node, right_node, weight=weight)
        return G


    def trace_source_for_one_site(self, chrm, pos, i_extend):
        xbam=XBamInfo(self.sf_bam, self.sf_bam_barcode, self.sf_reference)
        m_site_alignments=xbam.parse_alignments_for_one_site(chrm, pos, i_extend)
        # 3. construct the network for the alignments
        # 3.1. at the same time write out alignment to file
        sf_algnmt = self.s_working_folder + "{}_{}_algnmts.bam".format(chrm, pos)
        g_network = self._construct_simple_network_for_site(chrm, pos, m_site_alignments, sf_algnmt)
        # g_network = self._construct_weighted_network_for_site(chrm, pos, m_site_alignments, sf_algnmt)
        ##here remove the edges with weight<=2
        sf_gexf = "{0}{1}_{2}.gexf".format(self.s_working_folder, chrm, pos)
        nx.write_gexf(g_network, sf_gexf)  ##write to gexf file for visulization in Gephi


    ##for each insertion site, first find out those unique mapped reads, whose mates are mapped to other TE copies
    ##Then, get the barcodes of these reads
    ##Thrid, get the reads of the collected barcodes
    ##Four, check the relationship according to the barcodes
    ##"i_region_lth"
    def parse_alignments_for_insert_sites(self, i_region_lth):
        for chrm in self.m_sites:
            for insert_pos in self.m_sites[chrm]:
                return
####

    # import pandas as pd
    # import numpy as np
    # import vcf
    #
    # def vcf_input(path, c):
    #     '''read a vcf file and extract the information for only one Chromosome you choose.
    #     Puts the information into an hash and makes a DataFrame for all genotypes and samples.'''
    #     elem = ''
    #     dicchrom = {}  # Hash{'key' = [list]}  //dicchrom{['Key = samples'] = 'Value = genotype'
    #     help1 = []
    #
    #     vcf_reader = vcf.Reader(open(path, 'r'))
    #     for record in vcf_reader:
    #         for sample in record.samples:
    #             if sample['GT'] == '.|.' or sample['GT'] == './.':
    #                 elem = None
    #             elif sample['GT'] == '0|0' or sample['GT'] == '0/0':  # homozygot ref
    #                 elem = 0
    #             elif sample['GT'] == '1|0' or sample['GT'] == '1/0' or sample['GT'] == '0|1' or sample[
    #                 'GT'] == '0/1':  # heterozygot
    #                 elem = 1
    #             elif sample['GT'] == '1|1' or sample['GT'] == '1/1':  # homozygot alt
    #                 elem = 2
    #             if sample.sample in dicchrom:
    #                 dicchrom[sample.sample].append(elem)
    #             else:
    #                 help1 = [elem]
    #                 dicchrom[sample.sample] = help1
    #     data = pd.DataFrame(dicchrom)
    #     return data
    #
    # sf_vcf = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/NA12878_10x_v2/hybrid_true_positive.vcf"
    # c = "LINE1"
    # vcf_data = vcf_input(sf_vcf, c)
    # print vcf_data
