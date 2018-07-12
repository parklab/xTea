##11/27/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
from subprocess import *

NEARBY_REGION = 50
#PEAK_WINDOW = 100

class XFilter():
    def parse_sites_with_clip_cutoff(self, m_clip_pos_freq, cutoff_left_clip, cutoff_right_clip):
        i_clip_mate_in_rep = 2  ########################################################################################
        m_candidate_sites = {}
        for chrm in m_clip_pos_freq:
            for pos in m_clip_pos_freq[chrm]:
                ####here need to check the nearby region
                nearby_left_freq = 0
                nearby_right_freq = 0
                nearby_mate_in_rep = 0
                for i in range(-1 * NEARBY_REGION, NEARBY_REGION):
                    i_tmp_pos = pos + i
                    if i_tmp_pos in m_clip_pos_freq[chrm]:
                        nearby_left_freq += m_clip_pos_freq[chrm][i_tmp_pos][0]
                        nearby_right_freq += m_clip_pos_freq[chrm][i_tmp_pos][1]
                        nearby_mate_in_rep += m_clip_pos_freq[chrm][i_tmp_pos][2]

                # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip \
                #         and nearby_mate_in_rep >= i_clip_mate_in_rep:

                # for this version, doesn't check whether the clipped part is mapped to a repeat region
                if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip:
                    if chrm not in m_candidate_sites:
                        m_candidate_sites[chrm] = {}
                    i_left_cnt = m_clip_pos_freq[chrm][pos][0]
                    i_right_cnt = m_clip_pos_freq[chrm][pos][1]
                    i_mate_in_rep_cnt = m_clip_pos_freq[chrm][pos][2]
                    m_candidate_sites[chrm][pos] = (i_left_cnt, i_right_cnt, i_mate_in_rep_cnt)
        return m_candidate_sites


    def parse_sites_with_clip_cutoff_for_chrm(self, m_clip_pos_freq, cutoff_left_clip, cutoff_right_clip,
                                              cutoff_clip_mate_in_rep):
        m_candidate_sites = {}
        for pos in m_clip_pos_freq:
            ####here need to check the nearby region
            nearby_left_freq = 0
            nearby_right_freq = 0
            nearby_mate_in_rep = 0
            for i in range(-1 * NEARBY_REGION, NEARBY_REGION):
                i_tmp_pos = pos + i
                if i_tmp_pos in m_clip_pos_freq:
                    nearby_left_freq += m_clip_pos_freq[i_tmp_pos][0]
                    nearby_right_freq += m_clip_pos_freq[i_tmp_pos][1]
                    nearby_mate_in_rep += (
                    m_clip_pos_freq[i_tmp_pos][2] + m_clip_pos_freq[i_tmp_pos][3] + m_clip_pos_freq[i_tmp_pos][4])

            b_candidate=False
            # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip \
            #         and nearby_mate_in_rep >= cutoff_clip_mate_in_rep:
            #     b_candidate=True
            if (nearby_left_freq >= cutoff_left_clip or nearby_right_freq >= cutoff_right_clip) \
                    and nearby_mate_in_rep >= cutoff_clip_mate_in_rep:
                b_candidate=True

            if b_candidate==True:
                # if nearby_left_freq >= cutoff_left_clip and nearby_right_freq >= cutoff_right_clip:
                i_left_cnt = m_clip_pos_freq[pos][0]
                i_right_cnt = m_clip_pos_freq[pos][1]
                i_mate_in_rep_cnt = m_clip_pos_freq[pos][2]
                m_candidate_sites[pos] = (i_left_cnt, i_right_cnt, i_mate_in_rep_cnt)
        return m_candidate_sites

    ###output the candidate list in a file
    def output_candidate_sites(self, m_candidate_list, sf_out):
        with open(sf_out, "w") as fout_candidate_sites:
            for chrm in m_candidate_list:
                if self.is_decoy_contig_chrms(chrm):  ####decoy and other contigs are not interested!!!!
                    continue
                for pos in m_candidate_list[chrm]:
                    lth = len(m_candidate_list[chrm][pos])
                    fout_candidate_sites.write(chrm + "\t" + str(pos) + "\t")
                    for i in range(lth):
                        s_feature = str(m_candidate_list[chrm][pos][i])
                        fout_candidate_sites.write(s_feature + "\t")
                    fout_candidate_sites.write("\n")

    def is_decoy_contig_chrms(self, chrm):
        fields = chrm.split("_")
        if len(fields) > 1:
            return True
        elif chrm == "hs37d5":
            return True

        if chrm=="MT" or chrm=="chrMT":#doesn't consider the mitchrondrial DNA
            return True

        dot_fields = chrm.split(".")
        if len(dot_fields) > 1:
            return True
        else:
            return False


    def load_in_candidate_list(self, sf_candidate_list):
        m_list = {}
        with open(sf_candidate_list) as fin_candidate_sites:
            for line in fin_candidate_sites:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_list:
                    m_list[chrm] = {}
                if pos not in m_list[chrm]:
                    m_list[chrm][pos] = []
                for ivalue in fields[2:]:
                    m_list[chrm][pos].append(int(ivalue))
        return m_list

    # In the previous step (call_TEI_candidate_sites), some sites close to each other may be introduced together
    # If there are more than 1 site close to each other, than use the peak site as a representative
    def call_peak_candidate_sites(self, m_candidate_sites, peak_window):
        m_peak_candidate_sites = {}
        for chrm in m_candidate_sites:
            l_pos = list(m_candidate_sites[chrm].keys())
            l_pos.sort()  ###sort the candidate sites
            pre_pos = -1
            set_cluster = set()
            for pos in l_pos:
                if pre_pos == -1:
                    pre_pos = pos
                    set_cluster.add(pre_pos)
                    continue

                if pos - pre_pos > peak_window:  # find the peak in the cluster
                    max_clip = 0
                    tmp_candidate_pos = 0
                    for tmp_pos in set_cluster:
                        tmp_left_clip = int(m_candidate_sites[chrm][tmp_pos][0]) #left clip
                        tmp_right_clip = int(m_candidate_sites[chrm][tmp_pos][1]) #right clip
                        tmp_all_clip = tmp_left_clip + tmp_right_clip #all the clip
                        if max_clip < tmp_all_clip:
                            tmp_candidate_pos = tmp_pos
                            max_clip = tmp_all_clip
                    set_cluster.clear()
                    if chrm not in m_peak_candidate_sites:
                        m_peak_candidate_sites[chrm] = {}
                    if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                        m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
                pre_pos = pos
                set_cluster.add(pre_pos)
            # push out the last group
            max_clip = 0
            tmp_candidate_pos = 0
            for tmp_pos in set_cluster:
                tmp_left_clip = int(m_candidate_sites[chrm][tmp_pos][0])
                tmp_right_clip = int(m_candidate_sites[chrm][tmp_pos][1])
                tmp_all_clip = tmp_left_clip + tmp_right_clip
                if max_clip < tmp_all_clip:
                    tmp_candidate_pos = tmp_pos
                    max_clip = tmp_all_clip
            if chrm not in m_peak_candidate_sites:
                m_peak_candidate_sites[chrm] = {}
            if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                ##Here, use list in order to output the list (by call the output_candidate_sites function)
                m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
        return m_peak_candidate_sites

    def merge_clip_disc(self, sf_disc_tmp, sf_clip, sf_out):
        with open(sf_out, "w") as fout_list:
            m_disc={}
            with open(sf_disc_tmp) as fin_disc:
                for line in fin_disc:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    s_left_disc=fields[2]
                    s_right_disc=fields[3]
                    if chrm not in m_disc:
                        m_disc[chrm]={}
                    m_disc[chrm][pos]=(s_left_disc, s_right_disc)
            with open(sf_clip) as fin_clip:
                for line in fin_clip:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    if chrm not in m_disc:
                        print "Error happen at merge clip and disc feature step: {0} not exist".format(chrm)
                        continue
                    if pos not in m_disc[chrm]:
                        continue
                    s_left_disc=m_disc[chrm][pos][0]
                    s_right_disc=m_disc[chrm][pos][1]
                    fields.append(s_left_disc)
                    fields.append(s_right_disc)
                    fout_list.write("\t".join(fields) + "\n")

    ####This is to output all the candidates with all the clip, discord, barcode information in one single file
    def merge_clip_disc_barcode(self, sf_barcode_tmp, sf_disc, sf_out):
        with open(sf_out, "w") as fout_list:
            m_barcode={}
            with open(sf_barcode_tmp) as fin_barcode:
                for line in fin_barcode:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    s_nbarcode=fields[-1]
                    if chrm not in m_barcode:
                        m_barcode[chrm]={}
                    m_barcode[chrm][pos]=s_nbarcode
            with open(sf_disc) as fin_disc:
                for line in fin_disc:
                    fields=line.split()
                    chrm=fields[0]
                    pos=fields[1]
                    if chrm not in m_barcode:
                        print "Error happen at merge clip, disc and barcode step: {0} not exist".format(chrm)
                        continue
                    if pos not in m_barcode[chrm]:
                        continue
                    s_barcode=m_barcode[chrm][pos]
                    fields.append(s_barcode)
                    fout_list.write("\t".join(fields) + "\n")

    # In the previous step (call_TEI_candidate_sites), some sites close to each other may be introduced together
    # If there are more than 1 site close to each other, than use the peak site as a representative
    def call_peak_candidate_sites_all_features(self, m_candidate_sites, peak_window):
        m_peak_candidate_sites = {}
        for chrm in m_candidate_sites:
            l_pos = list(m_candidate_sites[chrm].keys())
            l_pos.sort()  ###sort the candidate sites
            pre_pos = -1
            set_cluster = set()
            for pos in l_pos:
                if pre_pos == -1:
                    pre_pos = pos
                    set_cluster.add(pre_pos)
                    continue

                if pos - pre_pos > peak_window:  # find the peak in the cluster
                    max_clip = 0
                    max_all=0
                    tmp_candidate_pos = 0
                    for tmp_pos in set_cluster:
                        tmp_clip = int(m_candidate_sites[chrm][tmp_pos][0])  # # ofclip related features
                        tmp_all = int(m_candidate_sites[chrm][tmp_pos][1])  # # of all features
                        if (max_clip < tmp_clip) or (max_clip==tmp_clip and max_all < tmp_all):
                            tmp_candidate_pos = tmp_pos
                            max_clip = tmp_clip
                            max_all=tmp_all
                    set_cluster.clear()
                    if chrm not in m_peak_candidate_sites:
                        m_peak_candidate_sites[chrm] = {}
                    if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                        m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
                pre_pos = pos
                set_cluster.add(pre_pos)
            # push out the last group
            max_clip = 0
            max_all = 0
            tmp_candidate_pos = 0
            for tmp_pos in set_cluster:
                tmp_clip = int(m_candidate_sites[chrm][tmp_pos][0])  # # ofclip related features
                tmp_all = int(m_candidate_sites[chrm][tmp_pos][1])  # # of all features
                if (max_clip < tmp_clip) or (max_clip == tmp_clip and max_all < tmp_all):
                    tmp_candidate_pos = tmp_pos
                    max_clip = tmp_clip
                    max_all = tmp_all
            set_cluster.clear()
            if chrm not in m_peak_candidate_sites:
                m_peak_candidate_sites[chrm] = {}
            if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                m_peak_candidate_sites[chrm][tmp_candidate_pos] = [max_clip]
        return m_peak_candidate_sites

    ###this function is used to combine some sites are close to each other
    def combine_closing_sites(self, sf_input, iwindow, sf_out):
        m_original_sites={}
        with open(sf_input) as fin_sites:
            for line in fin_sites:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                cur_sum_clip=int(fields[2]) + int(fields[3]) + int(fields[4])
                cur_sum_all = cur_sum_clip + int(fields[5]) + int(fields[6])

                if chrm not in m_original_sites:
                    m_original_sites[chrm]={}
                m_original_sites[chrm][pos]=(cur_sum_clip, cur_sum_all)

        m_peak_candidate_sites=self.call_peak_candidate_sites_all_features(m_original_sites, iwindow)
        with open(sf_input) as fin_sites, open(sf_out, "w") as fout_sites:
            for line in fin_sites:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                if (chrm in m_peak_candidate_sites) and (pos in m_peak_candidate_sites[chrm]):
                    fout_sites.write(line)

