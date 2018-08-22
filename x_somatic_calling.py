import sys
from x_clip_disc_filter import *


class SomaticMEICaller():
    def __init__(self, sf_ref):
        self.reference=sf_ref

    def process_chrm_name(self, schrm):
        if len(schrm) > 3 and schrm[:3] == "chr":
            return schrm[3:]
        else:
            return schrm

    # load the sites from file
    def load_sites(self, sf_ce):
        m_MEI = {}
        with open(sf_ce) as fin_ce:
            for line in fin_ce:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_MEI:
                    m_MEI[chrm] = {}
                if pos not in m_MEI[chrm]:
                    m_MEI[chrm][pos] = 1
        return m_MEI

    # get the candidate sites happen in both lists
    def get_overlap(self, m1, m2, slack):
        m_shared = {}
        for chrm in m1:
            if chrm in m2:
                for pos in m1[chrm]:
                    for i in range(-1 * slack, slack):
                        if (pos + i) in m2[chrm]:
                            if chrm not in m_shared:
                                m_shared[chrm] = {}
                            m_shared[chrm][pos] = 1
                            break
        return m_shared

    # find out the mosaic events with case, control and germline datasets
    def call_mosaic_from_case_control(self, sf_case, sf_control, sf_1kg, islack, sf_out):
        m_all_case = self.load_sites(sf_case)
        m_somatic = self.filter_by_given_list(sf_1kg, m_all_case, islack)
        m_mosaic = self.filter_by_given_list(sf_control, m_somatic, islack)

        ####
        with open(sf_out, "w") as fout_rslt:
            for chrm in m_mosaic:
                for pos in m_mosaic[chrm]:
                    s_info = "{0}\t{1}\n".format(chrm, pos)
                    fout_rslt.write(s_info)

    ####
    # s_info = "{0}\t{1}\n".format(chrm, pos)
    # fout_common.write(s_info)

    # #filter out the germline events from 1000G released results
    def filter_by_given_list(self, sf_1KG, m_MEIs, islack):
        m_candidates={}
        m_1kg=self.load_sites(sf_1KG)
        for chrm in m_MEIs:#doesn't have this chrm
            if chrm not in m_1kg:
                if chrm not in m_candidates:
                    m_candidates[chrm]={}
                for pos in m_MEIs[chrm]:
                    m_candidates[chrm][pos]=1
            else:###have this chrm
                for pos in m_MEIs[chrm]:
                    b_exist = False
                    for i in range(-1 * islack, islack):
                        if (pos + i) in m_1kg[chrm]:
                            b_exist = True
                            break
                    if b_exist==False:
                        if chrm not in m_candidates:
                            m_candidates[chrm]={}
                        m_candidates[chrm][pos]=1
        return m_candidates


    # filter out FP by allele frequency
    def filter_by_AF(self, sf_bam_list, sf_candidate_list, extnd, clip_slack, s_working_folder, n_jobs, af_cutoff,
                     sf_output):
        xcd = XClipDisc(sf_bam_list, s_working_folder, n_jobs, self.reference)
        xcd.calc_AF_by_clip_reads_of_given_list(sf_bam_list, sf_candidate_list, extnd, clip_slack, af_cutoff, sf_output)

    # m_black_list in format: {chrm:[]}
    def is_within_blacklist_region(self, m_black_list, chrm, pos):
        if chrm not in m_black_list:
            return False
        lo, hi = 0, len(m_black_list[chrm]) - 1
        while lo <= hi:
            mid = (lo + hi) / 2
            mid_start_pos = m_black_list[chrm][mid][0]
            mid_end_pos = m_black_list[chrm][mid][1]
            if pos >= mid_start_pos and pos <= mid_end_pos:
                return True
            elif pos > mid_start_pos:
                lo = mid + 1
            else:
                hi = mid - 1
        return False

    # filter out the candidates by black list
    def filter_by_black_list(self, sf_candidate_list, sf_black_list, sf_out):
        m_black_list = {}
        with open(sf_black_list) as fin_black_list:
            for line in fin_black_list:
                fields = line.split()
                chrm_1 = fields[0]
                chrm = self.process_chrm_name(chrm_1)  # remove the "chr" is exists
                start = int(fields[1])
                end = int(fields[2])
                if chrm not in m_black_list:
                    m_black_list[chrm] = []
                m_black_list[chrm].append((start, end))  # assume the list has been sorted

        with open(sf_candidate_list) as fin_list, open(sf_out, "w") as fout_list:
            for line in fin_list:
                fields = line.split()
                chrm_1 = fields[0]
                ins_chrm = self.process_chrm_name(chrm_1)
                ins_pos = int(fields[1])
                # now check whether hit the black list region
                if self.is_within_blacklist_region(m_black_list, ins_chrm, ins_pos) == True:
                    continue
                fout_list.write(line)

    # filter out the germline events from 1000G released results
    def filter_by_germline_list(self, sf_1KG, sf_candidate, islack, sf_out):
        m_1kg = self.load_sites(sf_1KG)
        with open(sf_candidate) as fin_mei, open(sf_out, "w") as fout_rslts:
            for line in fin_mei:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])

                if chrm not in m_1kg:
                    fout_rslts.write(line)
                else:  ###have this chrm
                    b_exist = False
                    for i in range(-1 * islack, islack):
                        if (pos + i) in m_1kg[chrm]:
                            b_exist = True
                            break
                    if b_exist == False:
                        fout_rslts.write(line)

#awk -F "," '{cnt=$2+$3+$4+$5; if(cnt==1){print $0}}' intersect_picked_neuron.csv > intersect_picked_neuron_single.txt
    #def call_mosaic_from_
##########################################################################
    #This version is to call
    def call_somatic_from_case_control(self, sf_case, sf_control, islack, sf_out):
        m_all_case = self.load_sites(sf_case)
