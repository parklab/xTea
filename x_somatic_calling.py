##09/05/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

####This is a stand alone module can run by itself;
####It has the "main" function at the bottom.

import os
from optparse import OptionParser
from x_clip_disc_filter import *

class CaseControlMode():
    def __init__(self, s_wfolder):
        self.s_wfolder=s_wfolder
        if os.path.exists(s_wfolder)==False:
            self.s_wfolder="./"
        elif s_wfolder[-1]!="/":
            s_wfolder+="/"

    ####load in the match between sample_id and the id parsed from bam
    def load_in_sample_id_vs_used_one(self, sf_id):
        m_indpdt_vs_bam_id={}
        m_bam_id_vs_indpdt={}
        with open(sf_id) as fin_id:
            for line in fin_id:
                fields=line.split()
                indpdt_id=fields[0]
                id_from_bam=fields[1]
                m_indpdt_vs_bam_id[indpdt_id]=id_from_bam
                m_bam_id_vs_indpdt[id_from_bam]=indpdt_id
        return m_indpdt_vs_bam_id, m_bam_id_vs_indpdt
####

    #for each line: "s_id sf_case sf_control"
    def load_in_samples(self, sf_idx):
        m_sf_case_control={}
        with open(sf_idx) as fin_idx:
            for line in fin_idx:
                fields=line.split()
                s_id=fields[0]
                sf_case=fields[1]
                sf_control=fields[2]
                m_sf_case_control[s_id]=(sf_case, sf_control)
        return m_sf_case_control

    # load the sites from file
    def load_sites(self, sf_ce):
        m_MEI = {}
        with open(sf_ce) as fin_ce:
            for line in fin_ce:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])
                sinfo="\t".join(fields[2:])
                if chrm not in m_MEI:
                    m_MEI[chrm] = {}
                if pos not in m_MEI[chrm]:
                    m_MEI[chrm][pos] = sinfo
        return m_MEI
####

    # #filter out the germline events from 1000G released results
    def cmp_case_control(self, m_case, m_control, islack):
        m_candidates={}
        for chrm in m_case:#doesn't have this chrm
            if chrm not in m_control:
                if chrm not in m_candidates:
                    m_candidates[chrm]={}
                for pos in m_case[chrm]:
                    m_candidates[chrm][pos]=m_case[chrm][pos]
            else:###have this chrm
                for pos in m_case[chrm]:
                    b_exist = False
                    for i in range(-1 * islack, islack):
                        if (pos + i) in m_control[chrm]:
                            b_exist = True
                            break
                    if b_exist==False:
                        if chrm not in m_candidates:
                            m_candidates[chrm]={}
                        m_candidates[chrm][pos]=m_case[chrm][pos]
        return m_candidates

####
    def _load_in_rslt_list(self, sf_rslt_list, b_process):
        m_rslt_list={}
        with open(sf_rslt_list) as fin_rslt_list:
            for line in fin_rslt_list:
                fields=line.split()
                sid=fields[0]
                s_rslt=fields[1]
                if b_process==True:
                    sid_fields=sid.split("_")
                    sid="_".join(sid_fields[:-1])

                m_rslt_list[sid]=s_rslt
        return m_rslt_list

####
    def prepare_case_control_idx(self, sf_id_vs_id, slist_case_rslt, slist_control_rslt, sf_out_idx):
        m_indpdt_vs_bam_id, m_bam_id_vs_indpdt=self.load_in_sample_id_vs_used_one(sf_id_vs_id)
        #print m_bam_id_vs_indpdt
        m_case_list=self._load_in_rslt_list(slist_case_rslt, False)
        #print m_case_list
        m_control_list=self._load_in_rslt_list(slist_control_rslt, True)
        print m_control_list
        with open(sf_out_idx,"w") as fout_idx:
            for s_bam_id in m_bam_id_vs_indpdt:
                if s_bam_id not in m_case_list:
                    print s_bam_id, "do not have a case result!!!"
                    continue
                sf_case=m_case_list[s_bam_id]

                s_sample_id=m_bam_id_vs_indpdt[s_bam_id]
                if s_sample_id not in m_control_list:
                    print "{0} do not have a control!!!".format(s_sample_id)
                    continue
                sf_control=m_control_list[s_sample_id]

                sinfo="{0}\t{1}\t{2}\n".format(s_bam_id, sf_case, sf_control)
                fout_idx.write(sinfo)

####
    def call_somatic_cases(self, sf_id_vs_id, slist_case_rslt, slist_control_rslt, islack):
        sf_idx=self.s_wfolder+"matched_case_control_rslt.list"
        self.prepare_case_control_idx(sf_id_vs_id, slist_case_rslt, slist_control_rslt, sf_idx)

        m_sf_case_control=self.load_in_samples(sf_idx)
        for s_id in m_sf_case_control:
            sf_case=m_sf_case_control[s_id][0]
            sf_control=m_sf_case_control[s_id][1]
            m_case=self.load_sites(sf_case)
            m_control=self.load_sites(sf_control)

            m_candidates=self.cmp_case_control(m_case, m_control, islack)
            sf_out=self.s_wfolder+s_id+".somatic"
####
            with open(sf_out, "w") as fout_rslt:
                for chrm in m_candidates:
                    for pos in m_candidates[chrm]:
                        sinfo="{0}\t{1}\t{2}\n".format(chrm, pos, m_candidates[chrm][pos])
                        fout_rslt.write(sinfo)
####

####
class SomaticMEICaller():#
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
                chrm1 = fields[0]
                chrm=self.process_chrm_name(chrm1)
                pos = int(fields[1])
                if chrm not in m_MEI:
                    m_MEI[chrm] = {}
                if pos not in m_MEI[chrm]:
                    m_MEI[chrm][pos] = line
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
    ##filter out the germline events from germline database
    def filter_by_given_list(self, sf_db, m_MEIs, islack):
        m_candidates={}
        m_db=self.load_sites(sf_db)
        for chrm in m_MEIs:#doesn't have this chrm
            if chrm not in m_db:
                if chrm not in m_candidates:
                    m_candidates[chrm]={}
                for pos in m_MEIs[chrm]:
                    m_candidates[chrm][pos]=1
            else:###have this chrm
                for pos in m_MEIs[chrm]:
                    b_exist = False
                    for i in range(-1 * islack, islack):
                        if (pos + i) in m_db[chrm]:
                            b_exist = True
                            break
                    if b_exist==False:
                        if chrm not in m_candidates:
                            m_candidates[chrm]={}
                        m_candidates[chrm][pos]=m_MEIs[chrm][pos]
        return m_candidates

    #export the selected candidates to a file
    def export_slcted_candidates(self, m_MEIs, sf_out):
        with open(sf_out, "w") as fout_mei:
            for chrm in m_MEIs:
                for pos in m_MEIs[chrm]:
                    s_line=m_MEIs[chrm][pos]
                    fout_mei.write(s_line)

####
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

# ####
# def parse_option():
#     parser = OptionParser()
#     parser.add_option("-i", "--input", dest="input",
#                       help="input file ", metavar="FILE")
#     parser.add_option("--case", dest="case",
#                       help="case results list", metavar="FILE")
#     parser.add_option("--control", dest="control",
#                       help="control results list", metavar="FILE")
#     parser.add_option("-p", "--wfolder", dest="wfolder", type="string",
#                       help="Working folder")
#     parser.add_option("-e", "--slack", dest="slack", type="int",
#                       help="slack value for comparing two positions")
#     (options, args) = parser.parse_args()
#     return (options, args)
#
# ####
# ####
# if __name__ == '__main__':
#     (options, args) = parse_option()
#     sf_vs_ids=options.input
#     sf_case_list=options.case
#     sf_control_list=options.control
#     s_wfolder=options.wfolder
#     i_slack=options.slack
#
#     case_control_mode=CaseControlMode(s_wfolder)
#     case_control_mode.call_somatic_cases(sf_vs_ids, sf_case_list, sf_control_list, i_slack)
