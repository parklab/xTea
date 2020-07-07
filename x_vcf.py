##11/04/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

#input is a list of xTEA raw output
#output is a merged vcf


class XVCF():
    # for output: each insertion is one record
    def cvt_raw_to_vcf(self, sf_rslt_list, i_peak_win, sf_vcf, min_ins_len=0):
        # first load in all the results from file
        m_samples, m_all_sites = self._load_in_raw_rslts(sf_rslt_list, min_ins_len)
        # each record in format: [chrm][pos]: [list of records]
        m_merged_sites = self.call_peak_candidate_sites_vcf(m_all_sites, i_peak_win)
        self.export_to_vcf(m_samples, m_merged_sites, sf_vcf)

    ####This is used to draw the pca by population
    def cvt_raw_rslt_for_pca_by_sample(self, sf_rslt_list, i_peak_win, sf_sample_info, sf_csv, min_ins_len=0):
        # first load in all the results from file
        m_samples, m_all_sites = self._load_in_raw_rslts(sf_rslt_list, min_ins_len)
        # each record in format: [chrm][pos]: [list of records]
        m_merged_sites = self.call_peak_candidate_sites_vcf(m_all_sites, i_peak_win)
        self.export_to_pca_input(sf_sample_info, m_merged_sites, sf_csv)

    ####
    def export_to_vcf(self, m_samples, m_merged_sites, sf_vcf):
        l_samples = list(m_samples.keys())
        with open(sf_vcf, "w") as fout_vcf:
            l_chrms = list(m_merged_sites.keys())
            l_chrms.sort()
            for chrm in l_chrms:  #
                l_pos = list(m_merged_sites[chrm].keys())
                l_pos.sort()  ###sort the candidate sites
                for pos in l_pos:
                    sinfo = "%s\t%d" % (chrm, pos)
                    m_hit_sample = {}
                    for rcd in m_merged_sites[chrm][pos]:
                        s_sample = rcd[0]
                        i_ins_len = int(rcd[1])
                        m_hit_sample[s_sample] = i_ins_len
                    for s_tmp_sample in l_samples:
                        if s_tmp_sample in m_hit_sample:
                            sinfo += "\t1"
                        else:
                            sinfo += "\t0"
                    sinfo += "\n"
                    fout_vcf.write(sinfo)
####
    ####
    def export_to_pca_input(self, sf_sample_info, m_merged_sites, sf_output):
        m_sample_info = self._load_in_sample_info(sf_sample_info)
        # print m_sample_info
        m_sample_sites = {}
        l_sites = []

        l_chrms = list(m_merged_sites.keys())
        l_chrms.sort()
        for chrm in l_chrms:  #
            l_pos = list(m_merged_sites[chrm].keys())
            l_pos.sort()  ###sort the candidate sites
            for pos in l_pos:
                l_rcd = m_merged_sites[chrm][pos]
                for (s_sample, ilen) in l_rcd:
                    if s_sample not in m_sample_sites:
                        m_sample_sites[s_sample] = {}
                    s_site = "%s_%d" % (chrm, pos)
                    l_sites.append(s_site)
                    m_sample_sites[s_sample][s_site] = 1

        with open(sf_output, "w") as fout_rslt:
            s_head = "sample"
            for s_site in l_sites:
                s_head += ("," + s_site)
            s_head + "\n"
            fout_rslt.write(s_head)
            for s_sample in m_sample_sites:
                if s_sample not in m_sample_info:
                    print s_sample, "not found in dict"
                    continue
                sample_info_rcd = m_sample_info[s_sample]
                s_region = sample_info_rcd[2]
                sinfo = s_sample + "," + s_region
                for s_tmp_site in l_sites:
                    if s_tmp_site in m_sample_sites[s_sample]:
                        sinfo += (",1")
                    else:
                        sinfo += (",0")
                sinfo += "\n"
                fout_rslt.write(sinfo)

                ####

    # If breakpoints are close each other, then should be merged or grouped
    # For each cluster, use the peak one as the representative one
    def call_peak_candidate_sites_vcf(self, m_candidate_sites, peak_window):
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
                    max_cnt = 0
                    tmp_candidate_pos = 0
                    l_tmp_rcds = []
                    for tmp_pos in set_cluster:
                        tmp_cnt = len(m_candidate_sites[chrm][tmp_pos])  # num of supported samples
                        for tmp_rcd in m_candidate_sites[chrm][tmp_pos]:
                            l_tmp_rcds.append(tmp_rcd)
                        if max_cnt < tmp_cnt:
                            tmp_candidate_pos = tmp_pos  #
                            max_cnt = tmp_cnt
                    set_cluster.clear()
                    if chrm not in m_peak_candidate_sites:
                        m_peak_candidate_sites[chrm] = {}
                    if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                        m_peak_candidate_sites[chrm][tmp_candidate_pos] = []
                    for tmp_rcd in l_tmp_rcds:
                        m_peak_candidate_sites[chrm][tmp_candidate_pos].append(tmp_rcd)
                pre_pos = pos
                set_cluster.add(pre_pos)

            # push out the last group
            max_cnt = 0
            tmp_candidate_pos = 0
            l_tmp_rcds = []
            for tmp_pos in set_cluster:
                tmp_cnt = len(m_candidate_sites[chrm][tmp_pos])  # num of supported samples
                for tmp_rcd in m_candidate_sites[chrm][tmp_pos]:
                    l_tmp_rcds.append(tmp_rcd)
                if max_cnt < tmp_cnt:
                    tmp_candidate_pos = tmp_pos  #
                    max_cnt = tmp_cnt
            set_cluster.clear()
            if chrm not in m_peak_candidate_sites:
                m_peak_candidate_sites[chrm] = {}
            if tmp_candidate_pos not in m_peak_candidate_sites[chrm]:
                m_peak_candidate_sites[chrm][tmp_candidate_pos] = []
            for tmp_rcd in l_tmp_rcds:
                m_peak_candidate_sites[chrm][tmp_candidate_pos].append(tmp_rcd)
        return m_peak_candidate_sites

    # this is for this SGD data
    def _load_in_sample_info(self, sf_info):
        m_info = {}
        with open(sf_info) as fin_info:
            for line in fin_info:
                fields = line.split(",")
                s_sample = fields[2]
                s_sex = fields[7]
                s_population = fields[8]
                s_region = fields[9]
                s_latitude = fields[12]
                s_longitude = fields[13]
                m_info[s_sample] = (s_sex, s_population, s_region, s_latitude, s_longitude)
        return m_info
####
    # given raw result list in format: sample-id, results-file-location (each line)
    # this is for exact breakpoints, no merg for nearby sites
    def _load_in_raw_rslts(self, sf_raw_list, min_ins_len=0):
        m_samples = {}
        m_merged_sites = {}
        with open(sf_raw_list) as fin_raw:
            for line in fin_raw:
                fields = line.split()
                s_sample = fields[0]
                sf_rslt = fields[1]
                m_samples[s_sample] = 1

                m_tmp = self._load_one_xtea_rslt(s_sample, sf_rslt, min_ins_len)
                self._combine_dict(m_tmp, m_merged_sites)
        return m_samples, m_merged_sites

    def _combine_dict(self, m_tmp, m_merged):
        for chrm in m_tmp:
            for pos in m_tmp[chrm]:
                for rcd in m_tmp[chrm][pos]:
                    if chrm not in m_merged:
                        m_merged[chrm] = {}
                    if pos not in m_merged[chrm]:
                        m_merged[chrm][pos] = []
                    m_merged[chrm][pos].append(rcd)

    # here only load the sample-id and ins-length for now
    def _load_one_xtea_rslt(self, s_sample, sf_rslt, min_ins_len=0):
        m_rslt = {}
        with open(sf_rslt) as fin_rslt:
            for line in fin_rslt:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])
                ins_lenth = int(fields[-2])
                if ins_lenth < min_ins_len:
                    continue

                if chrm not in m_rslt:
                    m_rslt[chrm] = {}
                if pos not in m_rslt[chrm]:
                    m_rslt[chrm][pos] = []

                m_rslt[chrm][pos].append((s_sample, ins_lenth))
        return m_rslt

####main function
####
##main function
# if __name__ == '__main__':

####