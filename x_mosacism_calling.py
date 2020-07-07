import sys


class MosacismTEICaller():
    def process_chrm_name(self, schrm):
        if len(schrm) > 3 and schrm[:3] == "chr":
            return schrm[3:]
        else:
            return schrm

    #load the sites from file
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

    #get the candidate sites happen in both lists
    def get_overlap(self, m1, m2, slack):
        m_shared={}
        for chrm in m1:
            if chrm in m2:
                for pos in m1[chrm]:
                    for i in range(-1 * slack, slack):
                        if (pos + i) in m2[chrm]:
                            if chrm not in m_shared:
                                m_shared[chrm]={}
                            m_shared[chrm][pos]=1
                            break
        return m_shared
####
                            #s_info = "{0}\t{1}\n".format(chrm, pos)
                            #fout_common.write(s_info)

    #filter out the germline events from 1000G released results
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

    #find out the mosaic events with case, control and germline datasets
    def call_mosaic_from_case_control(self, sf_case, sf_control, sf_1kg, islack, sf_out):
        m_all_case=self.load_sites(sf_case)
        m_somatic=self.filter_by_given_list(sf_1kg, m_all_case, islack)
        m_mosaic=self.filter_by_given_list(sf_control, m_somatic, islack)

        ####
        with open(sf_out, "w") as fout_rslt:
            for chrm in m_mosaic:
                for pos in m_mosaic[chrm]:
                    s_info="{0}\t{1}\n".format(chrm, pos)
                    fout_rslt.write(s_info)


if __name__ == '__main__':
    sf_case = sys.argv[1]
    sf_control = sys.argv[2]
    sf_1kg=sys.argv[3]
    islack = int(sys.argv[4])
    sf_out = sys.argv[5]
    mcaller = MosacismTEICaller()
    mcaller.call_mosaic_from_case_control(sf_case, sf_control, sf_1kg, islack, sf_out)

####