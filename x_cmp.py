import sys

EXTND = 200


class RsltCMP():
    def __init__(self, sf_input, sf_bcmk):
        self.sf_input = sf_input
        self.sf_bcmk = sf_bcmk

    #load the sites
    def load_sites(self, sf_input):
        m_sites = {}
        with open(sf_input) as fin_input:
            for line in fin_input:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_sites:
                    m_sites[chrm] = {}
                if pos not in m_sites[chrm]:
                    m_sites[chrm][pos] = []
                for svalue in fields[2:]:
                    m_sites[chrm][pos].append(svalue)
        return m_sites

    #process chrm name
    def process_chrm_name(self, chrm, b_with_chr):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        # print chrm, self.b_with_chr, b_chrm_with_chr ############################################################################
        if b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    #m_hit save the position of in m_dict1
    def cmp_dict(self, m_dict1, m_bcmk):
        n_overlap = 0
        m_unhit = {}
        m_hit={}

        b_contain_chr = False
        if "chr1" in m_bcmk:
            b_contain_chr = True

        for chrm_1 in m_bcmk:
            chrm = self.process_chrm_name(chrm_1, b_contain_chr)
            if chrm not in m_dict1:
                continue

            for pos in m_bcmk[chrm]:
                b_hit = False
                for i in range(-1 * EXTND, EXTND):
                    if (pos + i) in m_dict1[chrm]:
                        n_overlap += 1
                        b_hit = True
                        if chrm not in m_hit:
                            m_hit[chrm]={}
                        m_hit[chrm][pos+i]=1
                        break

                if b_hit == False:
                    if chrm not in m_unhit:
                        m_unhit[chrm] = {}
                    m_unhit[chrm][pos] = 1
        return n_overlap, m_unhit, m_hit

####
    #cmp results
    def cmp_rslts(self, sf_unhit):#
        m_bcmk = self.load_sites(self.sf_bcmk)
        m_input = self.load_sites(self.sf_input)

        n_bcmk, m_unhit1, m_hit1 = self.cmp_dict(m_bcmk, m_bcmk)#
        n_hit, m_unhit2, m_hit2 = self.cmp_dict(m_input, m_bcmk)
        n_hit_fp, m_unhit_fp, m_hit3 = self.cmp_dict(m_bcmk, m_input)

        print "{0}/{1} hits/all".format(n_hit, n_bcmk)
        with open(sf_unhit, "w") as fout_unhit:
            for chrm in m_unhit2:
                for pos in m_unhit2[chrm]:
                    fout_unhit.write(chrm + "\t" + str(pos) + "\n")

        with open(sf_unhit+".true_positive", "w") as fout_tp:
            for chrm in m_hit2:
                for pos in m_hit2[chrm]:
                    sinfo = "\t".join(m_input[chrm][pos])
                    fout_tp.write(chrm + "\t" + str(pos) + "\t" + sinfo + "\n")

        with open(sf_unhit+".false_positive", "w") as fout_fp:
            for chrm in m_unhit_fp:
                for pos in m_unhit_fp[chrm]:
                    sinfo="\t".join(m_input[chrm][pos])
                    fout_fp.write(chrm + "\t" + str(pos) + "\t" + sinfo+ "\n")

####
if __name__ == '__main__':
    sf_input = sys.argv[1]
    sf_bcmk = sys.argv[2]
    sf_unhit = sys.argv[3]
    rslt_cmp = RsltCMP(sf_input, sf_bcmk)
    rslt_cmp.cmp_rslts(sf_unhit)

