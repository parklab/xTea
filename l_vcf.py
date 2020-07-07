##11/27/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu, chong.simon.chu@gmail.com

import global_values


class L_Raw_Rslt():
    def __init__(self):
        pass

    def load_in_results(self, sf_in):
        m_seq={}
        with open(sf_in) as fin_rslt:
            #each in format:
            # chr13   98867635        LINE1   None    0:6067:+        None    TAGTCATTAGCT    not_transduction
            # not_transduction  seq  None    None    with_polyA
            for line in fin_rslt:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                s_seq=fields[9]
                s_id="%s%s%s" % (chrm, global_values.SEPERATOR, pos)
                m_seq[s_id]=s_seq
        return m_seq
####
    def load_in_results2(self, sf_in):
        m_rslt={}
        with open(sf_in) as fin_rslt:
            #each in format:
            # chr13   98867635        LINE1   None    0:6067:+        None    TAGTCATTAGCT    not_transduction
            # not_transduction  seq  None    None    with_polyA
            for line in fin_rslt:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                if chrm not in m_rslt:
                    m_rslt[chrm]={}
                m_rslt[chrm][pos]=line.rstrip()
        return m_rslt
####