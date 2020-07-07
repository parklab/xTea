import os
import pysam

class XPhylogeny():
    def construct_phylogeny(self):
        return

    ##construct the consensus from the muscle alignment
    def gnrt_consensus_from_algnmt(sefl, sf_algn):
        sf_afa = pysam.FastaFile(sf_algn)

        l_ids = sf_afa.references

        m_cns = {}
        for sid in l_ids:
            seq = sf_afa.fetch(sid)
            for (idx, val) in enumerate(seq):
                if idx not in m_cns:
                    m_cns[idx] = {}
                if val not in m_cns[idx]:
                    m_cns[idx][val] = 1
                else:
                    m_cns[idx][val] += 1

        s_cns = ""
        seq_len = len(m_cns)
        for idx in range(seq_len):
            i_max = 0
            c_tmp = ""
            for val in m_cns[idx]:
                if m_cns[idx][val] > i_max:
                    c_tmp = val
                    i_max = m_cns[idx][val]
            s_cns += c_tmp

        s_cns2 = ""
        for c in s_cns:
            if c == "-":
                continue
            s_cns2 += c
        sf_afa.close()
        return s_cns2
