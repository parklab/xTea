##12/03/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu, chong.simon.chu@gmail.com
import os
from l_rep_classification import *
from l_vcf import *

class IHaplotype():
    def __init__(self):
        pass

    #merge the results of two haplotypes
    def merge_xTEA_output_of_two_hap(self, sf_prefix1, sf_prefix2, i_type, swfolder, n_jobs, i_slack, sf_out_prefix):
        lrc = LRepClassification(swfolder, n_jobs)
        l_rslts1 = lrc.get_output_list_by_type(sf_prefix1, i_type)
        l_rslts1.append(sf_prefix1)
        l_rslts2 = lrc.get_output_list_by_type(sf_prefix2, i_type)
        l_rslts2.append(sf_prefix2)
        l_types = lrc.set_rep_type(i_type)
        l_types.append("all")

        n_rep=len(l_rslts1)
        for i in range(n_rep):
            sf_rslt1=l_rslts1[i]
            sf_rslt2=l_rslts2[i]
            # if one or both of the file is absent, then skip
            if os.path.isfile(sf_rslt1)==False or os.path.isfile(sf_rslt2)==False:
                print("Files %s or %s doesn't exist!" % (sf_rslt1, sf_rslt2))
                continue
            s_type=l_types[i]
            sf_merged=sf_out_prefix+"_merged_all_hap_"+s_type+".txt"
            self.merge_one_group(sf_rslt1, sf_rslt2, i_slack, sf_merged)
    #

    def merge_one_group(self, sf_rslt1, sf_rslt2, i_slack, sf_merged):
        raw_rslt=L_Raw_Rslt()
        m_rslt1=raw_rslt.load_in_results2(sf_rslt1)
        m_rslt2=raw_rslt.load_in_results2(sf_rslt2)
        with open(sf_merged,"w") as fout_merged:
            for chrm in m_rslt1:
                if chrm not in m_rslt2:
                    for pos in m_rslt1[chrm]:
                        fout_merged.write(m_rslt1[chrm][pos]+"\t"+"H1\n")
                    continue
                for pos in m_rslt1[chrm]:
                    b_hit=False
                    for tmp_pos in range(pos-i_slack, pos+i_slack):
                        if tmp_pos in m_rslt2[chrm]:
                            fout_merged.write(m_rslt1[chrm][pos] + "\t" + "H1,H2\n")
                            b_hit=True
                            break
                    if b_hit==False:
                        fout_merged.write(m_rslt1[chrm][pos] + "\t" + "H1\n")

            #find those in 2, but absent in 1
            for chrm in m_rslt2:
                if chrm not in m_rslt1:
                    for pos in m_rslt2[chrm]:
                        fout_merged.write(m_rslt2[chrm][pos]+"\t"+"H2\n")
                    continue
                for pos in m_rslt2[chrm]:
                    b_hit=False
                    for tmp_pos in range(pos-i_slack, pos+i_slack):
                        if tmp_pos in m_rslt1[chrm]:
                            b_hit=True
                            break
                    if b_hit==False:
                        fout_merged.write(m_rslt2[chrm][pos] + "\t" + "H2\n")

####