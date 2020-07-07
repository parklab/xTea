import os
import sys
import csv
from x_annotation import *

class XAnalysis():
    def cnt_fall_in_rep_extnd(self, sf_rmsk, sf_list, iextnd, sf_out):
        xanno = XAnnotation(sf_rmsk)
        xanno.load_rmsk_annotation_with_divgnt_extnd(iextnd)
        xanno.index_rmsk_annotation()

        l_in_rep = []
        l_out_rep = []
        with open(sf_list) as fin_list:
            for line in fin_list:
                fields = line.split()
                chrm_tmp = fields[0]
                chrm = chrm_tmp
                if len(chrm_tmp) < 3 or chrm_tmp[:3] != "chr":
                    chrm = "chr" + chrm_tmp
                pos = int(fields[1])
                b_in_rep, rep_pos = xanno.is_within_repeat_region(chrm, pos)
                s_rep = "{0}_{1}".format(chrm, pos)
                if b_in_rep == True: ###fall in the repeat
                    div_rate, sub_family, family, ref_start, ref_end = xanno.get_div_subfamily(chrm, rep_pos)
                    idist=abs(ref_start-pos)
                    if abs(ref_end-pos)<idist:
                        idist=abs(ref_end-pos)
                    l_rep = [s_rep, div_rate, idist, sub_family, family]
                    l_in_rep.append(l_rep)
                else:
                    l_out_rep.append(s_rep)
        self.dump_to_csv(l_in_rep, sf_out)
        return l_out_rep

    #dump the results to csv file
    def dump_to_csv(self, res, csvfile):
        # Assuming res is a list of lists
        with open(csvfile, "w") as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerows(res)

    #
    def cnt_fall_in_rep(self, sf_rmsk, sf_list, sf_out):
        l_extnd=[]
        l_extnd.append(0)
        l_extnd.append(50)
        l_extnd.append(150)
        l_extnd.append(300)

        for iextnd in l_extnd:
            sf_out_tmp="{0}.{1}.csv".format(sf_out, iextnd)
            l_out_rep=self.cnt_fall_in_rep_extnd(sf_rmsk, sf_list, iextnd, sf_out_tmp)
            print iextnd, len(l_out_rep)

