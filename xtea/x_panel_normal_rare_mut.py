##08/07/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

#this is to create an extra filter for somatic calling using a panel of normal
#Or use panel normal to define rare mutations

#Extra filters:
#1. filter out "not_hit_end_of_consensus"
#2. filter out "both hit end of consensus"
from optparse import OptionParser
#
class XPanelNormal():
    def __init__(self, sf_panel_normal):
        self.sf_panel_normal=sf_panel_normal
        self.m_panel_normal_sites={}
        self.l_pop_code=['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'WAS', 'OCN', 'CAS']
        #self.m_pop_map={}

    #def _set_pop_map(self):
    #    self.m_pop_map['SAN']='SAS'
    #    self.m_pop_map['ASN']='EAS'

    def load_panel_normal_vcf(self):
        with open(self.sf_panel_normal) as fin_vcf:
            for line in fin_vcf:
                if len(line)>0 and line[0]=="#":
                    continue
                fields=line.rstrip().split("\t")
                if len(fields)<8:
                    print("Error line fields is smaller than 8: ", line.rstrip())
                    continue
                s_chrm=fields[0]
                i_pos=int(fields[1])
                if s_chrm not in self.m_panel_normal_sites:
                    self.m_panel_normal_sites[s_chrm]={}
                self.m_panel_normal_sites[s_chrm][i_pos]={}

                sinfo = fields[7]
                info_fields = sinfo.split(";")
####
                for sterm in info_fields:
                    info_sub_fields = sterm.split("=")
                    s_term_id = info_sub_fields[0]
                    s_term_value = info_sub_fields[1]

                    for s_pop_code in self.l_pop_code:
                        if s_pop_code + "_AF" == s_term_id:
                            f_term_value=float(s_term_value)
                            self.m_panel_normal_sites[s_chrm][i_pos][s_pop_code]=f_term_value

####
    def call_rare_mutations(self, sf_in_vcf, i_slack, f_af_cutoff, sf_rare_vcf):
        with open(sf_in_vcf) as fin_in_vcf, open(sf_rare_vcf,"w") as fout_rare:
            for line in fin_in_vcf:
                if len(line)>0 and line[0]=="#":
                    fout_rare.write(line.rstrip()+"\n")
                    continue
                fields=line.rstrip().split("\t")
                s_chrm=fields[0]
                i_pos=int(fields[1])
                ####
                if s_chrm not in self.m_panel_normal_sites:
                    continue
                i_tmp_start=i_pos-i_slack
                i_tmp_end=i_pos+i_slack
                b_hit=True

                for i_tmp in range(i_tmp_start, i_tmp_end):
                    if i_tmp in self.m_panel_normal_sites[s_chrm]:
                        for s_tmp_pop_code in self.m_panel_normal_sites[s_chrm][i_tmp]:
                            f_tmp_af=self.m_panel_normal_sites[s_chrm][i_tmp][s_tmp_pop_code]
                            if f_tmp_af>f_af_cutoff:
                                b_hit=False
                                break
                        if b_hit==False:
                            break
                if b_hit==True:
                    fout_rare.write(line.rstrip()+"\n")
####
####
def parse_option():
    parser = OptionParser()
    parser.add_option("--panel", dest="panel", type="string",
                      help="Reference panel VCF")
    parser.add_option("-i", "--input", dest="input", type="string",
                      help="The merged VCF")
    parser.add_option("-c", "--cutoff", dest="cutoff", type="float", default=0.01,
                      help="Lower bound of allel frequency")
    parser.add_option("--slack", dest="slack", type="int", default=150,
                      help="Slack value")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    (options, args) = parser.parse_args()
    return (options, args)

####
if __name__ == '__main__':
    (options, args) = parse_option()
    sf_ref_panel=options.panel
    sf_vcf=options.input
    sf_out=options.output
    f_af_cutoff=options.cutoff
    i_slack=options.slack

    panel_normal=XPanelNormal(sf_ref_panel)
    panel_normal.load_panel_normal_vcf()
    panel_normal.call_rare_mutations(sf_vcf, i_slack, f_af_cutoff, sf_out)
####