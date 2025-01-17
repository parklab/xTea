##08/27/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

from optparse import OptionParser

l_pop_code=['AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'WAS', 'OCN', 'CAS']
f_base=0.00001
n_level=int(1.0/f_base)+1
#n_high_af=0.1
#n_low_af=0.005

def add_a_site(f_af, s_pop, m_af):
    for i in range(n_level):
        if (float(i)*f_base) > f_af:
            break
        m_af[s_pop][i]+=1

####
def cnt_population_AF(sf_vcf, n_high_af, n_low_af, sf_out_csv):
    m_pop_af={}
    for pop_code in l_pop_code:
        m_pop_af[pop_code]=[]
        for i in range(n_level):
            m_pop_af[pop_code].append(0)

    sf_sites = sf_out_csv + ".high_af.sites"
    with open(sf_vcf) as fin_vcf, open(sf_sites, "w") as fout_sites:
        for line in fin_vcf:
            if line[0]=="#":
                continue
            fields=line.split('\t')
            sinfo=fields[7]
            info_fields=sinfo.split(";")
####
            n_cnt = 0
            n_cnt_larger_upper_bound=0
            n_cnt_larger_lower_bound=0
            s_tmp_pop_code = ""
            s_tmp_af = ""
            for sterm in info_fields:#
                info_sub_fields=sterm.split("=")
                s_term_id=info_sub_fields[0]
                s_term_value=info_sub_fields[1]

                for s_pop_code in l_pop_code:
                    if s_pop_code+"_AF" == s_term_id:
                        n_cnt+=1
                        s_tmp_pop_code=s_pop_code
                        s_tmp_af=s_term_value
                        if float(s_term_value)>n_low_af:
                            n_cnt_larger_lower_bound+=1
                        if float(s_term_value)> n_high_af:
                            n_cnt_larger_upper_bound+=1

            if n_cnt==1:##
                f_tmp_af=float(s_tmp_af)
                add_a_site(f_tmp_af, s_tmp_pop_code, m_pop_af)

            if n_cnt_larger_upper_bound ==1 and n_cnt_larger_lower_bound==1:
                #if f_tmp_af > n_high_af:
                fout_sites.write(line.rstrip()+"\n")

    with open(sf_out_csv, "w") as fout_csv:
        fout_csv.write("AF,Population,Count\n")
        for pop_code in m_pop_af:
            for i in range(n_level):#
                f_tmp=float(i)*f_base
                n_cnt=m_pop_af[pop_code][i]
                sinfo="{0},{1},{2}\n".format(f_tmp, pop_code, n_cnt)
                fout_csv.write(sinfo)
####

def parse_option():
    parser = OptionParser()
    parser.add_option("--high", dest="high_af", type="float", default=0.1,
                      help="Upper bound of allel frequency")
    parser.add_option("--low", dest="low_af", type="float", default=0.005,
                      help="Lower bound of allel frequency")
    parser.add_option("-i", "--input", dest="input", type="string",
                      help="The merged VCF")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == '__main__':
    (options, args) = parse_option()
    sf_vcf=options.input
    sf_out=options.output
    f_upper=options.high_af
    f_lower=options.low_af
    cnt_population_AF(sf_vcf, f_upper, f_lower, sf_out)
####

####