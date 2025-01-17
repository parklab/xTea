##06/15/2022
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

####
#get the cleavage site sequence with given postions and the reference genome
####
import global_values
import pysam
import sys

class CleavageSite():
    def __init__(self, sf_ref):
        self.sf_ref=sf_ref

    #
    def collect_cleavage_site_seqs_for_sites(self, sf_sites, sf_out):
        with open(sf_sites) as fin_sites:
            l_sites_sense=[]
            l_sites_antisense=[]
            for line in fin_sites:
                #chr1    1375678 .       two_side_tprt_both
                if global_values.TWO_SIDE_TPRT_BOTH not in line:
                    continue
                fields=line.rstrip().split()
                s_chrm=fields[0]
                s_pos=fields[1]
                s_ori=fields[2]
                t_rcd=(s_chrm, s_pos, s_ori)
                if s_ori=="+":
                    l_sites_sense.append(t_rcd)
                elif s_ori=="-":
                    l_sites_antisense.append(t_rcd)

        m_seqs_sense, m_seq_sense_freq=self.parse_cleavage_site_seq_given_pos(l_sites_sense)
        m_seqs_antisense, m_seq_antisense_freq=self.parse_cleavage_site_seq_given_pos(l_sites_antisense)

        #merge the positive and negative ones
        m_seq_freq={}
        for s_seq in m_seq_antisense_freq:
            i_rc_freq=m_seq_antisense_freq[s_seq]
            if s_seq not in m_seq_freq:
                m_seq_freq[s_seq]=i_rc_freq
            else:
                m_seq_freq[s_seq] += i_rc_freq
        ##
        for s_seq in m_seq_sense_freq:
            i_freq=m_seq_sense_freq[s_seq]
            s_seq_rc = self._gnrt_reverse_complementary(s_seq)
            if s_seq_rc not in m_seq_freq:####
                m_seq_freq[s_seq_rc]=i_freq
            else:
                m_seq_freq[s_seq_rc] += i_freq

        sf_out_positive=sf_out+".positive"
        with open(sf_out_positive,"w") as fout_positive:
            for s_seq in m_seq_sense_freq:
                fout_positive.write(s_seq+"\t"+str(m_seq_sense_freq[s_seq])+"\n")
        sf_out_negative=sf_out+".negative"
        with open(sf_out_negative,"w") as fout_negative:
            for s_seq in m_seq_antisense_freq:
                fout_negative.write(s_seq+"\t"+str(m_seq_antisense_freq[s_seq])+"\n")
        with open(sf_out,"w") as fout:
            for s_seq in m_seq_freq:
                fout.write(s_seq+"\t"+str(m_seq_freq[s_seq])+"\n")

        #output sites cleavage_seqs
        sf_out_sites_positive=sf_out+".positive_sites"
        with open(sf_out_sites_positive,"w") as fout_sites_positive:
            for s_site in m_seqs_sense:
                s_seq=m_seqs_sense[s_site][0]
                fout_sites_positive.write(s_site+"\t"+s_seq+"\n")
        sf_out_sites_negative = sf_out + ".negative_sites"
        with open(sf_out_sites_negative, "w") as fout_sites_negative:
            for s_site in m_seqs_antisense:
                s_seq = m_seqs_antisense[s_site][0]
                fout_sites_negative.write(s_site + "\t" + s_seq + "\n")
    ####
    #generate the position weighted matrix as input of
    ##    V1  V2  V3  V4  V5  V6  V7  V8
    ## 1 0.0 0.0 0.0 0.3 0.2 0.0 0.0 0.0
    ## 2 0.8 0.2 0.8 0.3 0.4 0.2 0.8 0.2
    ## 3 0.2 0.8 0.2 0.4 0.3 0.8 0.2 0.8
    ## 4 0.0 0.0 0.0 0.0 0.1 0.0 0.0 0.0
    def gnrt_pwm(self, sf_seq_freq, i_seq_len, sf_out):
        with open(sf_seq_freq) as fin_seq_freq, open(sf_out,"w") as fout_matrix:
            l_pos_freq=[]
            l_pos_total_cnt=[]##
            for i in range(i_seq_len):
                m_one_pos_freq={'A':0, 'C':0, 'G':0, 'T':0}
                l_pos_freq.append(m_one_pos_freq)
                l_pos_total_cnt.append(0)#
            for line in fin_seq_freq:
                fields=line.rstrip().split()
                s_seq=fields[0]
                i_freq=int(fields[1])
                if len(s_seq)!=i_seq_len:
                    continue
                for i in range(i_seq_len):
                    l_pos_freq[i][s_seq[i]]+=i_freq
                    l_pos_total_cnt[i]+=i_freq
            ####
            for s_char in ['A','C','G','T']:
                s_info=""
                for i in range(i_seq_len):
                    s_char_freq=str(float(l_pos_freq[i][s_char])/float(l_pos_total_cnt[i]))
                    if s_info=="":
                        s_info=s_char_freq
                    else:
                        s_info += (" "+s_char_freq)
                fout_matrix.write(s_info+"\n")
####

    ##
    def parse_cleavage_site_seq_given_pos(self, l_sites):
        m_cleavage_seq={}
        m_seq_freq={}
        f_fa = pysam.FastaFile(self.sf_ref) #
        m_ref_chrms = {}
        for tmp_chrm in f_fa.references:
            m_ref_chrms[tmp_chrm] = 1
        b_with_chr = False
        if "chr1" in m_ref_chrms:
            b_with_chr = True
        for rcd in l_sites:#
            s_cleavage_seq = "NULL"
            ins_chrm = rcd[0]
            ipos = int(rcd[1])
            s_ori = rcd[2]
            istart=ipos-2
            iend=ipos+5

            ref_chrm = self.process_chrm_name(ins_chrm, b_with_chr)
            if s_ori=="-":#if reverse-complementary
                istart=ipos-5
                iend=ipos+2

            s_cleavage_seq = f_fa.fetch(ref_chrm, istart, iend)
            s_id=ins_chrm+"~"+str(ipos)
            m_cleavage_seq[s_id]=(s_cleavage_seq, s_ori)
            if s_cleavage_seq not in m_seq_freq:
                m_seq_freq[s_cleavage_seq]=1
            else:
                m_seq_freq[s_cleavage_seq] += 1
        f_fa.close()
        return m_cleavage_seq, m_seq_freq
####

    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "b_with_chr"
    def process_chrm_name(self, chrm, b_with_chr):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True
        # print chrm, self.b_with_chr, b_chrm_with_chr ###############################################################
        if b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    ####
    def _set_chr_map(self, chr_map):
        chr_map['A'] = 'T'
        chr_map['a'] = 'T'
        chr_map['T'] = 'A'
        chr_map['t'] = 'A'
        chr_map['C'] = 'G'
        chr_map['c'] = 'G'
        chr_map['G'] = 'C'
        chr_map['g'] = 'C'
        chr_map['N'] = 'N'
        chr_map['n'] = 'N'

    def _gnrt_reverse_complementary(self, s_seq):
        chr_map = {}
        self._set_chr_map(chr_map)
        s_rc = ""
        for s in s_seq[::-1]:
            if s not in chr_map:
                s_rc += "N"
            else:
                s_rc += chr_map[s]
        return s_rc

####
if __name__ == '__main__':
    sf_ref=sys.argv[1]
    sf_sites=sys.argv[2]
    sf_out=sys.argv[3]
    cs=CleavageSite(sf_ref)
    cs.collect_cleavage_site_seqs_for_sites(sf_sites, sf_out)
    i_seq_len=7
    sf_out_matrix=sf_out+".freq_matrix"
    cs.gnrt_pwm(sf_out, i_seq_len, sf_out_matrix)
####

####