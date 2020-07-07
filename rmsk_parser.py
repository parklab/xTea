##11/17/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu, chong.simon.chu@gmail.com


class RMSK_Parser():
    def __init__(self, sf_rmsk):
        self.sf_rmsk=sf_rmsk

    '''
    Sample intput format:
    SW   perc perc perc  query              position in query    matching  repeat            position in repeat
    score   div. del. ins.  sequence           begin end   (left)   repeat    class/family    begin  end    (left)  ID
    
    1000   18.8  5.2 11.1  10_110453470_ctg1      1   445  (828) + L1M4c     LINE/L1            737   1246 (6271)   1 *
    4241    3.9  1.1 10.1  10_110453470_ctg1    420  1241   (32) + L1HS      LINE/L1           5334   6088   (67)   2
    6913    6.4  6.3  2.8  12_117814446_ctg1     59  2011    (0) C L1HS      LINE/L1            (1)   6154   4143  23
    '''
    def parse_rmsk(self, s_slct_type=""):
        m_hits_by_site={}
        m_hits_by_ctg={}
        with open(self.sf_rmsk) as fin_rmsk:
            for line in fin_rmsk:
                fields=line.split()
                if len(fields)<1:
                    continue
                if fields[0]=="score" or fields[0]=="SW":
                    continue

                ####
                s_contig=fields[4]
                ctg_fields=s_contig.split("_") ####Note: potential bug here
                #print s_contig
                ins_chrm=ctg_fields[0]
                ins_pos=int(ctg_fields[1])

                i_contig_start=int(fields[5])
                i_contig_end=int(fields[6])
                s_rc=fields[8]
                s_sub_family=fields[9]
                s_rep_family=fields[10]

                if ("" != s_slct_type) and (s_slct_type not in s_sub_family):
                    continue

                i_cns_start=-1
                i_cns_end=-1
                b_rc = False
                if s_rc=="+":
                    i_cns_start=int(fields[11])
                    i_cns_end=int(fields[12])
                else:
                    i_cns_start=int(fields[13])
                    i_cns_end=int(fields[12])
                    b_rc=True

                if ins_chrm not in m_hits_by_site:
                    m_hits_by_site[ins_chrm]={}
                if ins_pos not in m_hits_by_site[ins_chrm]:
                    m_hits_by_site[ins_chrm][ins_pos]=[]
                hit_rcd=(i_contig_start, i_contig_end, b_rc, s_sub_family, s_rep_family,
                                               i_cns_start, i_cns_end)
                m_hits_by_site[ins_chrm][ins_pos].append(hit_rcd)

                if ins_chrm not in m_hits_by_ctg:
                    m_hits_by_ctg[ins_chrm]={}
                if ins_pos not in m_hits_by_ctg[ins_chrm]:
                    m_hits_by_ctg[ins_chrm][ins_pos]={}

                if s_contig not in m_hits_by_ctg[ins_chrm][ins_pos]:
                    m_hits_by_ctg[ins_chrm][ins_pos][s_contig]=[]
                m_hits_by_ctg[ins_chrm][ins_pos][s_contig].append(hit_rcd)
        return m_hits_by_site, m_hits_by_ctg
####
####
    ####
    ####1000   18.8  5.2 11.1  10_110453470_ctg1      1   445  (828) + L1M4c     LINE/L1            737   1246 (6271)   1 *
    def slct_copies_with_min_len(self, i_min_len, b_tmplt_with_chr):
        m_copies={}
        with open(self.sf_rmsk) as fin_rmsk:
            for line in fin_rmsk:
                fields=line.split()
                if len(fields)<1:
                    continue
                if fields[0]=="score" or fields[0]=="SW":
                    continue
                ori_chrm=fields[4]
                istart=int(fields[5])
                iend=int(fields[6])
                if abs(iend-istart)<i_min_len:
                    continue
                ####
                chrm = self._process_chrm_name(b_tmplt_with_chr, ori_chrm)
                b_rc=False
                if fields[8]=="C":
                    b_rc=True
                if chrm not in m_copies:
                    m_copies[chrm] = {}
                m_copies[chrm][istart] = (iend, b_rc)
        return m_copies

####
## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, b_tmplt_with_chr, chrm):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True
        # print chrm, self.b_with_chr, b_chrm_with_chr #################################################################

        if b_tmplt_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_tmplt_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_tmplt_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm
####