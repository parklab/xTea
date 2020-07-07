##12/20/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu, chong.simon.chu@gmail.com

# this module is designed to call out the dimorphic HERV copies.
# This is a supplementary module of the "l_complex_sv.py" module, as most of the full-length HERV copy are long
# thus, it is hard to cover the whole insertion with the two flanking regions in one contig (maybe >10K)
# Here, we use two extra steps to check this:
# 1. Use the reference annotation, and find out all the "LTR" solo copies, and do intersection
# with the collected breakpoints (from the TE-complex-sv step);
# 2. If both breakpoints found, then check the previous assembled and flank trimmed seq.
import global_values
from optparse import OptionParser
from l_contig import *
from rmsk_TE_ins_classifer import *

####
class DimorphicHERV():
    def __init__(self):
        self.i_herv_tail = 9500
        self.i_min_herv_len = 150
        self.HERV_ratio=0.05
        self.s_l_suffix="L"
        self.s_r_suffix="R"
        self.s_LTR="LTR"
        self.s_HERV="HERV"

    '''
    #example annotation type:
    667   3.9  0.0  0.0  chr1        46416   46493 (248909929) C  LTR12F         LTR/ERV1             (115)  404    327     53
    350  31.0  8.0  0.0  chr1        70504   70712 (248885710) C  LTR89          LTR/ERVL?            (259)  620    336     75
    '''

    def slct_LTR_from_rmsk(self, sf_anno, sf_lbrkpnts, sf_rbrkpnts, i_win):
        m_lbrkpnts = self.load_in_brkpnts(sf_lbrkpnts)
        m_rbrkpnts = self.load_in_brkpnts(sf_rbrkpnts)
        m_slcted = {}
        with open(sf_anno) as fin_rmsk:
            for line in fin_rmsk:
                fields = line.split()
                if len(fields) < 1:
                    continue
                if fields[0] == "score" or fields[0] == "SW":
                    continue
                    ####
                ####potential issue here: chrm may not match!!!!!
                chrm = fields[4]
                i_start = int(fields[5])
                i_end = int(fields[6])
                s_rc = fields[8]
                s_sub_family = fields[9]
                s_rep_family = fields[10]
                # if self.s_LTR not in s_sub_family:
                #    continue
                if ("int" in s_sub_family) or ("ERV" in s_sub_family):
                    continue
                i_cns_start = -1
                i_cns_end = -1
                b_rc = False
                if s_rc == "+":
                    i_cns_start = int(fields[11])
                    i_cns_end = int(fields[12])
                else:
                    i_cns_start = int(fields[13])
                    i_cns_end = int(fields[12])
                    b_rc = True

                # check whether hit the brkpoints
                # 1. check whether i_start in rbrkpnts
                b_lhit, i_r_brkpnt = self.find_a_hit(chrm, i_start, i_win, m_rbrkpnts)
                if b_lhit == False:
                    continue
                # 2. check whether i_end in lbrkpnts
                b_rhit, i_l_brkpnt = self.find_a_hit(chrm, i_end, i_win, m_lbrkpnts)
                if b_rhit == False:
                    continue

                s_id = chrm + global_values.SEPERATOR + str(i_l_brkpnt) + global_values.SEPERATOR + str(i_r_brkpnt)
                m_slcted[s_id] = (chrm, i_start, i_end, s_sub_family, s_rep_family, b_rc, i_cns_start, i_cns_end)
        return m_slcted

    ####
    ####
    ####Here sf_ltrimed indicates right-clipped reads
    def check_slct_seq_of_slcted_LTR(self, m_slcted, sf_ltrimed, sf_rtrimed, sf_sites, sf_out):
        lrd_ctg = LRD_Contig("./", 1)
        # load in the left trimmed seqs, each id in format: chrm~pos
        m_l_ctgs = lrd_ctg.load_fasta_to_dict(sf_ltrimed)  # left-breakpoint assembled contig (after trim left-flank)
        # load in the right trimmed seqs
        m_r_ctgs = lrd_ctg.load_fasta_to_dict(sf_rtrimed)  # right-breakpoints assembled contig (after trim right-flank)

        sf_fa=sf_out + ".fa"
        with open(sf_sites, "w") as fout_sites, open(sf_fa, "w") as fout_fa:
            for s_id in m_slcted:
                (chrm, i_rmsk_start, i_rmsk_end, s_sub_family, s_family, b_rc, i_cns_start, i_cns_end) = m_slcted[s_id]
                l_id_fields = s_id.split(global_values.SEPERATOR)
                i_lbrkpnt = int(l_id_fields[1])
                i_rbrkpnt = int(l_id_fields[2])

                s_l_id = chrm + global_values.SEPERATOR + str(i_lbrkpnt)
                if s_l_id not in m_l_ctgs:
                    continue
                s_r_id = chrm + global_values.SEPERATOR + str(i_rbrkpnt)
                if s_r_id not in m_r_ctgs:
                    continue
                s_l_seq = m_l_ctgs[s_l_id][:self.i_herv_tail]#for left breakpoint, take the front part
                s_r_seq = m_r_ctgs[s_r_id][-1*self.i_herv_tail:]#for right breakpoint, take the tail part

                # s_lr_id=chrm+global_values.SEPERATOR+str(i_lbrkpnt)+global_values.SEPERATOR+str(i_rbrkpnt)
                fout_sites.write(chrm + "\t" + str(i_lbrkpnt) + "\t" + str(i_rbrkpnt) + "\t" + s_sub_family + "\n")
                fout_fa.write(">" + s_l_id + global_values.SEPERATOR + self.s_l_suffix+"\n" + s_l_seq + "\n")
                fout_fa.write(">" + s_r_id + global_values.SEPERATOR + self.s_r_suffix+"\n" + s_r_seq + "\n")
        # run RepeatMasker on the seqs
        herv_clsfer = HERVClassifier("./", 1)
        herv_clsfer.run_RMSK(sf_fa)
        sf_fa_rmsk = sf_fa + ".out"
        m_all_ctgs = lrd_ctg.load_fasta_to_dict(sf_fa)
        f_hit_cutoff = self.HERV_ratio
        i_min_len = self.i_min_herv_len
        i_max_len = self.i_herv_tail
        m_candidates = herv_clsfer.classify_from_rmsk_output(sf_fa_rmsk, m_all_ctgs, {}, f_hit_cutoff, i_min_len,
                                                         i_max_len, sf_out, self.s_HERV)
        return m_candidates
####
    ####
    ####
    # Here, for right breakpoints, the tail assembled part should have HERV seq
    # for left breakpoints, the front assembled part should also have HERV seq
    def slct_final_dimorphic_HERV_from_rmsk(self, sf_sites, m_slcted, sf_out):
        # load in the sites (which will pair the brkpnts)
        m_brkpnt_pairs = {}
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields = line.split()
                chrm = fields[0]
                if chrm not in m_brkpnt_pairs:
                    m_brkpnt_pairs[chrm] = {}
                i_lbrk = int(fields[1])
                i_rbrk = int(fields[2])
                s_sub_family=fields[3]
                m_brkpnt_pairs[chrm][i_lbrk] = (i_rbrk,s_sub_family)

        #print m_brkpnt_pairs,"pairs"

        m_left={}#left breakpoints (right clip indicates "left" brkpnts)
        m_right={}
        for s_id in m_slcted:#s_id in format: chrm~pos~L/R
            (chrm, pos, s_sub_family, cns_start, cns_end, seq_start, seq_end, s_contig_seq)=m_slcted[s_id]
            if s_id[-1] == self.s_l_suffix:
                if chrm not in m_left:
                    m_left[chrm]={}
                m_left[chrm][pos]=(chrm, pos, s_sub_family, cns_start, cns_end, seq_start, seq_end, s_contig_seq)
            else:
                if chrm not in m_right:
                    m_right[chrm]={}
                m_right[chrm][pos]=(chrm, pos, s_sub_family, cns_start, cns_end, seq_start, seq_end, s_contig_seq)

        #print m_left,"left"
        #print m_right,"right"

        with open(sf_out, "w") as fout:
            for chrm in m_left:
                if chrm not in m_brkpnt_pairs:
                    continue
                if chrm not in m_right:
                    continue
                for i_lbrk in m_left[chrm]:
                    l_l_rcd=m_left[chrm][i_lbrk]
                    if i_lbrk not in m_brkpnt_pairs[chrm]:
                        continue
                    i_l_seq_start=l_l_rcd[5]
                    i_l_seq_end=l_l_rcd[6]
                    s_l_contig_seq = l_l_rcd[-1][i_l_seq_start:i_l_seq_end]

                    i_rbrk=m_brkpnt_pairs[chrm][i_lbrk][0]
                    s_ltr_sub_family=m_brkpnt_pairs[chrm][i_lbrk][1]

                    if i_rbrk not in m_right[chrm]:
                        continue
                    l_r_rcd=m_right[chrm][i_rbrk]
                    i_r_seq_start = l_r_rcd[5]
                    i_r_seq_end = l_r_rcd[6]
                    s_r_contig_seq = l_r_rcd[-1][i_r_seq_start:i_r_seq_end]

                    s_contig_seq=s_l_contig_seq
                    s_sub_family=l_l_rcd[2]
                    cns_start=l_l_rcd[3]
                    cns_end=l_l_rcd[4]
                    seq_start=l_l_rcd[5]
                    seq_end=l_l_rcd[6]
                    if len(s_r_contig_seq)>len(s_l_contig_seq):
                        s_contig_seq = s_r_contig_seq
                        s_sub_family = l_r_rcd[2]
                        cns_start = l_r_rcd[3]
                        cns_end = l_r_rcd[4]
                        seq_start = l_r_rcd[5]
                        seq_end = l_r_rcd[6]

                    fout.write(chrm + "\t" + str(i_rbrk) + "\t.\t.\t" + s_contig_seq +
                               "\t60\t.\tEND="+str(i_lbrk)+";SUB_FAMILY=" + s_sub_family +
                               ";LTR_SUB_FAMILY=" + s_ltr_sub_family +
                               ";CNS=" + str(cns_start) + "-" + str(cns_end) +
                               ";SEQ_REGION=" + str(seq_start) + "-" + str(seq_end) + "\n")
####

    def load_in_brkpnts(self, sf_brkpnts):
        m_brkpnts = {}
        with open(sf_brkpnts) as fin_brkpnts:
            for line in fin_brkpnts:
                fields = line.rstrip().split()
                chrm = fields[0]
                if len(chrm) < 3 or chrm[:3] != "chr":
                    chrm = "chr" + chrm
                pos = int(fields[1])

                if chrm not in m_brkpnts:
                    m_brkpnts[chrm] = {}
                m_brkpnts[chrm][pos] = 1
        return m_brkpnts

    def find_a_hit(self, chrm, pos, i_win, m_brkpnts):
        if chrm not in m_brkpnts:
            return False, pos
        i_tmp_start = pos - i_win
        i_tmp_end = pos + i_win
        b_hit = False
        i_hit_pos = pos
        for i_tmp in range(i_tmp_start, i_tmp_end):
            if i_tmp in m_brkpnts[chrm]:
                b_hit = True
                i_hit_pos = i_tmp
                break
        return b_hit, i_hit_pos


def parse_option():
    parser = OptionParser()

    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-a", "--anno", dest="anno",
                      help="LTR annotation file", metavar="FILE")
    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("--i2", dest="input2", default="null",
                      help="input file ", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")
    (options, args) = parser.parse_args()
    return (options, args)


####
if __name__ == '__main__':
    (options, args) = parse_option()
    dm = DimorphicHERV()
    sf_anno = options.anno
    sf_wfolder = options.wfolder
    i_extend = -1000
    sf_out = options.output

    sf_asm_folder = sf_wfolder + "l_asm_tmp" + "/neg_{0}/".format(abs(i_extend))
    sf_l_trimmed = sf_asm_folder + "left_trimmed_contig.fa"  # trimmed segments from left flank alignments
    sf_r_trimmed = sf_asm_folder + "right_trimmed_contig.fa"  # trimmed segments from right flank alignments

    i_win = 150
    sf_lbrkpnts = options.input
    sf_rbrkpnts = options.input2

    m_slcted = dm.slct_LTR_from_rmsk(sf_anno, sf_lbrkpnts, sf_rbrkpnts, i_win)
    sf_sites = sf_out + ".sites"
    sf_tmp=sf_out+".tmp"
    m_candidates=dm.check_slct_seq_of_slcted_LTR(m_slcted, sf_l_trimmed, sf_r_trimmed, sf_sites, sf_tmp)
    dm.slct_final_dimorphic_HERV_from_rmsk(sf_sites, m_candidates, sf_out)
    ####
####