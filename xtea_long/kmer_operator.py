##12/16/2020
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu or chong.simon.chu@gmail.com

class KmerOperator():
    def __init__(self, i_k):##
        self.k_len=i_k
        self.m_alpha_rc={}
        self.m_alpha = {}#
        self._set_rc_alpha_map()
        self._set_upper_case_map()

    def set_k(self, i_k):
        self.k_len=i_k

    def gnrt_reverse_complementary(self, s_seq):
        s_rc = ""
        for s in s_seq[::-1]:
            if s not in self.m_alpha_rc:
                s_rc += "N"
            else:
                s_rc += self.m_alpha_rc[s]
        return s_rc

    #convert the seq to upper case format
    def cvt_to_upper_case_seq(self, s_seq):
        s_new_seq=""
        for s in s_seq:
            if s not in self.m_alpha:
                s_new_seq+="N"
            else:
                s_new_seq+=self.m_alpha[s]
        return s_new_seq
    ####
    ####generate one institution distance
    def gnrt_one_substitution_seqs(self, s_seq, b_rc=True):
        m_seqs={}
        l_alpha=["A","C","G","T"]
        n_len=len(s_seq)
        for i in range(n_len):
            for c_new in l_alpha:
                s_seq_new=""
                if i==0:
                    s_seq_new=c_new+s_seq[1:]
                elif i==n_len-1:
                    s_seq_new=s_seq[0:n_len-1]+c_new
                else:
                    s_seq_new=s_seq[0:i]+c_new+s_seq[i+1:]
                m_seqs[s_seq_new] = 1
                if b_rc==True:
                    s_seq_new_rc=self.gnrt_reverse_complementary(s_seq_new)
                    m_seqs[s_seq_new_rc]=1
        return m_seqs

####
    ####generate sequence with indel distance of 1
    ####s_seq is the orginal sequence
    def gnrt_one_insertion_distance_seqs(self, s_seq):
        m_seqs = {}
        l_alpha = ["A", "C", "G", "T"]  #
        n_len = len(s_seq)
        for i in range(n_len):
            for c_new in l_alpha:
                s_new_seq=""
                if i==0:
                    s_new_seq=c_new+s_seq
                elif i==n_len-1:
                    s_new_seq=s_seq+c_new
                else:
                    s_new_seq=s_seq[:i]+c_new+s_seq[i:]
                m_seqs[s_new_seq]=1

        return m_seqs

    ####
    def gnrt_one_deletion_distance_seqs(self, s_seq):
        m_seqs = {}
        n_len = len(s_seq)
        for i in range(n_len):##
            s_new_seq = ""
            if i == 0:
                s_new_seq = s_seq[1:]
            elif i == n_len - 1:
                s_new_seq = s_seq[0:n_len-1]
            else:
                s_new_seq = s_seq[:i] + s_seq[i+1:]
            m_seqs[s_new_seq] = 1
        return m_seqs

    ####
    def _set_upper_case_map(self):
        self.m_alpha['A'] = 'A'
        self.m_alpha['a'] = 'A'
        self.m_alpha['T'] = 'T'
        self.m_alpha['t'] = 'T'
        self.m_alpha['C'] = 'C'
        self.m_alpha['c'] = 'C'
        self.m_alpha['G'] = 'G'
        self.m_alpha['g'] = 'G'
        self.m_alpha['N'] = 'N'
        self.m_alpha['n'] = 'N'

    ####
    def _set_rc_alpha_map(self):
        self.m_alpha_rc['A'] = 'T'
        self.m_alpha_rc['a'] = 'T'
        self.m_alpha_rc['T'] = 'A'
        self.m_alpha_rc['t'] = 'A'
        self.m_alpha_rc['C'] = 'G'
        self.m_alpha_rc['c'] = 'G'
        self.m_alpha_rc['G'] = 'C'
        self.m_alpha_rc['g'] = 'C'
        self.m_alpha_rc['N'] = 'N'
        self.m_alpha_rc['n'] = 'N'