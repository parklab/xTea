##11/01/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

####canonical signals: AATAAA / TTTATT, ATTAAA/TTTAAT
####variations: AGTAAA, TATAAA, CATAAA, GATAAA, AATATA, AATACA, AATAGA, AAAAAG, ACTAAA, AAGAAA, AATGAA, TTTAAA, AAAACA,
# GGGGCT
import re

class PolyA():
    ####
    def is_poly_A_T(self, seq):  ###for a temp version here
        cnt_A = 0
        cnt_T = 0
        for ch in seq:
            if ch == 'A' or ch == 'a':
                cnt_A += 1
            elif ch == 'T' or ch == 't':
                cnt_T += 1
        cnt_cutoff = (len(seq) * 0.75)
        if (cnt_A > cnt_cutoff) or (cnt_T > cnt_cutoff):
            return True
        return False

    def contain_poly_A_T(self, seq, n_min_cnt):
        max_A = 0
        max_T = 0
        cum_A = 0
        cum_T = 0
        for ch in seq:
            if ch == 'A' or ch == 'a':
                cum_A += 1
            else:
                if cum_A > max_A:
                    max_A = cum_A
                cum_A = 0

            if ch == 'T' or ch == 't':
                cum_T += 1
            else:
                if cum_T > max_T:
                    max_T = cum_T
                cum_T = 0
        if cum_A > max_A:
            max_A = cum_A
        if cum_T > max_T:
            max_T = cum_T

        if max_A >= n_min_cnt or max_T >= n_min_cnt:
            return True
        else:
            return False

    def is_consecutive_polyA_T(self, seq):
        if ("AAAAA" in seq) or ("TTTTT" in seq) or ("AATAA" in seq) or ("TTATT" in seq):
            return True
        else:
            return False
            ####

    def is_consecutive_polyA(self, seq):
        if ("AAAAA" in seq) or ("AATAA" in seq):
            return True
        else:
            return False

    def is_consecutive_polyA_T_with_ori(self, seq):
        b_polyA=True
        if ("AAAAA" in seq) or  ("AATAA" in seq) :#polyA
            return True, b_polyA
        elif ("TTTTT" in seq) or ("TTATT" in seq):#polyT
            b_polyA=False
            return True, b_polyA
        else:
            return False, False

    #this is not used for now, but should be used to replace "is_consecutive_polyA_T"
    #Right clipped reads should be polyT only, while left clipped reads should be polyA only.
    def is_consecutive_polyA_T_with_oritation(self, seq, b_left):
        if b_left==True:#left clip
            if ("AAAAA" in seq) or ("AATAA" in seq):
                return True
        else:
            if ("TTTTT" in seq) or ("TTATT" in seq):
                return True
        # if ("AAAAA" in seq) or ("TTTTT" in seq) or ("AATAA" in seq) or ("TTATT" in seq):
        #     return True
        # else:
        return False

    #check whether contain enough A or T
    def contain_enough_A_T(self, seq, n_min_A_T):
        n_A=0
        n_T=0
        for ch in seq:
            if ch == 'A' or ch == 'a':
                n_A+=1
            if ch == 'T' or ch == 't':
                n_T+=1

        if n_A>=n_min_A_T or n_T>=n_min_A_T:
            return True
        return False
