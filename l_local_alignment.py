##03/14/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

#The code is revised from an implementation released here:
#https://nbviewer.jupyter.org/github/BenLangmead/comp-genomics-class/blob/master/notebooks/CG_DP_Local.ipynb
import numpy

class Local_alignment():
    def exampleCost(self, xc, yc):
        ''' Cost function: 2 to match, -6 to gap, -4 to mismatch '''
        if xc == yc: return 2  # match
        if xc == '-' or yc == '-': return -6  # gap
        return -4

    def smithWaterman(self, x, y, s):
        ''' Calculate local alignment values of sequences x and y using
            dynamic programming.  Return maximal local alignment value. '''
        V = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                V[i, j] = max(V[i - 1, j - 1] + s(x[i - 1], y[j - 1]),  # diagonal
                              V[i - 1, j] + s(x[i - 1], '-'),  # vertical
                              V[i, j - 1] + s('-', y[j - 1]),  # horizontal
                              0)  # empty
        #argmax = numpy.where(V == V.max())
        #return V, int(V[argmax])
        return V

    def traceback(self, V, x, y, s):
        """ Trace back from given cell in local-alignment matrix V """
        # get i, j for maximal cell
        i, j = numpy.unravel_index(numpy.argmax(V), V.shape)
        xscript, alx, aly, alm = [], [], [], []
        while (i > 0 or j > 0) and V[i, j] != 0:
            diag, vert, horz = 0, 0, 0
            if i > 0 and j > 0:
                diag = V[i - 1, j - 1] + s(x[i - 1], y[j - 1])
            if i > 0:
                vert = V[i - 1, j] + s(x[i - 1], '-')
            if j > 0:
                horz = V[i, j - 1] + s('-', y[j - 1])
            if diag >= vert and diag >= horz:
                match = x[i - 1] == y[j - 1]
                xscript.append('M' if match else 'R')
                alm.append('|' if match else ' ')
                alx.append(x[i - 1])
                aly.append(y[j - 1])
                i -= 1
                j -= 1
            elif vert >= horz:
                xscript.append('D')
                alx.append(x[i - 1])
                aly.append('-')
                alm.append(' ')
                i -= 1
            else:
                xscript.append('I')
                aly.append(y[j - 1])
                alx.append('-')
                alm.append(' ')
                j -= 1
        xscript = (''.join(xscript))[::-1]
        seqs = ''.join(aly[::-1])
        # alignment = '\n'.join(map(lambda x: ''.join(x), [alx[::-1], alm[::-1], aly[::-1]]))
        return xscript, seqs

    def editDistanceLikeCost(self, xc, yc):
        return 1 if xc == yc else -1

    def is_seqs_matched(self, s_q, s_t, f_ratio):
        V = self.smithWaterman(s_q, s_t, self.editDistanceLikeCost)
        xscript, seqs = self.traceback(V, s_q, s_t, self.exampleCost)
        if len(s_t) == 0:
            return False
        if float(len(seqs)) / float(len(s_t)) > f_ratio:
            return True
        return False
####