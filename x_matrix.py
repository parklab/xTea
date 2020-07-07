import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
from x_barcode import XBarcode

class XMatrix():
    def __init__(self, sf_bam, sf_matrix):
        self.sf_bam=sf_bam
        self.sf_matrix=sf_matrix

    def gnrt_matrix(self):
        #TBD
        return