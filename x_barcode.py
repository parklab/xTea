##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import sys
import pysam



class XBarcode():
    def __init__(self, sf_bam, sf_barcode_bam):
        self.sf_bam = sf_bam
        self.sf_barcode_bam = sf_barcode_bam

    ##collect all the barcodes of a region
    def collect_barcode_in_region(self, chrm, istart, iend):
        set_barcode = set()
        samfile = pysam.AlignmentFile(self.sf_bam, "rb")
        for alignmt in samfile.fetch(chrm, istart, iend):
            if alignmt.has_tag("BX"):
                s_barcode = alignmt.get_tag("BX")
                set_barcode.add(s_barcode)
        samfile.close()

    def get_shared_barcode(self, set1, set2):
        return set1 & set2

    # for given barcode, get all the related alignments
    # filter out those unreleated to L1 repeats
    def _get_alignments_for_given_barcode(self, s_barcode):
        bamfile = pysam.AlignmentFile(self.sf_barcode_bam, "rb")
        iter_algnmts = bamfile.fetch(s_barcode)
        l_selected_algnmts = []
        for alnmt in iter_algnmts:  # each in "pysam.AlignedSegment" format
            if alnmt.is_duplicate == True or alnmt.is_supplementary == True:  ##skip duplicate and supplementary ones
                continue
            if alnmt.is_unmapped == True:  #### for now, just skip the unmapped reads ???????????????????????????????
                continue
            if alnmt.is_secondary == True:  ##skip secondary alignment
                continue
            l_selected_algnmts.append(alnmt)
        bamfile.close()
        return l_selected_algnmts

####