##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
##This code is used to merge the bam files of the sam individuals
##In order to avoid the mixture of barcodes, we

import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool

SAMBAMBA_PATH = "sambamba"


def unwrap_self_merge(arg, **kwarg):
    return BamMerger._run_merge_bam_chrm(*arg, **kwarg)


class BamMerger():
    def __init__(self, sf_file_list, sf_head_bam, n_jobs, working_folder, sf_merged):
        self.sf_file_list = sf_file_list
        self.sf_merged = sf_merged
        self.sf_head_bam = sf_head_bam
        self.n_jobs = n_jobs
        self.working_folder = working_folder
        if self.working_folder[-1] != "/":
            self.working_folder += "/"

    # merge two bam files (aligned from long ranger)
    # Need to change the "BX" tage for each alignment, will skip the first one
    def _run_merge_bam_chrm(self, chrm_record):
        chrm = chrm_record[0]
        sf_head_bam = chrm_record[1]
        sf_file_list = chrm_record[2]
        working_folder = chrm_record[3]
        sf_merged = working_folder + "tempfile_" + str(chrm)

        l_bams = []
        with open(sf_file_list) as fin_list:
            for line in fin_list:
                l_bams.append(line.rstrip())

        head_bam = pysam.AlignmentFile(sf_head_bam, "rb")
        merged_bam = pysam.AlignmentFile(sf_merged, "wb", template=head_bam)
        b_first_bam = True
        i_sample = 2
        for sf_bam in l_bams:
            samfile = pysam.AlignmentFile(sf_bam, "rb")
            for alignmt in samfile.fetch(chrm):
                if b_first_bam == True:
                    merged_bam.write(alignmt)
                else:
                    ##change the barcode
                    new_bx = "empty"
                    if alignmt.has_tag("BX"):
                        old_bx = alignmt.get_tag('BX')  ##e.g. GCGGTGATGGGGAAAA-1
                        fields = old_bx.split("-")
                        new_bx = fields[0] + "-" + str(i_sample)
                    alignmt.set_tag('BX', new_bx)
                    merged_bam.write(alignmt)
            b_first_bam = False
            i_sample += 1
            samfile.close()
        merged_bam.close()
        head_bam.close()

        self.sort_index_bam_with_sambamba(sf_merged, working_folder)


    def merge_bam(self):
        samfile = pysam.AlignmentFile(self.sf_head_bam, "rb")
        references = samfile.references
        l_chrm_records = []

        for chrm in references:
            l_chrm_records.append((chrm, self.sf_head_bam, self.sf_file_list, self.working_folder))
        samfile.close()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_merge, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        cmd = "{0} merge -t {1} {2} ".format(SAMBAMBA_PATH, self.n_jobs, self.sf_merged)
        ##merge the temp files to generate the final merged bam
        for record in l_chrm_records:
            chrm = record[0]
            sf_tmp_chrm = self.working_folder + "tempfile_" + str(chrm)+".sorted.bam"
            if os.path.isfile(sf_tmp_chrm):
                cmd = cmd + sf_tmp_chrm + " "
        Popen(cmd, shell=True, stdout=PIPE).communicate()

        ##delete the temp files
        for record in l_chrm_records:
            chrm = record[0]
            sf_tmp_chrm = self.working_folder + "tempfile_" + str(chrm)+".sorted.bam"
            if os.path.isfile(sf_tmp_chrm):
                os.remove(sf_tmp_chrm) ###################################################################################


    # sort and index the merged bam
    def sort_index_bam_with_sambamba(self, sf_bam, sf_tmp_folder):
        cmd = "{0} sort -m 10G --tmpdir={1} {2}".format(SAMBAMBA_PATH, sf_tmp_folder, sf_bam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()
        os.remove(sf_bam)  ###################################################################################
