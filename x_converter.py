##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool

'''
Function: convert original bam to barcode indexed bam
1. save the "RNAME" field to "ch" tag
2. save the "RNEXT" field to "nh" tag
3. save the "POS" field to "cp" tag
4. set the "RNAME" field to the "BX" tag (barcode, with "-" changed to "_")
5. set the "POS" field to 0
6. set the "RNEXT" field to "*" (just to save space)
7. Reheader the bam (change the chromosomes to "barcodes"), to allow randome access
8. Sort and index the final bam
'''

SAMTOOLS_PATH = "samtools"  # by default, assume samtools is global available

def unwrap_self_header(arg, **kwarg):
    return BarcodeHeader._run_parse_header_by_chrom(*arg, **kwarg)


class BarcodeHeader():
    def __init__(self, sf_bam, working_folder, n_jobs):
        self._sf_bam = sf_bam
        self._working_folder = working_folder
        self._n_jobs = n_jobs

    def _run_parse_header_by_chrom(self, record):
        chrm = record[0]
        sf_bam = record[1]
        working_folder = record[2]

        samfile = pysam.AlignmentFile(sf_bam, "rb")
        set_barcode = set()
        for alignmt in samfile.fetch(chrm):
            if alignmt.has_tag("BX"):  # if there is a barcode
                s_barcode = alignmt.get_tag("BX")
                fields = s_barcode.split('-')
                s_new_barcode = "_".join(fields)
                set_barcode.add(s_new_barcode)
            else:  # if no barcode
                set_barcode.add("empty")
        samfile.close()

        # here also save the barcode into a file, prepare for creating the header
        # later, the barcodes will be collected to generate the final header
        bam_fields = sf_bam.split("/")
        sf_barcode_file = working_folder + str(bam_fields[-1]) + "_" + str(chrm) + ".barcode"
        with open(sf_barcode_file, "w") as fout_barcode:
            for barcode in set_barcode:
                fout_barcode.write(barcode + "\n")

    # Return the map of barcodes, and each record in format: {barcode: index}
    # Also the barcode are saved into the file
    def parse_barcode_from_alignment(self, sf_merge_barcode):
        samfile = pysam.AlignmentFile(self._sf_bam, "rb")
        references = samfile.references
        l_chrm_records = []

        for chrm in references:
            l_chrm_records.append((chrm, self._sf_bam, self._working_folder))
        samfile.close()

        pool = Pool(self._n_jobs)
        pool.map(unwrap_self_header, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        # merge the files
        cnt_barcode = 0
        bam_fields = self._sf_bam.split("/")
        m_barcodes = {}
        with open(sf_merge_barcode, "w") as fout_merged_barcode:
            for chrm in references:
                sf_barcode_file = self._working_folder + str(bam_fields[-1]) + "_" + str(chrm) + ".barcode"
                if os.path.isfile(sf_barcode_file) == False:
                    print "Error: Barcode file for chrom {0} doesn't exist!!!!".format(chrm)
                    continue
                with open(sf_barcode_file) as fin_barcode_file:
                    for line in fin_barcode_file:
                        sbarcode = line.rstrip()
                        if sbarcode not in m_barcodes:
                            m_barcodes[sbarcode] = cnt_barcode
                            fout_merged_barcode.write(sbarcode + "\t" + str(cnt_barcode) + "\n")
                            cnt_barcode += 1
                # remove the file
                os.remove(sf_barcode_file)

    def _get_old_bam_header(self):
        samfile = pysam.AlignmentFile(self._sf_bam, "rb")
        old_header = samfile.header  # it's a dictionary
        samfile.close()
        return old_header

    def gnrt_barcode_bam_header(self, sf_merge_barcode):
        old_header = self._get_old_bam_header()
        new_header = {}
        if "HD" in old_header:
            new_header["HD"] = old_header["HD"]
        if "PG" in old_header:
            new_header["PG"] = old_header["PG"]
        if "RG" in old_header:
            new_header["RG"] = old_header["RG"]
        new_header["SQ"] = []
        with open(sf_merge_barcode) as fin_barcode:
            for line in fin_barcode:
                fields = line.split()
                sbarcode = fields[0]
                # index=fields[1]
                m_tmp_barcode = {'LN': 1}
                m_tmp_barcode['SN'] = sbarcode
                new_header["SQ"].append(m_tmp_barcode)
        return new_header

    def _write_header_to_sam(self, header, sf_header_sam_file):
        with pysam.AlignmentFile(sf_header_sam_file, "w", header=header) as outf:  ##write the new header into file
            a = pysam.AlignedSegment()
            a.query_name = "read_28833_29006_6945_tmp"
            a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            a.flag = 99
            a.reference_id = 0
            a.reference_start = 32
            a.mapping_quality = 20
            a.cigar = ((0, 10), (2, 1), (0, 25))
            a.next_reference_id = 0
            a.next_reference_start = 199
            a.template_length = 167
            a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
            outf.write(a)

    def write_old_header_to_file(self, sf_header):
        old_header=self._get_old_bam_header()
        self._write_header_to_sam(old_header, sf_header)

##Function: pool.map doesnt' accept function in the class
##So, use this wrapper to solve this
def unwrap_self_converter(arg, **kwarg):
    return XConverter._run_cvt_bam_by_chrm(*arg, **kwarg)


class XConverter():
    # sf_bam: original input bam
    # sf_barcode_ba: final output (sorted and indexed)
    # n_jobs: number of cores
    def __init__(self, sf_bam, sf_barcode_bam, n_jobs):
        self._sf_bam = sf_bam
        self._sf_barcode_bam = sf_barcode_bam
        self._n_jobs = n_jobs
        self.working_folder = "./"

    def _gnrt_sam_header(self, bh, sf_barcode_file, sf_header_sam):
        bh.parse_barcode_from_alignment(sf_barcode_file)  # write the barcode and index into file
        new_header = bh.gnrt_barcode_bam_header(sf_barcode_file)
        bh._write_header_to_sam(new_header, sf_header_sam)

    def _run_cvt_bam_by_chrm(self, record):
        chrm = record[0]
        sf_bam = record[1]
        working_folder = record[2]
        sf_barcode_file = record[3]

        m_barcodes = {}  # dictionary of all barcodes
        with open(sf_barcode_file) as fin_barcode:  # load barcode from file
            for line in fin_barcode:
                fields = line.split()
                barcode = fields[0]
                index = fields[1]
                m_barcodes[barcode] = int(index)

        samfile = pysam.AlignmentFile(sf_bam, "rb")
        # get the reference id and name 1:1 map
        m_chrm = {}
        references = samfile.references
        for schrm in references:
            chrm_id = samfile.get_tid(schrm)
            m_chrm[chrm_id] = schrm
        m_chrm[-1] = "*"

        bam_fields = sf_bam.split("/")
        sf_cvted_sam = working_folder + bam_fields[-1] + "_" + chrm
        cvted_sam = pysam.AlignmentFile(sf_cvted_sam, "wb", template=samfile)
        # test_cnt=0
        for alignmt in samfile.fetch(chrm):
            s_new_barcode = "empty"
            if alignmt.has_tag("BX"):
                s_barcode = alignmt.get_tag("BX")
                fields = s_barcode.split('-')
                s_new_barcode = "_".join(fields)
                if s_new_barcode not in m_barcodes:
                    print "Error Happen!!! Barcode {0} doesnt' match the list".format(s_new_barcode)
            old_ref_id = alignmt.reference_id
            old_pos = alignmt.reference_start
            old_nref_id = alignmt.next_reference_id
            i_index = m_barcodes[s_new_barcode]

            a = pysam.AlignedSegment()  # new pysam.AlignedSegment()
            a.query_name = alignmt.query_name
            a.query_sequence = alignmt.query_sequence
            a.flag = alignmt.flag
            a.reference_id = i_index
            a.reference_start = 0
            a.mapping_quality = alignmt.mapping_quality
            a.cigar = alignmt.cigar
            a.next_reference_id = -1
            a.next_reference_start = alignmt.next_reference_start
            a.template_length = alignmt.template_length
            a.query_qualities = alignmt.query_qualities
            a.set_tags(alignmt.get_tags())
            a.set_tag("BX", None)
            a.set_tag("ch", m_chrm[old_ref_id], 'Z') ##original mapped chromosome
            a.set_tag("cp", old_pos, 'i') ##original mapped position
            a.set_tag("nh", m_chrm[old_nref_id], 'Z') ##original mate mapped chromosome
            cvted_sam.write(a)  ##perhaps use a buffer to speed up

            # alignmt.reference_id = i_index
            # alignmt.reference_start = 0
            # alignmt.next_reference_id=-1  ############??????????????????????????????
            # alignmt.set_tag("BX", None)
            # alignmt.set_tag("ch", m_chrm[old_ref_id], 'Z')
            # alignmt.set_tag("cp", old_pos, 'i')
            # alignmt.set_tag("nh", m_chrm[old_nref_id], 'Z')
            # cvted_sam.write(alignmt) ##perhaps use a buffer to speed up

            # test_cnt += 1  ###############################################################################################3
            # if test_cnt > 1000:
            #     break
        cvted_sam.close()
        samfile.close()

    def cvt_bam_to_barcode_bam(self):
        # open bam, read one by one
        # convert to barcode index, and then write to the new bam
        bh = BarcodeHeader(self._sf_bam, self.working_folder, self._n_jobs)
        sf_barcode_file = self.working_folder + "merged_barcodes.txt"  # this file to save all the barcodes
        sf_header_sam = self.working_folder + "header.sam"  ##this file to save the header bam file
        self._gnrt_sam_header(bh, sf_barcode_file, sf_header_sam)

        ##first run this chrom by chrom
        # get the head of the alignment file
        samfile = pysam.AlignmentFile(self._sf_bam, "rb")
        references = samfile.references
        l_chrm_records = []
        for chrm in references:
            l_chrm_records.append((chrm, self._sf_bam, self.working_folder, sf_barcode_file))
        samfile.close()

        pool = Pool(self._n_jobs)
        pool.map(unwrap_self_converter, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        ##then merge and reheader the bam files
        sf_merged_reheaded_bam = self.working_folder + "merged_reheaded.bam"
        cmd = "{0} cat -h {1} -o {2} ".format(SAMTOOLS_PATH, sf_header_sam, sf_merged_reheaded_bam)
        bam_fields = self._sf_bam.split("/")
        for chrm in references:
            sf_cvted_sam = self.working_folder + bam_fields[-1] + "_" + chrm
            cmd += "{0} ".format(sf_cvted_sam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

        # remove the temporary bam files
        for chrm in references:
            sf_cvted_sam = self.working_folder + bam_fields[-1] + "_" + chrm
            os.remove(sf_cvted_sam)

        # sort and index
        #pysam.sort("-m", "1G", "-o", self._sf_barcode_bam, sf_merged_reheaded_bam)
        #pysam.index(self._sf_barcode_bam)
        ##clean the temporary files
        os.remove(sf_header_sam)
        # os.remove(sf_merged_reheaded_bam)
        return

    def set_working_folder(self, working_folder):
        self.working_folder = working_folder
        if self.working_folder[-1] != "/":
            self.working_folder += "/"

    def set_samtools_path(self, samtools_path):
        global SAMTOOLS_PATH
        SAMTOOLS_PATH = samtools_path
