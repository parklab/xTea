##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
from x_converter import *
from x_merger import *
from x_trimmer import *
from optparse import OptionParser

####
##parse the options
def parse_option():
    parser = OptionParser()
    parser.add_option("-c", "--convert",
                      action="store_true", dest="convert", default=False,
                      help="convert the bam to barcode indexed bam")
    parser.add_option("-m", "--merge",
                      action="store_true", dest="merge", default=False,
                      help="merge 10X bams")
    parser.add_option("-s", "--statistic",
                      action="store_true", dest="statistic", default=False,
                      help="basic statistic from 10X bam")
    parser.add_option("-t", "--trim",
                      action="store_true", dest="trim", default=False,
                      help="Trim the clipped reads caused by errors from 10X bam")
    parser.add_option("-i", "--input", dest="input",
                      help="input file) ", metavar="FILE")
    parser.add_option("-r", "--reference", dest="reference",
                      help="The reference file) ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")

    (options, args) = parser.parse_args()
    return (options, args)


##main function
if __name__ == '__main__':
    (options, args) = parse_option()
    if options.convert:
        sf_ori_bam = options.input
        sf_barcode_bam = options.output
        s_working_folder = options.wfolder
        n_jobs = options.cores
        xcvtr = XConverter(sf_ori_bam, sf_barcode_bam, n_jobs)
        xcvtr.set_working_folder(s_working_folder)
        xcvtr.cvt_bam_to_barcode_bam()
    elif options.merge:
        sf_file_list= options.input
        sf_head_bam=options.bam
        sf_merged=options.output
        s_tmp_folder = options.wfolder
        n_jobs=options.cores

        bam_merger=BamMerger(sf_file_list, sf_head_bam, n_jobs, s_tmp_folder, sf_merged)
        bam_merger.merge_bam()
        #bam_merger.sort_index_bam_with_sambamba(sf_merged, s_tmp_folder)
    elif options.trim:
        sf_input_bam=options.bam
        sf_ref=options.reference
        s_tmp_folder = options.wfolder
        n_jobs = options.cores
        sf_trimmed_bam=options.output
        rt = ReadTrimmer(sf_input_bam, n_jobs)
        rt.trim_reads(s_tmp_folder, sf_ref, sf_trimmed_bam)
##TDList:
##1. given a region, get out how many barcodes inside
##2. given a site, get out how many barcodes cover
##3. 