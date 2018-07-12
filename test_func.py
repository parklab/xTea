import os
import sys
from optparse import OptionParser
import pysam
# from x_filter import *
#
# sf_test=sys.argv[1]
# xfilter = XFilter()
# m_original_sites = xfilter.load_in_candidate_list(sf_test)
# sf_peak_sites = "test_clip_peak_candidate.list"
# m_sites_clip_peak = xfilter.call_peak_candidate_sites(m_original_sites)  # get the peak sites
# xfilter.output_candidate_sites(m_sites_clip_peak, sf_peak_sites)  # output the sites

# sf_sam=sys.argv[1]
# sf_out=sys.argv[2]
# with open(sf_out,"w") as fout_fa:
#     samfile = pysam.AlignmentFile(sf_sam, "r")
#     for algnmt in samfile.fetch():
#         l_cigar=algnmt.cigar
#         cnt_m=0
#         cnt_total=0
#         for grp in l_cigar:
#             if int(grp[0])==0:
#                 cnt_m+=int(grp[1])
#             cnt_total+=int(grp[1])
#         if cnt_total==0:
#             continue
#         if float(cnt_m)/float(cnt_total) < 0.95:
#             continue
#         qname=algnmt.query_name
#         seq=algnmt.query_sequence
#         fout_fa.write(">"+qname+"\n")
#         fout_fa.write(seq + "\n")


# extract the reads by read names
def extract_reads_by_name_list(sf_bam):
    #l_names = self.load_read_names(sf_names)
    bamfile = pysam.AlignmentFile(sf_bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)  # here use hashing to save the read names in the memory
    name_indexed.build()



#
# def test_index_name(sf_bam):
#     bamfile = pysam.AlignmentFile(sf_bam, 'rb')
#     name_indexed = pysam.IndexedReads(bamfile)  # here use hashing to save the read names in the memory
#     name_indexed.build()

def check_process_asm_contig(sf_contig):
    b_valid = True
    if os.path.isfile(sf_contig) == True and os.stat(sf_contig).st_size > 0:
        with open(sf_contig, "rb") as f:
            f.seek(-2, os.SEEK_END)  # Jump to the second last byte.
            while f.read(1) != b"\n":  # Until EOL is found...
                f.seek(-2, os.SEEK_CUR)  # ...jump back the read byte plus one more.
            last = f.readline()  # Read last line.
            if len(last) > 0 and last[0] == ">":
                b_valid = False

    if b_valid == False:
        lines = []
        with open(sf_contig) as fin_contig:
            lines = fin_contig.readlines()
        with open(sf_contig, "w") as fout_contig:
            fout_contig.writelines([item for item in lines[:-1]])

def parse_option():
    parser = OptionParser()
    parser.add_option("--Insert",
                      action="store_true", dest="insert", default=False,
                      help="calculate insert size")
    parser.add_option("-C", "--clip",
                      action="store_true", dest="clip", default=False,
                      help="Call candidate TEI sites from clipped reads")
    parser.add_option("-D", "--discordant",
                      action="store_true", dest="discordant", default=False,
                      help="Filter with discordant paired end reads")

    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("-r", "--reference", dest="reference",
                      help="The reference file ", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")

    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == '__main__':
    # ss=sys.argv[1]
    # check_process_asm_contig(ss)
    # (options, args) = parse_option()
    # if options.insert:
    #     print "insert I1"

    sf_bam=sys.argv[1]
    extract_reads_by_name_list(sf_bam)
    