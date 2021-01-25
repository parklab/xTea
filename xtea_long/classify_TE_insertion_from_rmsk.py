##12/17/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

from rmsk_TE_ins_classifer import *
from optparse import OptionParser
####
####
def parse_option():#
    parser = OptionParser()
    parser.add_option("-H", "--herv",
                      action="store_true", dest="herv", default=False,
                      help="For HERV only")
    parser.add_option("-R", "--ref_sva",
                      action="store_true", dest="ref_sva", default=False,
                      help="For reference SVA copies")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")
    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("--rmsk", dest="rmsk",
                      help="The repeatmasker annotation file ", metavar="FILE")
    parser.add_option("--hexam", dest="hexam",
                      help="The hexamer list file ", metavar="FILE")
    parser.add_option("--mast2", dest="mast2",
                      help="The mast2 list file ", metavar="FILE")
    parser.add_option("--ref", dest="reference",
                      help="The mast2 list file ", metavar="FILE")
    parser.add_option("--seq", dest="seq", default="",#the sequence file
                      help="The fasta file ", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    (options, args) = parser.parse_args()
    return (options, args)

####
if __name__ == '__main__':
    (options, args) = parse_option()
    b_herv=options.herv
    b_ref_sva=options.ref_sva

    sf_wfolder=options.wfolder
    n_jobs=options.cores
    sf_rmsk=options.rmsk
    sf_mast2=options.mast2
    sf_hexamer=options.hexam
    sva_clsfer=SVAClassifier(sf_hexamer, sf_wfolder, n_jobs)
    sf_out=options.output
    sf_seq = options.seq  #the assembled contig
    sf_reference=options.reference
    m_contigs={}#
    if os.path.isfile(sf_seq)==True:
        m_contigs=sva_clsfer.load_in_contigs(sf_seq)

    if b_herv==False and b_ref_sva==False:
        f_mast2_cutoff=0.7#at least 75% of the sequence (either TD or internal segment) is algned
        #sf_out_sva=sf_out+"_SVA.txt"
        f_hit_cutoff_sva=0.6#SVA at least take 50% of all the sequences
        m_sva=sva_clsfer.classify_SVA_from_rmsk_output(sf_rmsk, m_contigs, sf_mast2, sf_reference, f_hit_cutoff_sva,
                                                       f_mast2_cutoff, sf_out)

        alu_clsfer=AluClassifier(sf_wfolder, n_jobs)#
        f_hit_cutoff=0.85
        i_min_len=100
        i_max_len=400
        sf_out_Alu=sf_out+"_Alu.txt"
        b_sample_id=True
####note need to change this if needed
        #alu_clsfer.set_separator("~") #
        m_alu=alu_clsfer.classify_from_rmsk_output(sf_rmsk, m_contigs, m_sva, f_hit_cutoff, i_min_len, i_max_len,
                                                   sf_out_Alu, b_sample_id)

        l1_clsfer=LINE1Classifier(sf_wfolder, n_jobs)
        f_hit_cutoff = 0.55
        i_min_len = 150
        i_max_len = 6500
        sf_out_L1 = sf_out + "_L1_cutoff_0.55.txt"
####note need to change this if needed
        #l1_clsfer.set_separator("~")
        l1_clsfer.classify_from_rmsk_output(sf_rmsk, m_contigs, m_sva, f_hit_cutoff, i_min_len, i_max_len,
                                            sf_out_L1, b_sample_id)
    ####
    #### -R --ref {}
    elif b_ref_sva==True:
        i_max_dist=50#max distance between two records in rmsk annotation (if smaller than this, then view as one rcd)
        sva_clsfer.collect_full_length_ref_SVA(sf_rmsk, i_max_dist, sf_reference, sf_out)

    else:
        herv_clsfer=HERVClassifier(sf_wfolder, n_jobs)
        f_hit_cutoff = 0.85
        i_min_len = 150
        i_max_len = 9500
        sf_out_herv = sf_out + "_herv_cutoff_0.85.txt"
        m_already={}
        b_sample_id=False #not in format: HG002_chrm_pos
        herv_clsfer.classify_from_rmsk_output(sf_rmsk, m_contigs, m_already, f_hit_cutoff,
                                              i_min_len, i_max_len, sf_out_herv, b_sample_id)
####
    #check the transduction part, and the internal region
####
####