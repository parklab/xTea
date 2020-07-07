import pysam
import sys
from multiprocessing import Pool
from x_alignments import BamInfo
from x_annotation import *


####
def unwrap_self_extract_rna_reads_for_region(arg, **kwarg):
    return RNAignment.extract_reads_for_chrm(*arg, **kwarg)

class RNAignment():
    def __init__(self, sf_bam, sf_ref):
        self.sf_bam = sf_bam
        self.out_header = None
        self.chrm_id_name = {}
        self.sf_reference=sf_ref
        self.bam_info=BamInfo(sf_bam, sf_ref)

    def cnt_unique_mapped_reads_within_rep(self, sf_rmsk, sf_working_folder, n_jobs):

        m_chrms=self.bam_info.get_all_chrom_name_length()

        l_chrms=[]
        for chrm in m_chrms:
            record=(chrm, 0, m_chrms[chrm], sf_rmsk, sf_working_folder)
            l_chrms.append(record)
            #self.extract_reads_for_chrm(record)
        pool = Pool(n_jobs)
        pool.map(unwrap_self_extract_rna_reads_for_region, zip([self] * len(l_chrms), l_chrms), 1)
        pool.close()
        pool.join()

        #merge the results
        m_rep_cnt={}
        for rcd in l_chrms:
            chrm=rcd[0]
            sf_unq_cnt_chrm = sf_working_folder + "{0}_unique_cnt_within_rep.txt".format(chrm)
            with open(sf_unq_cnt_chrm) as fin_cnt:
                for line in fin_cnt:
                    fields=line.split()
                    sub_family=fields[0]
                    if sub_family not in m_rep_cnt:
                        m_rep_cnt[sub_family]=int(fields[1])
                    else:
                        m_rep_cnt[sub_family] += int(fields[1])

        sf_merged=sf_working_folder + "all_unique_cnt_within_rep.txt"
        with open(sf_merged, "w") as fout_merged:
            n_all=int(m_rep_cnt["all"])
            fout_merged.write("all\t{0}\n".format(m_rep_cnt["all"]))
            for sub_family in m_rep_cnt:
                if sub_family=="all":
                    continue
                n_sub_family=int(m_rep_cnt[sub_family])
                ratio=float(n_sub_family)/float(n_all)
                fout_merged.write("{0}\t{1}\t{2}\n".format(sub_family, n_sub_family, ratio))


    def extract_reads_for_chrm(self, record):
        chrm=record[0]
        sf_rmsk=record[3]
        sf_wfolder=record[4]

        x_annot = XAnnotation(sf_rmsk)
        x_annot.load_rmsk_annotation_with_divgnt_extnd(10)  # extend a little bit
        x_annot.index_rmsk_annotation()

        sf_unq_cnt_chrm=sf_wfolder+"{0}_unique_cnt_within_rep.txt".format(chrm)
        sf_L1HS_fa_chrm = sf_wfolder + "{0}_L1HS.fa".format(chrm)

        bamfile = pysam.AlignmentFile(self.sf_bam, 'rb', reference_filename=self.sf_reference)

        m_anno=x_annot.get_rmsk_annotation()
        n_all_reads=0
        with open(sf_unq_cnt_chrm, "w") as fout_unq, open(sf_L1HS_fa_chrm, "w") as fout_fa:
            m_rep_cnt={}#num of unique rna read within different rep families
            m_rep_cnt['L1HS_all']=0
            for alignment in bamfile.fetch(chrm):
                if alignment.is_duplicate == True or alignment.is_supplementary == True:  ##skip duplicate and supplemty ones
                    continue
                if alignment.is_unmapped == True:  #### for now, just skip the unmapped reads
                    continue
                if alignment.is_secondary == True:  ##skip secondary alignment
                    continue
                n_all_reads+=1

                chrm_rmsk="chr"+chrm
                pos=alignment.next_reference_start
                bhit, i_mid_pos = x_annot.is_within_repeat_region(chrm_rmsk, pos)
                if bhit==False:
                    continue

                read_id=alignment.query_name
                s_pair="_2"
                if alignment.is_read1 == True:
                    s_pair = "_1"
                seq=alignment.query_sequence

                sub_family=m_anno[chrm_rmsk][i_mid_pos][0][3]
                if sub_family=="L1HS":
                    m_rep_cnt['L1HS_all']+=1
                    fout_fa.write(">{0}{1}\n".format(read_id, s_pair))
                    fout_fa.write(seq+"\n")

                mapq = alignment.mapping_quality
                if mapq < 255:  # unique mapped for star
                    continue

                if sub_family not in m_rep_cnt:
                    m_rep_cnt[sub_family]=0
                m_rep_cnt[sub_family] += 1

            fout_unq.write("all\t{0}\n".format(n_all_reads)) #save all the reads
            for sub_family in m_rep_cnt:
                fout_unq.write("{0}\t{1}\n".format(sub_family, m_rep_cnt[sub_family]))

        bamfile.close()

#
# def XExpression():
#     def __init__(self, s_working_folder, n_jobs):
#         self.working_folder = s_working_folder
#         if self.working_folder[-1] != "/":
#             self.working_folder += "/"
#         self.n_jobs = n_jobs
#
#     def run_cmd(self, cmd):
#         Popen(cmd, shell=True, stdout=PIPE).communicate()
#


if __name__ == '__main__':
    sf_algnmt=sys.argv[1]
    sf_ref=sys.argv[2]
    sf_rmsk=sys.argv[3]
    sf_working_folder=sys.argv[4]
    n_jobs=int(sys.argv[5])
    rna_alnmt=RNAignment(sf_algnmt, sf_ref)
    rna_alnmt.cnt_unique_mapped_reads_within_rep(sf_rmsk, sf_working_folder, n_jobs)

