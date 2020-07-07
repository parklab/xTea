import pysam
from bwa_align import *
from l_rep_masker import *
from l_contig import *
from x_polyA import *

#To-do-list: check copied from which transcriptome, which means find the exon combination. Not only report the gene

class PsudogeneMasker(LRepMasker):
    def __init__(self, swfolder, n_jobs):
        LRepMasker.__init__(self, swfolder, n_jobs)
        self._psudogene = "PSUDOGENE"
        self._polyA_ck_win=50

    ####
    def align_exon_2_contigs(self, sf_contigs, sf_exons, sf_out_bam):
        bwa_algn = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, self.n_jobs)
        bwa_algn.index_ref_file(sf_contigs)#index the contig
        bwa_algn.align_exon_to_contig(sf_contigs, sf_exons, sf_out_bam)

    #parse the algnmt (by aligning exons to contigs)
    #bam is sorted
    def parse_exon_2_contig_algnmt(self, sf_bam, sf_contig, f_max_clip_cutoff, f_min_cov, sf_out):
        m_candidates={}
        bam_info = BamInfo(sf_bam, sf_contig)
        m_contigs_cov={}
        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=sf_contig)
        for algnmt in samfile.fetch():
            if algnmt.is_secondary or algnmt.is_supplementary == True:#skip supplementary, but keep the secondary alignment
                continue
            if algnmt.is_duplicate == True:  ##duplciate
                continue
            if algnmt.is_unmapped == True:  ##unmapped
                continue
            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue
            if algnmt.mapping_quality < global_values.MINIMUM_DISC_MAPQ:#by default this is set to 20
                continue

            i_exon_len = len(algnmt.query_sequence)
            i_start = algnmt.reference_start##
            i_end=i_start+i_exon_len#
            i_clip_len=0
            #first make sure the exon is fully (or almost) mapped
            if l_cigar[0][0] == 4 or l_cigar[0][0] ==5:  # left clipped, right breakpoints
                i_clip_len+=(l_cigar[0][1])
            if l_cigar[-1][0] == 4 or l_cigar[-1][0] ==5:  # right clipped, left breakpoints
                i_clip_len += (l_cigar[-1][1])
                i_end-=i_clip_len

            if i_exon_len<=0:
                continue
            if float(i_clip_len)/float(i_exon_len)>f_max_clip_cutoff:
                continue

            s_contig_id = algnmt.reference_name
            if s_contig_id not in m_contigs_cov:
                m_contigs_cov[s_contig_id]={}
            # example id: chr20:57495966:57496364:exon:ENST00000423479.7:12:CTCF
            query_name = algnmt.query_name
            tmp_fields = query_name.split(":")
            gene_name=tmp_fields[-1]
            if gene_name not in m_contigs_cov[s_contig_id]:
                m_contigs_cov[s_contig_id][gene_name]=[]
            m_contigs_cov[s_contig_id][gene_name].append((i_start, i_end))

        #second, make sure the contig is fully (or almost) covered by the same gene
        m_contig_lenth = bam_info.get_all_chrom_name_length()
        for s_contig in m_contigs_cov:
            if s_contig not in m_contig_lenth:
                print("% is not found in alignment header!!! Abnormal!" % s_contig)
                continue
            i_contig_len=m_contig_lenth[s_contig]
            for gene_name in m_contigs_cov[s_contig]:
                l_cov=m_contigs_cov[s_contig][gene_name]
                i_tmp_cov=self.calc_cov_length(l_cov)
                if float(i_tmp_cov)/float(i_contig_len)>f_min_cov:
                    m_candidates[s_contig]=(gene_name, l_cov)
                    break
####
        #third, request polyA/T
        lrd_ctg=LRD_Contig(self.swfolder, self.n_jobs)
        m_contig_seqs=lrd_ctg.load_fasta_to_dict(sf_contig)

        #is_consecutive_polyA_T(self, seq)
        m_new_candidates=self.check_polyA_T(m_candidates, m_contig_seqs, self._polyA_ck_win)#
        samfile.close()
        self.dump_rslts(m_new_candidates, m_contig_seqs, sf_out)

    def check_polyA_T(self, m_candidates, m_contigs, n_chk_size):
        m_new_candidates={}
        plyA = PolyA()
        for s_contig in m_candidates:
            if s_contig not in m_contigs:
                print("% is not found in m_contigs!!! Abnormal!" % s_contig)
                continue
            s_seq=m_contigs[s_contig]
            if plyA.is_consecutive_polyA_T(s_seq[-1*n_chk_size:])==True or \
                            plyA.is_consecutive_polyA_T(s_seq[:n_chk_size])==True:
                m_new_candidates[s_contig]=m_candidates[s_contig]
        return m_new_candidates

    #calculate the total covered length
    def calc_cov_length(self, l_cov):
        l_cov.sort(key = lambda x: x[0])
        i_beg=0
        i_tail=0
        i_total=0
        for (i_start, i_end) in l_cov:
            if i_start>=i_tail:
                i_total+=(i_end-i_start+1)
                i_beg=i_start
                i_tail=i_end
            elif i_end>i_tail:
                    i_total += (i_end - i_tail+1)
                    i_tail=i_end
        return i_total

    #dump the results
    def dump_rslts(self, m_candidates, m_contig_seqs, sf_out):
        with open(sf_out,"w") as fout:
            for s_contig in m_candidates:
                s_seq=m_contig_seqs[s_contig]#contig sequence
                (gene_name, l_cov)=m_candidates[s_contig]
                fields=s_contig.split(global_values.SEPERATOR)
                fout.write(fields[0]+"\t"+fields[1]+"\t.\t.\t"+s_seq+"\t60\t.\tSRC_GENE="+gene_name+";SEGMTS=")
                b_flag=True
                m_cov={}
                for (i_start, i_end) in l_cov:
                    s_id=str(i_start)+"-"+str(i_end)
                    if s_id in m_cov:
                        continue
                    if b_flag==True:
                        fout.write(s_id)
                        m_cov[s_id]=1
                        b_flag=False
                    else:
                        fout.write(","+s_id)
                        m_cov[s_id] = 1
                fout.write("\n")

####