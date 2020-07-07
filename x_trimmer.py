##11/22/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
from x_converter import BarcodeHeader

# pass out clipped reads
# 1. check how many reads are clipped (# of clipped by position),
# 2. generate a new bam file, with the "bad" clipped reads trimmed out

##Function: pool.map doesnt' accept function in the class
##So, use this wrapper to solve this

BWA_T_CUTOFF = 32  ##the minimum clipped length
BWA_REALIGN_CUTOFF = 11
NEARBY_REGION = 5
CLIP_FREQ = 10
TRIM_CLIP_FREQ = 2

BWA_PATH = "bwa"
SAMTOOLS_PATH = "samtools"
STATISTIC_SUFFIX = ".statistic"
CLIP_FQ_SUFFIX = ".clipped.fq"
CLIP_BAM_SUFFIX = ".clipped.sam"
CLIP_POS_SUFFIX = ".clip_pos"
OUTPUT_BAM_SUFFIX = ".out_bam"
OUTPUT_BAM_HEADER = ".bam_header.sam"


# SAMTOOLS_PATH="samtools"

##Function: pool.map doesnt' accept function in the class
##So, use this wrapper to solve this
def unwrap_self_cnt_clip(arg, **kwarg):
    return ReadTrimmer.collect_clip_by_chrm(*arg, **kwarg)


def unwrap_self_parse_clip(arg, **kwarg):
    return ReadTrimmer.filter_out_low_freq_clip_by_chrm(*arg, **kwarg)


class ReadTrimmer():
    def __init__(self, sf_bam, n_jobs):
        self.sf_bam = sf_bam
        self.working_folder = "./"
        self.n_jobs = n_jobs

    def set_working_folder(self, working_folder):
        self.working_folder = working_folder
        if working_folder[-1] != "/":
            self.working_folder += "/"

    def collect_clip_by_chrm(self, record):
        chrm = record[0]
        sf_bam = record[1]

        sf_clip_fq = self.working_folder + chrm + CLIP_FQ_SUFFIX  # this is to save the clipped part for re-alignment
        f_clip_fq = open(sf_clip_fq, "w")

        cnt_unmapped_read1 = 0
        cnt_unmapped_read2 = 0
        cnt_soft_clipped_read1 = 0
        cnt_soft_clipped_read2 = 0
        cnt_hard_clipped_read1 = 0
        cnt_hard_clipped_read2 = 0
        cnt_soft_clipped_len30_read1 = 0  ##the mapped part should at least 30M
        cnt_soft_clipped_len30_read2 = 0
        cnt_duplicate = 0
        cnt_not_barcoded = 0
        cnt_all_reads = 0
        samfile = pysam.AlignmentFile(sf_bam, "rb")
        m_clip_pos = {}
        for algnmt in samfile.fetch(chrm):  ##fetch reads mapped to "chrm"
            ##here need to skip the secondary and supplementary alignments?
            # if algnmt.is_secondary or algnmt.is_supplementary:
            #     continue
            cnt_all_reads += 1
            if algnmt.is_duplicate == True:  ##duplciate
                cnt_duplicate += 1
                continue

            b_first = True
            if algnmt.is_read2 == True:
                b_first = False
            if algnmt.is_unmapped == True:  # unmapped
                if b_first:
                    cnt_unmapped_read1 += 1
                else:
                    cnt_unmapped_read2 += 1
                continue

            if algnmt.has_tag("BX") == False:
                cnt_not_barcoded += 1

            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue
            if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
                continue

            query_name = algnmt.query_name
            query_seq = algnmt.query_sequence
            query_quality = algnmt.query_qualities  ##this is different from the one saved in the fastq/sam, no offset 33 to subtract
            map_pos = algnmt.reference_start
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                if map_pos not in m_clip_pos:
                    m_clip_pos[map_pos] = 1
                else:
                    m_clip_pos[map_pos] += 1

                if algnmt.is_supplementary or algnmt.is_secondary:  ###secondary and supplementary are not considered
                    continue
                if l_cigar[0][0] == 4:  # soft-clip
                    clipped_seq = query_seq[:l_cigar[0][1]]
                    clipped_qulity = self._cvt_to_Ascii_quality(query_quality[:l_cigar[0][1]])
                    clipped_rname = query_name + "_2"
                    if b_first:
                        cnt_soft_clipped_read1 += 1
                        if l_cigar[0][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read1 += 1
                        clipped_rname = query_name + "_1"
                    else:
                        cnt_soft_clipped_read2 += 1
                        if l_cigar[0][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read2 += 1
                    f_clip_fq.write(clipped_rname + "\n")
                    f_clip_fq.write(clipped_seq + "\n+\n")
                    f_clip_fq.write(clipped_qulity + "\n")
                elif l_cigar[0][0] == 5:  # hard-clip
                    if b_first:
                        cnt_hard_clipped_read1 += 1
                    else:
                        cnt_hard_clipped_read2 += 1

            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                ##calculate the exact clip position
                for (type, lenth) in l_cigar[:-1]:
                    if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                        continue
                    else:
                        map_pos += lenth

                if map_pos not in m_clip_pos:
                    m_clip_pos[map_pos] = 1
                else:
                    m_clip_pos[map_pos] += 1

                if algnmt.is_supplementary or algnmt.is_secondary:  ###secondary and supplementary are not considered
                    continue
                if l_cigar[-1][0] == 4:  #soft-clip
                    clipped_rname = query_name + "_2"
                    start_pos = -1 * l_cigar[-1][1]
                    clipped_seq = query_seq[start_pos:]
                    clipped_qulity = self._cvt_to_Ascii_quality(query_quality[start_pos:])
                    if b_first:
                        cnt_soft_clipped_read1 += 1
                        if l_cigar[-1][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read1 += 1
                        clipped_rname = query_name + "_1"
                    else:
                        cnt_soft_clipped_read2 += 1
                        if l_cigar[-1][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read2 += 1
                    f_clip_fq.write(clipped_rname + "\n")
                    f_clip_fq.write(clipped_seq + "\n+\n")
                    f_clip_fq.write(clipped_qulity + "\n")
                elif l_cigar[-1][0] == 5:  # hard-clip
                    if b_first:
                        cnt_hard_clipped_read1 += 1
                    else:
                        cnt_hard_clipped_read2 += 1
        f_clip_fq.close()

        # write to files
        sf_statistic = self.working_folder + chrm + STATISTIC_SUFFIX
        with open(sf_statistic, "w") as fout_statistic:
            fout_statistic.write(str(cnt_unmapped_read1) + "\t")
            fout_statistic.write(str(cnt_unmapped_read2) + "\t")
            fout_statistic.write(str(cnt_soft_clipped_read1) + "\t")
            fout_statistic.write(str(cnt_soft_clipped_read2) + "\t")
            fout_statistic.write(str(cnt_hard_clipped_read1) + "\t")
            fout_statistic.write(str(cnt_hard_clipped_read2) + "\t")
            fout_statistic.write(str(cnt_soft_clipped_len30_read1) + "\t")
            fout_statistic.write(str(cnt_soft_clipped_len30_read2) + "\t")
            fout_statistic.write(str(cnt_duplicate) + "\t")
            fout_statistic.write(str(cnt_not_barcoded) + "\t")
            fout_statistic.write(str(cnt_all_reads) + "\n")

        sf_clip_pos = self.working_folder + chrm + CLIP_POS_SUFFIX
        with open(sf_clip_pos, "w") as fout_clip_pos:
            for pos in m_clip_pos:
                fout_clip_pos.write(str(pos) + "\t" + str(m_clip_pos[pos]) + "\n")
        samfile.close()

    def _cvt_to_Ascii_quality(self, l_score):
        new_score = [x + 33 for x in l_score]
        return ''.join(map(chr, new_score))

    def collect_clipped_reads(self, sf_clip_statistic, sf_all_clip_fq):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb")
        references = samfile.references
        l_chrm_records = []
        for chrm in references:
            l_chrm_records.append((chrm, self.sf_bam, self.working_folder))
        samfile.close()

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_cnt_clip, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        cnt_unmapped_read1 = 0
        cnt_unmapped_read2 = 0
        cnt_soft_clipped_read1 = 0
        cnt_soft_clipped_read2 = 0
        cnt_hard_clipped_read1 = 0
        cnt_hard_clipped_read2 = 0
        cnt_soft_clipped_len30_read1 = 0  ##the mapped part should at least 30M
        cnt_soft_clipped_len30_read2 = 0
        cnt_duplicate = 0
        cnt_not_barcoded = 0
        cnt_all_reads = 0
        # merge the statistic
        for chrm in references:
            sf_statistic = self.working_folder + chrm + STATISTIC_SUFFIX
            if os.path.isfile(sf_statistic) == False:
                print "Error: Statistic file for chrom {0} doesn't exist!!!!".format(chrm)
                continue
            with open(sf_statistic) as fin_statistic:
                for line in fin_statistic:
                    fields = line.split()
                    cnt_unmapped_read1 += int(fields[0])
                    cnt_unmapped_read2 += int(fields[1])
                    cnt_soft_clipped_read1 += int(fields[2])
                    cnt_soft_clipped_read2 += int(fields[3])
                    cnt_hard_clipped_read1 += int(fields[4])
                    cnt_hard_clipped_read2 += int(fields[5])
                    cnt_soft_clipped_len30_read1 += int(fields[6])  ##the mapped part should at least 30M
                    cnt_soft_clipped_len30_read2 += int(fields[7])
                    cnt_duplicate += int(fields[8])
                    cnt_not_barcoded += int(fields[9])
                    cnt_all_reads += int(fields[10])
            os.remove(sf_statistic)

        with open(sf_clip_statistic, "w") as fout_statistic:
            fout_statistic.write(str(cnt_unmapped_read1) + "\t" + "unmapped_read1\n")
            fout_statistic.write(str(cnt_unmapped_read2) + "\t" + "unmapped_read2\n")
            fout_statistic.write(str(cnt_soft_clipped_read1) + "\t" + "soft_clipped_read1\n")
            fout_statistic.write(str(cnt_soft_clipped_read2) + "\t" + "soft_clipped_read2\n")
            fout_statistic.write(str(cnt_hard_clipped_read1) + "\t" + "hard_clipped_read1\n")
            fout_statistic.write(str(cnt_hard_clipped_read2) + "\t" + "hard_clipped_read2\n")
            fout_statistic.write(str(cnt_soft_clipped_len30_read1) + "\t" + "soft_clipped_len30_read1\n")
            fout_statistic.write(str(cnt_soft_clipped_len30_read2) + "\t" + "soft_clipped_len30_read2\n")
            fout_statistic.write(str(cnt_duplicate) + "\t" + "duplicate\n")
            fout_statistic.write(str(cnt_not_barcoded) + "\t" + "non_barcoded\n")
            fout_statistic.write(str(cnt_all_reads) + "\tall_reads\n")

        # m_low_freq_pos = {}
        m_clip_freq = {}
        for i in range(1, CLIP_FREQ + 1):
            m_clip_freq[i] = 0
        for chrm in references:
            # if chrm not in m_low_freq_pos:
            #    m_low_freq_pos[chrm] = {}
            sf_clip_pos = self.working_folder + chrm + CLIP_POS_SUFFIX
            if os.path.isfile(sf_clip_pos) == False:
                print "Error: Position file for chrom {0} doesn't exist!!!!".format(chrm)
                continue

            m_clip_pos_freq = {}
            with open(sf_clip_pos) as fin_clip_pos:
                for line in fin_clip_pos:
                    fields = line.split()
                    cur_pos = int(fields[0])
                    cur_cnt = int(fields[1])
                    m_clip_pos_freq[cur_pos] = cur_cnt
            # os.remove(sf_clip_pos) ########################################################################################
            ##check nearby region,whether there is a clip
            for pos in m_clip_pos_freq:
                freq = m_clip_pos_freq[pos]
                if freq > CLIP_FREQ:
                    continue
                nearby_freq = 0
                for i in range(-1 * NEARBY_REGION, NEARBY_REGION):
                    if (pos + i) in m_clip_pos_freq:
                        nearby_freq += m_clip_pos_freq[pos + i]
                if nearby_freq < CLIP_FREQ:
                    m_clip_freq[freq] += 1
                    # m_low_freq_pos[chrm][pos] = freq

        # with open(sf_low_freq_clip_pos, "w") as fout_low_freq_clip_pos:
        #     for chrm in m_low_freq_pos:
        #         for pos in m_low_freq_pos[chrm]:
        #             freq = m_low_freq_pos[chrm][pos]
        #             fout_low_freq_clip_pos.write(chrm + "\t" + str(pos) + "\t" + str(freq) + "\n")

        with open(sf_clip_statistic, "a") as fout_statistic:
            for freq in m_clip_freq:
                num_of_reads = m_clip_freq[freq]
                fout_statistic.write(str(freq) + "\t" + str(num_of_reads) + "\n")

        ##merge the clipped reads
        # sf_all_clip_fq = self.working_folder + ".clipped.fq"
        with open(sf_all_clip_fq, "w") as fout_all:
            for chrm in references:
                sf_clip_fq = self.working_folder + chrm + CLIP_FQ_SUFFIX
                with open(sf_clip_fq) as fin_clip:
                    for line in fin_clip:
                        fout_all.write(line)

    ##re-align the soft-clipped reads
    ##this is limited to primary alignments and non supplementary alignments
    def realign_clipped_reads_to_reference(self, sf_ref, sf_all_clip_fq, sf_clip_bam):
        cmd = "{0} mem -T {1} -k {2} -t {3} {4} {5} > {6}".format(BWA_PATH, BWA_REALIGN_CUTOFF, BWA_REALIGN_CUTOFF,
                                                                  self.n_jobs, sf_ref, sf_all_clip_fq, sf_clip_bam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    ####This only load the mapped reads, and keep in a set
    def _load_clipped_algnmt(self, sf_clip_algnmt):
        set_mapped_clip_part = set()
        sam_file = pysam.AlignmentFile(sf_clip_algnmt, "r")
        for algnmt in sam_file.fetch():
            if algnmt.is_unmapped == True:
                continue
            query_name = algnmt.query_name
            set_mapped_clip_part.add(query_name)
        sam_file.close()
        return set_mapped_clip_part

    def _load_low_freq_clip_pos(self, sf_low_freq_clip_pos):
        m_low_freq_clip_pos = {}
        with open(sf_low_freq_clip_pos) as fin_low_freq_clip_pos:
            for line in fin_low_freq_clip_pos:
                fields = line.split()
                pos = int(fields[0])
                freq = int(fields[1])
                m_low_freq_clip_pos[pos] = freq
        return m_low_freq_clip_pos

    ##return true if the position has less than TRIM_CLIP_FREQ reads clipped there
    def _is_low_freq_clip_pos(self, pos, m_clip_pos_freq):
        ##check nearby region,whether there is a clip
        freq = m_clip_pos_freq[pos]
        if freq > TRIM_CLIP_FREQ:
            return False
        nearby_freq = 0
        for i in range(-1 * NEARBY_REGION, NEARBY_REGION):
            if (pos + i) in m_clip_pos_freq:
                nearby_freq += m_clip_pos_freq[pos + i]
        if nearby_freq > TRIM_CLIP_FREQ:
            return False
        return True

    ##trim the clipped alignment and keep the region within [start, end]
    ##note: seq, quality and cigar will be changed
    def _trim_clipped_alignment(self, istart, iend, cigar_tuples, alignmt):
        clipped_seq = alignmt.query_sequence[istart:iend]
        clipped_qulity = alignmt.query_qualities[istart:iend]
        if iend == 0:
            clipped_seq = alignmt.query_sequence[istart:]
            clipped_qulity = alignmt.query_qualities[istart:]

        a = pysam.AlignedSegment()  # new pysam.AlignedSegment()
        a.query_name = alignmt.query_name
        a.query_sequence = clipped_seq
        a.flag = alignmt.flag
        a.reference_id = alignmt.reference_id
        a.reference_start = alignmt.reference_start
        a.mapping_quality = alignmt.mapping_quality
        a.cigartuples = cigar_tuples
        a.next_reference_id = alignmt.next_reference_id
        a.next_reference_start = alignmt.next_reference_start
        a.template_length = alignmt.template_length
        a.query_qualities = clipped_qulity
        a.set_tags(alignmt.get_tags())
        return a

    ####
    def filter_out_low_freq_clip_by_chrm(self, record):
        chrm = record[0]
        sf_bam = record[1]
        sf_clip_alignmt = record[2]
        working_folder = record[3]
        sf_clip_pos = working_folder + chrm + CLIP_POS_SUFFIX
        sf_out_bam_chrm = working_folder + chrm + OUTPUT_BAM_SUFFIX

        print "Load clipped alignments for ", chrm
        set_mapped_clip_part = self._load_clipped_algnmt(sf_clip_alignmt)
        print "Load clip positions for ", chrm
        m_low_freq_clip_pos = self._load_low_freq_clip_pos(sf_clip_pos)

        samfile = pysam.AlignmentFile(sf_bam, "rb")
        merged_bam = pysam.AlignmentFile(sf_out_bam_chrm, "wb", template=samfile)
        #print "Working on ", chrm, " now ..."
        for algnmt in samfile.fetch(chrm):
            # first check whether read is clipped or not
            if algnmt.is_duplicate == True:  ##remove duplicate reads
                continue
            if algnmt.is_unmapped == True:  # unmapped
                merged_bam.write(algnmt)
                continue
            if algnmt.is_supplementary or algnmt.is_secondary:  ###secondary and supplementary are not considered
                merged_bam.write(algnmt)
                continue

            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue
            if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
                merged_bam.write(algnmt)
                continue

            if l_cigar[0][0] == 5 or l_cigar[-1][0] == 5:  # hard clip
                merged_bam.write(algnmt)
                continue

            if len(l_cigar) > 1 and l_cigar[0][0] != 4 and l_cigar[-1][0] != 4: #both end not clipped
                merged_bam.write(algnmt)
                continue
            query_name = algnmt.query_name
            map_pos = int(algnmt.reference_start)
            cigar_tuples = algnmt.cigartuples
            clipped_rname = query_name + "_1"
            b_first = True  ##first in pair by default
            if algnmt.is_read2 == True:  ##second in pair
                b_first = False
                clipped_rname = query_name + "_2"
            istart = 0
            iend = 0
            # if l_cigar[0][0] == 4 and l_cigar[-1][0] == 4:#both end are soft-clipped
            if l_cigar[0][0] == 4:  # left soft-clipped
                ##first, check the frequency at the clip position
                istart = int(l_cigar[0][1])
                del cigar_tuples[0]

            if l_cigar[-1][0] == 4:  # right soft-clipped
                ##calculate the exact clip position
                for (type, lenth) in l_cigar[:-1]:
                    if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                        continue
                    else:
                        map_pos += int(lenth)
                iend = (-1 * l_cigar[-1][1])
                del cigar_tuples[-1]
            if map_pos not in m_low_freq_clip_pos:
                print "Error happen: {0}:{1} not in library!!!!!".format(chrm, map_pos)################################
                print cigar_tuples, clipped_rname
                break ##################################################################################################
            if self._is_low_freq_clip_pos(map_pos, m_low_freq_clip_pos) == False:
                merged_bam.write(algnmt)
                continue
            # then, check whether the clipped part is mappable, if yes, then keep it in the alignment
            ##otherwise, trim the original read
            if clipped_rname in set_mapped_clip_part:
                merged_bam.write(algnmt)
                continue
            # trim the alignment
            new_alignmt = self._trim_clipped_alignment(istart, iend, cigar_tuples, algnmt)
            merged_bam.write(new_alignmt)

        merged_bam.close()
        samfile.close()

    ####
    ##Input: 1)Original bam, 2)clip_pos_freq file, 3)re-alignment of the clipped regions
    ##Output: processed bam file
    def filter_out_low_freq_clip_reads(self, sf_clip_algnmt, sf_out_bam):
        samfile = pysam.AlignmentFile(self.sf_bam, "rb")
        references = samfile.references
        l_chrm_records = []
        for chrm in references:
            l_chrm_records.append((chrm, self.sf_bam, sf_clip_algnmt, self.working_folder))
        samfile.close()

        # print "Begin to map ..."
        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_parse_clip, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()
        # print "Begin to reduce ..."
        #self.filter_out_low_freq_clip_by_chrm(l_chrm_records[0]) ######################################################

        bh = BarcodeHeader(self.sf_bam, self.working_folder, self.n_jobs)
        sf_old_bam_header = self.working_folder + OUTPUT_BAM_HEADER
        bh.write_old_header_to_file(sf_old_bam_header)
        ##merge the alignments
        cmd = "{0} cat -h {1} -o {2} ".format(SAMTOOLS_PATH, sf_old_bam_header, sf_out_bam)
        for chrm in references:
            sf_out_bam_chrm = self.working_folder + chrm + OUTPUT_BAM_SUFFIX
            if os.path.isfile(sf_out_bam_chrm):
                cmd += "{0} ".format(sf_out_bam_chrm)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    ####
    def trim_reads(self, s_working_folder, sf_ref, sf_trimmed_bam):
        self.set_working_folder(s_working_folder)
        sf_clip_algnmt = self.working_folder + CLIP_BAM_SUFFIX
        sf_all_clip_fq = self.working_folder + CLIP_FQ_SUFFIX
        sf_clip_statistic = self.working_folder + STATISTIC_SUFFIX

        self.collect_clipped_reads(sf_clip_statistic, sf_all_clip_fq)
        self.realign_clipped_reads_to_reference(sf_ref, sf_all_clip_fq, sf_clip_algnmt)
        self.filter_out_low_freq_clip_reads(sf_clip_algnmt, sf_trimmed_bam)

####
###TDList:
###1.for bam sorted by name: check the clipped, also has a supplementary alignment alignment
###2. get all the clipped reads (filter out those sites only have one read clip)
###3. re-align the parsed reads
####
