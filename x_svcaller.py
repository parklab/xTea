import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool

# pass out clipped reads
# 1. check how many reads are clipped (# of clipped by position),
# 2. generate a new bam file, with the "bad" clipped reads trimmed out

##Function: pool.map doesnt' accept function in the class
##So, use this wrapper to solve this

BWA_T_CUTOFF = 32  ##the minimum clipped length
NEARBY_REGION = 5
CLIP_FREQ = 10


##Function: pool.map doesnt' accept function in the class
##So, use this wrapper to solve this
def unwrap_self_cnt_clip(arg, **kwarg):
    return ClipReads.cnt_clip_by_chrm(*arg, **kwarg)

class ClipReads():
    def __init__(self, sf_bam, n_jobs):
        self.sf_bam = sf_bam
        self.working_folder = "./"
        self.n_jobs = n_jobs

    def set_working_folder(self, working_folder):
        self.working_folder = working_folder
        if working_folder[-1] != "/":
            self.working_folder += "/"

    def cnt_clip_by_chrm(self, record):
        chrm = record[0]
        sf_bam = record[1]

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

            map_pos = algnmt.reference_start
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                if l_cigar[0][0] == 4:  # soft-clip
                    if b_first:
                        cnt_soft_clipped_read1 += 1
                        if l_cigar[0][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read1 += 1
                    else:
                        cnt_soft_clipped_read2 += 1
                        if l_cigar[0][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read2 += 1
                elif l_cigar[0][0] == 5:  # hard-clip
                    if b_first:
                        cnt_hard_clipped_read1 += 1
                    else:
                        cnt_hard_clipped_read2 += 1
                if map_pos not in m_clip_pos:
                    m_clip_pos[map_pos] = 1
                else:
                    m_clip_pos[map_pos] += 1

            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                if l_cigar[-1][0] == 4:  # soft-clip
                    if b_first:
                        cnt_soft_clipped_read1 += 1
                        if l_cigar[-1][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read1 += 1
                    else:
                        cnt_soft_clipped_read2 += 1
                        if l_cigar[-1][1] >= BWA_T_CUTOFF:
                            cnt_soft_clipped_len30_read2 += 1
                elif l_cigar[-1][0] == 5:  # hard-clip
                    if b_first:
                        cnt_hard_clipped_read1 += 1
                    else:
                        cnt_hard_clipped_read2 += 1
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

        # write to files
        sf_statistic = self.working_folder + chrm + ".statistic"
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

        sf_clip_pos = self.working_folder + chrm + ".clip_pos"
        with open(sf_clip_pos, "w") as fout_clip_pos:
            for pos in m_clip_pos:
                fout_clip_pos.write(str(pos) + "\t" + str(m_clip_pos[pos])+"\n")
        samfile.close()


    def cnt_clipped_reads(self, sf_clip_statistic, sf_low_freq_clip_pos):
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
            sf_statistic = self.working_folder + chrm + ".statistic"
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

        m_low_freq_pos = {}
        m_clip_freq = {}
        for i in range(1, CLIP_FREQ + 1):
            m_clip_freq[i] = 0
        for chrm in references:
            if chrm not in m_low_freq_pos:
                m_low_freq_pos[chrm] = {}
            sf_clip_pos = self.working_folder + chrm + ".clip_pos"
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
            os.remove(sf_clip_pos)
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
                    m_low_freq_pos[chrm][pos] = freq

        with open(sf_low_freq_clip_pos, "w") as fout_low_freq_clip_pos:
            for chrm in m_low_freq_pos:
                for pos in m_low_freq_pos[chrm]:
                    freq=m_low_freq_pos[chrm][pos]
                    fout_low_freq_clip_pos.write(chrm+"\t"+str(pos)+"\t"+str(freq)+"\n")


        with open(sf_clip_statistic, "a") as fout_statistic:
            for freq in m_clip_freq:
                num_of_reads=m_clip_freq[freq]
                fout_statistic.write(str(freq)+"\t"+str(num_of_reads)+"\n")

if __name__ == '__main__':
    sf_bam=sys.argv[1]
    sf_clip_statistic=sys.argv[2]
    sf_low_freq_clip_pos=sys.argv[3]
    n_jobs=int(sys.argv[4])

    cr=ClipReads(sf_bam, n_jobs)
    cr.cnt_clipped_reads(sf_clip_statistic, sf_low_freq_clip_pos)