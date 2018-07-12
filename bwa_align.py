from subprocess import *

class BWAlign():
    def __init__(self, BWA_PATH, BWA_REALIGN_CUTOFF, n_jobs):
        self.BWA_PATH=BWA_PATH
        self.BWA_REALIGN_CUTOFF=BWA_REALIGN_CUTOFF
        self.n_jobs = n_jobs

    # re-align the collected clipped and discordant reads
    def realign_clipped_reads(self, sf_ref, sf_reads, sf_out_sam):
        cmd = "{0} mem -t {1} -T {2} -k {3} {4} {5} > {6}".format(self.BWA_PATH, self.n_jobs, self.BWA_REALIGN_CUTOFF,
                                                                  self.BWA_REALIGN_CUTOFF, sf_ref, sf_reads, sf_out_sam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    # re-align the collected clipped and discordant reads
    def realign_disc_reads(self, sf_ref, sf_reads, sf_out_sam):
        cmd = "{0} mem -t {1} {2} {3} > {4}".format(self.BWA_PATH, self.n_jobs, sf_ref, sf_reads, sf_out_sam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    # re-align the collected reads
    def realign_reads(self, sf_ref, sf_reads, sf_out_sam):
        cmd = "{0} mem -t {1} {2} {3} > {4}".format(self.BWA_PATH, self.n_jobs, sf_ref, sf_reads, sf_out_sam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    # re-align the collected reads
    def realign_reads_to_bam(self, SAMTOOLS, sf_ref, sf_reads, sf_out_bam):
        cmd = "{0} mem -t {1} {2} {3} | {4} view -hSb - | {5} sort -o {6} - && {7} index {8}".format(
            self.BWA_PATH, self.n_jobs, sf_ref, sf_reads, SAMTOOLS, SAMTOOLS, sf_out_bam, SAMTOOLS, sf_out_bam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    ####