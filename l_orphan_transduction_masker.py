import global_values
import pysam

class OrphanTransduction():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        self.n_jobs = n_jobs


    def slct_unqiue_algnmt_from_ref_algnmt(self, sf_algnmt):
        s_open_fmt = "rb"
        samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
        for algnmt in samfile.fetch():  # check each alignment
            mapq = algnmt.mapping_quality
            if mapq < global_values.TRANSDCT_UNIQ_MAPQ:
                continue
        samfile.close()
####