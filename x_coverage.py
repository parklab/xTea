import os
import pysam
from x_alignments import *
from multiprocessing import Pool

def unwrap_self_calc_depth_for_site(arg, **kwarg):
    return ReadDepth._calc_depth_one_site(*arg, **kwarg)

class ReadDepth():
    def __init__(self, working_folder, n_jobs, sf_ref):
        self.working_folder = working_folder
        self.n_jobs = n_jobs
        self.sf_reference = sf_ref

    def _get_map_end_pos(self, l_cigar, start_pos):
        end_pos=start_pos
        for (type, lenth) in l_cigar:
            if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                continue
            else:
                end_pos += lenth
        return end_pos

    def __calc_coverage(self, mcov, win_size):
        isum=0
        for tmp_pos in mcov:
            isum+=mcov[tmp_pos]
        fcov=float(isum)/float(win_size)
        return fcov

    #calc the read depth for one site
    def _calc_depth_one_site(self, record):
        chrm=record[0][0]
        insertion_pos=record[0][1]
        search_win=record[0][2]
        focal_win=record[0][3]

        mlcov={}#save the coverage of left bases
        mrcov={}#save the coverage of right bases

        start_pos = insertion_pos - search_win
        if start_pos <= 0:
            start_pos = 1
        end_pos = insertion_pos + search_win

        sf_bam=record[1]
        working_folder = record[2]

        bam_info = BamInfo(sf_bam, self.sf_reference)
        b_with_chr = bam_info.is_chrm_contain_chr()
        chrm_in_bam = self._process_chrm_name(b_with_chr, chrm)
        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        for algnmt in samfile.fetch(chrm_in_bam, start_pos, end_pos):  ##fetch reads mapped to "chrm:start_pos-end_pos"
            ##here need to skip the secondary and supplementary alignments
            if algnmt.is_secondary or algnmt.is_supplementary:
                continue
            if algnmt.is_duplicate == True:  ##duplciate
                continue
            if algnmt.is_unmapped == True:  ##unmapped
                continue
            map_pos = algnmt.reference_start #mapping position
            l_cigar = algnmt.cigar
            tmp_end_pos=self._get_map_end_pos(l_cigar, map_pos)
            lfocal_start=insertion_pos-focal_win
            lfocal_end=insertion_pos
            rfocal_start=insertion_pos+1
            rfocal_end=insertion_pos+focal_win
            #print map_pos, tmp_end_pos, insertion_pos
            for i in range(map_pos, tmp_end_pos):
                if i>=lfocal_start and i<=lfocal_end:
                    if i not in mlcov:
                        mlcov[i]=1
                    else:
                        mlcov[i]+=1
                if i>=rfocal_start and i<=rfocal_end:
                    if i not in mrcov:
                        mrcov[i]=1
                    else:
                        mrcov[i]+=1
        samfile.close()

        #calc the coverage
        flcov=self.__calc_coverage(mlcov, focal_win)
        frcov=self.__calc_coverage(mrcov, focal_win)

        return (chrm, insertion_pos, flcov, frcov)

####

    #calculate the local read depth of given sites
    #sf_sites: the interested sites
    #search_win: collect reads in this range
    #focal_win: calc the depth in this range
    def calc_coverage_of_sites(self, sf_sites, sf_bam_list, search_win, focal_win):
        m_site_cov = {}
        with open(sf_bam_list) as fin_bams:
            for line in fin_bams:
                sf_bam=line.rstrip()

                l_records = []
                with open(sf_sites) as fin_list:
                    for line in fin_list:
                        fields = line.split()
                        chrm = fields[0]
                        pos = int(fields[1])  # candidate insertion site
                        l_records.append(((chrm, pos, search_win, focal_win), sf_bam, self.working_folder))

                pool = Pool(self.n_jobs)
                l_sites_cov=pool.map(unwrap_self_calc_depth_for_site, zip([self] * len(l_records), l_records), 1)
                pool.close()
                pool.join()

                for record in l_sites_cov:
                    ins_chrm=record[0]
                    ins_pos=record[1]
                    flcov=record[2]
                    frcov=record[3]

                    if ins_chrm not in m_site_cov:
                        m_site_cov[ins_chrm]={}
                    if ins_pos not in m_site_cov[ins_chrm]:
                        m_site_cov[ins_chrm][ins_pos] = []
                        m_site_cov[ins_chrm][ins_pos].append(0)
                        m_site_cov[ins_chrm][ins_pos].append(0)

                    m_site_cov[ins_chrm][ins_pos][0] += flcov
                    m_site_cov[ins_chrm][ins_pos][1] += frcov
        return m_site_cov

    ####
    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, b_tmplt_with_chr, chrm):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True
        # print chrm, self.b_with_chr, b_chrm_with_chr #################################################################

        if b_tmplt_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_tmplt_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_tmplt_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

####for test only
# if __name__ == '__main__':
#     sf_sites=sys.argv[1]
#     sf_bam_list=sys.argv[2]
#     sf_ref=sys.argv[3]
#     n_jobs=int(sys.argv[4])
#     wfolder=sys.argv[5]
#
#     rd=ReadDepth(wfolder, n_jobs, sf_ref)
#     search_win=200
#     focal_win=100
#     m_cov=rd.calc_coverage_of_sites(sf_sites, sf_bam_list, search_win, focal_win)
#     print m_cov
