##01/18/2022
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

#Given a list of sites, get those TE insertion related reads only, put them in a separate bam file for bamsnap/igv

import os
from x_alignments import *
from x_basic_info import *
from global_values import *
import pysam

class XScreenshotBAM():#####
    def __init__(self, sf_sites, sf_bam_list, sf_ref, n_jobs, sf_wfolder):
        self.sf_sites=sf_sites
        self.sf_bam_list=sf_bam_list
        self.sf_ref=sf_ref
        self.sf_wfolder=sf_wfolder
        if len(self.sf_wfolder)>0 and self.sf_wfolder[-1]!="/":
            self.sf_wfolder+="/"#
        if len(self.sf_wfolder)<=0:
            self.sf_wfolder = "./"
        self.n_jobs=n_jobs

    def chk_index_bam(self, sf_bam_list):
        b_exist=True
        s_sample_name, l_ori_bams = self._load_in_bams(sf_bam_list)
        for sf_bam in l_ori_bams:
            if os.path.isfile(sf_bam)==False:##
                b_exist=False
                break
            #index the bam
            if os.path.isfile(sf_bam+".bai")==False:
                self._index_algnmt(sf_bam)
        return b_exist

    def _index_algnmt(self, sf_bam):#
        cmd_runner = CMD_RUNNER()
        #sort the bam first, and then index it
        cmd1= "%s sort -o %s.sorted %s" % (global_values.SAMTOOLS_PATH, sf_bam, sf_bam)
        cmd_runner.run_cmd_small_output(cmd1)
        #copy the file to cover the ori
        cmd2="mv %s.sorted %s" % (sf_bam, sf_bam)
        cmd_runner.run_cmd_small_output(cmd2)
        cmd = "{0} index {1}".format(global_values.SAMTOOLS_PATH, sf_bam)
        cmd_runner.run_cmd_small_output(cmd)
#

    def _gnrt_new_bam_name(self, s_old_name):
        sf_screenshot_bam = self.sf_wfolder + "screenshot_" + s_old_name + ".bam"
        return sf_screenshot_bam

    ####
    ####
    def _load_in_bams(self, sf_bam_list):
        l_ori_bams = []
        s_sample_name = ""
        with open(sf_bam_list) as fin_bams:
            for line in fin_bams:
                fields=line.rstrip().split()
                if len(fields)<2:
                    return
                s_sample_name=fields[0]
                for s_tmp_bam in fields[1:]:
                    l_ori_bams.append(s_tmp_bam)
                break
        return s_sample_name, l_ori_bams

    def prepare_new_screenbam_name_list(self, sf_new_list):#
        with open(sf_new_list, "w") as fout_bam_list:  #
            s_sample_name, l_ori_bams = self._load_in_bams(self.sf_bam_list)
            s_new_list = s_sample_name
            for sf_bam in l_ori_bams:
                s_fname = os.path.basename(sf_bam)  #
                sf_screenshot_bam = self._gnrt_new_bam_name(s_fname)  #
                s_new_list += (" " + sf_screenshot_bam)
            fout_bam_list.write(s_new_list+"\n")


    ####
    def gnrt_screenshot_bam(self, i_clip_extnd):
        #first load in all the sites
        l_sites=[]
        with open(self.sf_sites) as fin_sites:
            for line in fin_sites:
                fields=line.rstrip().split()
                s_chrm=fields[0]
                i_pos=int(fields[1])
                l_sites.append((s_chrm, i_pos))

        s_sample_name, l_ori_bams=self._load_in_bams(self.sf_bam_list)
        #process the bam(s) one by one
        #with open(sf_new_bam_list,"w") as fout_bam_list:#
        #s_new_list=s_sample_name
        for sf_bam in l_ori_bams:#
            #sf_bam=line.rstrip()
            ori_samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_ref)
            m_chrm_id = self._get_chrm_id_name(ori_samfile)

            s_fname=os.path.basename(sf_bam)#
            sf_screenshot_bam=self._gnrt_new_bam_name(s_fname)#
            #s_new_list+=(" "+sf_screenshot_bam)
            scrn_bam = pysam.AlignmentFile(sf_screenshot_bam, "wb", template=ori_samfile)

            bam_info = BamInfo(sf_bam, self.sf_ref)
            #for each bam file, we need to calculate the basic information, especially the insert size
            x_basic_info = X_BasicInfo(self.sf_wfolder, self.n_jobs, self.sf_ref)
            search_win = 500
            f_ave_cov, rlth, mean_is, std_var = x_basic_info.collect_basic_info_one_sample(sf_bam, self.sf_ref,search_win )
            i_max_is=int(mean_is+3*std_var)
            for (s_chrm, i_ins_pos) in l_sites:
                b_with_chr = bam_info.is_chrm_contain_chr()
                chrm_in_bam = self._process_chrm_name(b_with_chr, s_chrm)
                start_pos = i_ins_pos - i_max_is
                if start_pos <= 0:
                    start_pos = 1
                end_pos = i_ins_pos + i_max_is
                for algnmt in ori_samfile.fetch(chrm_in_bam, start_pos, end_pos):#
                    ##here need to skip the secondary and supplementary alignments?
                    if algnmt.is_secondary or algnmt.is_supplementary:
                        continue
                    if algnmt.is_duplicate == True:  ##duplciate
                        continue
                    if algnmt.is_unmapped == True:  ##unmapped
                        continue
                    l_cigar = algnmt.cigar
                    if len(l_cigar) < 1:  # wrong alignment
                        continue
                    #
                    mapq = algnmt.mapping_quality
                    map_pos = algnmt.reference_start
                    mate_chrm = '*'
                    mate_pos = 0
                    if (algnmt.next_reference_id in m_chrm_id) and (algnmt.mate_is_unmapped == False) \
                            and (algnmt.next_reference_id >= 0):
                        mate_chrm = algnmt.next_reference_name
                        mate_pos = algnmt.next_reference_start

                    ####
                    b_qualified_rd=False
                    #if clipped reads, and clip position is within +/-20bp, then keep the reads
                    if l_cigar[0][0] == 4:  #left clipped
                        if mapq<global_values.MINIMUM_CLIP_MAPQ:
                            continue
                        if abs(i_ins_pos-map_pos)<=i_clip_extnd:
                            b_qualified_rd=True
                    elif l_cigar[-1][0] == 4:  # right clipped #
                        ##calculate the exact clip position
                        if mapq<global_values.MINIMUM_CLIP_MAPQ:
                            continue
                        for (type, lenth) in l_cigar[:-1]:
                            if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                                continue
                            else:
                                map_pos += lenth
                        if abs(i_ins_pos-map_pos)<=i_clip_extnd:
                            b_qualified_rd=True
                    ####
                    if mate_chrm == "*":
                        continue
                    xchrom = XChromosome()  ###decoy seuqence and contigs are not interested
                    if xchrom.is_decoy_contig_chrms(mate_chrm) == True:
                        continue

                    b_anchor_hard_clip = False
                    if l_cigar[-1][0] == 5 or l_cigar[0][0] == 5:
                        b_anchor_hard_clip = True
                    #next check discordant pairs
                    if (mapq > global_values.MINIMUM_DISC_MAPQ) and (b_anchor_hard_clip == False):
                        ## here only collect the read names for discordant reads, later will re-align the discordant reads
                        if self.is_discordant(chrm_in_bam, map_pos, mate_chrm, mate_pos, global_values.DISC_THRESHOLD) == True:
                            if abs(map_pos - i_ins_pos) <= i_max_is:#
                                b_qualified_rd=True
####
                    if b_qualified_rd==True:
                        scrn_bam.write(algnmt)#write the alignmt to the screenshot bam
            scrn_bam.close()
            if os.path.isfile(sf_screenshot_bam+".bai")==False:
                self._index_algnmt(sf_screenshot_bam)
            ori_samfile.close()
####
####
####
    ####Need to be merged later~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ####This is drectly copies from x_clip_disc_filter.py,
    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, b_tmplt_with_chr, chrm):#
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        if b_tmplt_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_tmplt_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_tmplt_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

####
    def _get_chrm_id_name(self, samfile):
        m_chrm = {}
        references = samfile.references
        for schrm in references:
            chrm_id = samfile.get_tid(schrm)
            m_chrm[chrm_id] = schrm
        m_chrm[-1] = "*"
        return m_chrm

    ####Need to be merged later~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ####This is drectly copies from x_clip_disc_filter.py
    ##whether a pair of read is discordant (for TEI only) or not
    def is_discordant(self, chrm, map_pos, mate_chrm, mate_pos, is_threshold):
        if chrm != mate_chrm:  ###of different chroms
            return True
        else:
            if abs(mate_pos - map_pos) > is_threshold:  # Of same chrom, but insert size are quite large
                return True
        return False
