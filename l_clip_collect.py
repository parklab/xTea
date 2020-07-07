##1/2/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

##
##03/04/19: Add function -- collect all the clip positions
##

import pysam
from multiprocessing import Pool
import global_values

def unwrap_self_collect_clip_parts_lrd(arg, **kwarg):
    return LRDClipReadInfo.collect_clip_contain_reads_one_site2(*arg, **kwarg)
def unwrap_self_collect_reads_cover_region_lrd(arg, **kwarg):
    return LRDClipReadInfo.collect_reads_cover_region_one_site(*arg, **kwarg)


class LRDClipReadInfo():
    def __init__(self, n_jobs, sf_bam="", sf_ref="", s_working_folder=""):
        self.n_jobs = n_jobs

        self.sf_bam = sf_bam
        self.sf_reference = sf_ref
        self.working_folder = s_working_folder

    ####Given sites, collect the clipped sequence (or contained insertions) for the sites
    def collect_seqs_for_sites_in_parallel(self, l_sites):
        #each record in format: (sf_bam, ins_chrm, ins_pos, sf_tmp)
        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_clip_parts_lrd, zip([self] * len(l_sites), l_sites), 1)
        pool.close()
        pool.join()

    #give regions, collect the reads that cover each region
    def collect_seqs_cover_regions_in_parallel(self, l_sites):#
        #each record in format: (sf_bam, ins_chrm, copy_pos_start, copy_pos_end, sf_tmp)
        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_reads_cover_region_lrd, zip([self] * len(l_sites), l_sites), 1)
        pool.close()
        pool.join()
####

    ####For each site, collect the clipped and contained sequences
    def collect_clip_contain_reads_one_site(self, record):
        sf_bam=record[0]
        ins_chrm=record[1]
        ins_pos=record[2]
        wfolder=record[3]

        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)##
        s_site = "{0}_{1}".format(ins_chrm, ins_pos)
        sf_out = wfolder + "{0}.fa".format(s_site)  # to save the clipped long reads
        print "Collect clipped long reads for {0}:{1}".format(ins_chrm, ins_pos)
        with open(sf_out, "w") as fout_clip:
            ins_pos_start=ins_pos - global_values.LRD_CLIP_SEARCH_WIN
            ins_pos_end=ins_pos + global_values.LRD_CLIP_SEARCH_WIN
            if ins_pos_start<0:
                print ins_chrm, ins_pos, "very small start!"
                ins_pos_start=0
            for algnmt in samfile.fetch(ins_chrm, ins_pos_start, ins_pos_end):
                if algnmt.is_unmapped == True:  #unmapped
                    continue
####
                l_cigar = algnmt.cigar
                map_pos = algnmt.reference_start  # the mapping position
                query_seq = algnmt.query_sequence
                query_name = algnmt.query_name

                if query_seq == None:
                    continue
                seq_lenth = len(query_seq)
                lclip_len = 0
                rclip_len = 0
                if l_cigar[0][0] == 4:# left-clip
                    brkpnt = ins_pos - global_values.BRKPNT_CK_WIN
                    lclip_len = l_cigar[0][1]
                    clip_pos = map_pos
                    if lclip_len > global_values.LRD_MIN_CLIP_LTH:
                        # print "left", algnmt.query_name, brkpnt, clip_pos, (clip_pos-lclip_len)
                        if brkpnt < clip_pos and brkpnt > (clip_pos - lclip_len):
                            clip_seq = query_seq[:lclip_len+global_values.LRD_EXTND_LEN]
                            shead = ins_chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR \
                                    + query_name + "L"
                            fout_clip.write(">" + shead + "\n")
                            fout_clip.write(clip_seq + "\n")

                if l_cigar[-1][0] == 4:  # right clipped
                    brkpnt = ins_pos + global_values.BRKPNT_CK_WIN
                    rclip_len = l_cigar[-1][1]
                    if rclip_len > global_values.LRD_MIN_CLIP_LTH:
                        clip_pos = map_pos + seq_lenth - lclip_len - rclip_len #this is not accuracy
                        # print "right", algnmt.query_name, brkpnt, map_pos, seq_lenth, lclip_len, clip_pos, (clip_pos+rclip_len)
                        if brkpnt > clip_pos and brkpnt < (clip_pos + rclip_len):
                            # print "right-clip", algnmt.query_name
                            clip_seq = query_seq[-1 * (rclip_len+global_values.LRD_EXTND_LEN):]
                            shead = ins_chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR + \
                                    query_name + "R"
                            fout_clip.write(">" + shead + "\n")
                            fout_clip.write(clip_seq + "\n")

                #also check whether there is a long insertion
                l_ctn_seqs=self._check_contained_insertion_from_cigar(map_pos, l_cigar, query_seq, ins_pos)
                i_cnt=1
                for s_ctn_seq in l_ctn_seqs:
                    shead = ins_chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR \
                            + query_name + "N" +str(i_cnt)
                    fout_clip.write(">" + shead + "\n")
                    fout_clip.write(s_ctn_seq + "\n")
                    i_cnt+=1
        samfile.close()

####
####
    def _check_contained_insertion_from_cigar(self, map_pos, l_cigar, s_seq, ins_pos):
        ref_pos=map_pos
        seq_pos=0
        l_seqs=[]
        for (opn, lth) in l_cigar:
            if opn == 0:  # alignment match (M)
                ref_pos += lth
                seq_pos+=lth
            elif opn==2 or opn==3:#deletion (D) or skipped region (N)
                ref_pos+=lth
            elif opn==4: #soft-clip (S)
                seq_pos+=lth
            elif opn==5 or opn==6:#hard-clip (H) or padding (P)
                continue
            elif opn==7 or opn==8:#sequence match (=) or sequence mismatch (X)
                ref_pos+=lth
                seq_pos+=lth
            elif opn==1:#insertion
                if lth >= global_values.LRD_MIN_INS_LTH:
                    if abs(ref_pos-ins_pos)<=global_values.LRD_MAX_BRKPNT_CHK_REGION:
                        ####here we extend both end (2kbp by default)!!!!!!!!!!
                        istart = seq_pos - global_values.LRD_EXTND_LEN
                        if istart < 0:
                            istart = 0
                        iend = seq_pos + lth + global_values.LRD_EXTND_LEN
                        # print istart, iend
                        s_tmp = s_seq[istart:iend]
                        l_seqs.append(s_tmp)
                seq_pos+=lth
        return l_seqs

####
    ####this version, we search those clipped reads with clip postion within a region[a,b]
    ####For each site, collect the clipped and contained sequences
    def collect_clip_contain_reads_one_site2(self, record):#
        sf_bam = record[0]
        ins_chrm = record[1]
        ins_pos = record[2]
        wfolder = record[3]

        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        s_site = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, ins_pos)
        i_extend = global_values.LRD_EXTND_LEN
        sf_out = wfolder + "{0}_{1}.fa".format(s_site, i_extend)  # to save the clipped long reads
        print "Collect clipped long reads for {0}:{1}".format(ins_chrm, ins_pos)
        with open(sf_out, "w") as fout_clip:
            ins_pos_start = ins_pos - global_values.LRD_CLIP_SEARCH_WIN #by default 20k
            ins_pos_end = ins_pos + global_values.LRD_CLIP_SEARCH_WIN
            if ins_pos_start < 0:
                print ins_chrm, ins_pos, "very small start!"
                ins_pos_start = 0
            for algnmt in samfile.fetch(ins_chrm, ins_pos_start, ins_pos_end):
                if algnmt.is_unmapped == True:  # unmapped
                    continue
                if algnmt.mapping_quality < global_values.LRD_MIN_MAPQ:  # by default this is set to 20
                    continue
                if algnmt.is_duplicate == True:  ##duplciate
                    continue
                if algnmt.is_secondary or algnmt.is_supplementary:#one read, will only collect the sequence once
                    continue
####
                l_cigar = algnmt.cigar
                if len(l_cigar) < 1:  # wrong alignment
                    continue
                map_pos = algnmt.reference_start  # the mapping position
                query_seq = algnmt.query_sequence
                query_name1 = algnmt.query_name #in some cases, there are space or "\t" within the query name
                query_name_fields=query_name1.split()
                query_name="/".join(query_name_fields)

                if query_seq == None:
                    continue

                # also check whether there is a long insertion
                l_ctn_seqs, i_rclip_pos, n_mapped = self._check_contained_insertion_from_cigar2(map_pos, l_cigar,
                                                                                                query_seq, ins_pos)
                if n_mapped < global_values.LRD_MIN_MAP_LEN:  # require at least the mapped region is 1k (by default)
                    continue
                i_cnt = 1
                for s_ctn_seq in l_ctn_seqs:
                    shead = ins_chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR \
                            + query_name + "N" + str(i_cnt)
                    fout_clip.write(">" + shead + "\n")
                    fout_clip.write(s_ctn_seq + "\n")
                    i_cnt += 1

                seq_lenth = len(query_seq)
                lclip_len = 0
                rclip_len = 0
                if l_cigar[0][0] == 4:  # left-clip
                    #brkpnt = ins_pos - global_values.BRKPNT_CK_WIN
                    lclip_len = l_cigar[0][1]
                    clip_pos = map_pos
                    if lclip_len > global_values.LRD_MIN_CLIP_COLLECT_LTH:
                        # print "left", algnmt.query_name, brkpnt, clip_pos, (clip_pos-lclip_len)
                        #if brkpnt < clip_pos and brkpnt > (clip_pos - lclip_len):
                        if abs(clip_pos-ins_pos)<=global_values.LRD_MAX_BRKPNT_CHK_REGION:
                            shead = ins_chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR \
                                    + query_name + "L"

                            clip_seq = query_seq[:lclip_len + global_values.LRD_EXTND_LEN]
                            if global_values.LRD_EXTND_LEN<=0:
                                clip_seq=query_seq
                            fout_clip.write(">" + shead + "\n")
                            fout_clip.write(clip_seq + "\n")

                if l_cigar[-1][0] == 4:  # right clipped
                    #brkpnt = ins_pos + global_values.BRKPNT_CK_WIN
                    rclip_len = l_cigar[-1][1]
                    if rclip_len > global_values.LRD_MIN_CLIP_COLLECT_LTH:
                        #clip_pos = map_pos + seq_lenth - lclip_len - rclip_len  # this is not accuracy
                        # print "right", algnmt.query_name, brkpnt, map_pos, seq_lenth, lclip_len, clip_pos, (clip_pos+rclip_len)
                        #if brkpnt > clip_pos and brkpnt < (clip_pos + rclip_len):
                        if abs(i_rclip_pos-ins_pos)<=global_values.LRD_MAX_BRKPNT_CHK_REGION:
                            # print "right-clip", algnmt.query_name
                            clip_seq = query_seq[-1 * (rclip_len + global_values.LRD_EXTND_LEN):]
                            if global_values.LRD_EXTND_LEN<=0:
                                clip_seq=query_seq
                            shead = ins_chrm + global_values.SEPERATOR + str(ins_pos) + global_values.SEPERATOR + \
                                    query_name + "R"
                            fout_clip.write(">" + shead + "\n")
                            fout_clip.write(clip_seq + "\n")
        samfile.close()
    ####
    ####
    def collect_reads_cover_region_one_site(self, record):####
        sf_bam = record[0]
        ins_chrm = record[1]
        i_regin_start = record[2]
        i_regin_end = record[3]
        wfolder = record[4]

        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        s_site = "{0}{1}{2}".format(ins_chrm, global_values.SEPERATOR, i_regin_start)
        i_extend = global_values.LRD_EXTND_LEN
        sf_out = wfolder + "{0}_{1}.fa".format(s_site, i_extend)  # to save the clipped long reads
        print "Collect clipped long reads for {0}:{1}_{2}".format(ins_chrm, i_regin_start, i_regin_end)
        with open(sf_out, "w") as fout_reads:
            ins_pos_start = i_regin_start - global_values.LRD_CLIP_SEARCH_WIN  # by default 20k
            ins_pos_end = i_regin_end + global_values.LRD_CLIP_SEARCH_WIN
            if ins_pos_start < 0:
                print ins_chrm, i_regin_start, "very small start!"
                ins_pos_start = 0
            for algnmt in samfile.fetch(ins_chrm, ins_pos_start, ins_pos_end):#
                if algnmt.is_unmapped == True:  # unmapped
                    continue
                if algnmt.mapping_quality < global_values.LRD_MIN_MAPQ:  # by default this is set to 20
                    continue
                if algnmt.is_duplicate == True:  ##duplciate
                    continue
                if algnmt.is_secondary or algnmt.is_supplementary:  # one read, will only collect the sequence once
                    continue
                    ####
                l_cigar = algnmt.cigar
                if len(l_cigar) < 1:  # wrong alignment
                    continue
                map_pos = algnmt.reference_start  # the mapping position
                query_seq = algnmt.query_sequence
                query_name1 = algnmt.query_name  # in some cases, there are space or "\t" within the query name
                query_name_fields = query_name1.split()
                query_name = "/".join(query_name_fields)

                if query_seq == None:
                    continue

                seq_lenth = len(query_seq)
                lclip_len = 0
                rclip_len = 0

                map_end=map_pos+seq_lenth
                b_save=False
                i_regin_start1=i_regin_start-global_values.LRD_MAX_BRKPNT_CHK_REGION
                i_region_end1=i_regin_end+global_values.LRD_MAX_BRKPNT_CHK_REGION
                if map_pos>i_regin_start1 and map_pos<i_region_end1:
                    b_save=True
                if map_end>i_regin_start1 and map_end<i_region_end1:
                    b_save=True

                i_total_clip=0
                i_right_clip_pos = self._get_right_clip_pos(map_pos, l_cigar)
                i_left_clip_pos =-1
                if l_cigar[0][0] == 4:  # left-clip
                    # brkpnt = ins_pos - global_values.BRKPNT_CK_WIN
                    lclip_len = l_cigar[0][1]
                    i_total_clip+=lclip_len
                    i_left_clip_pos=map_pos

                if l_cigar[-1][0] == 4:  # right clipped
                    # brkpnt = ins_pos + global_values.BRKPNT_CK_WIN
                    rclip_len = l_cigar[-1][1]
                    i_total_clip+=rclip_len

                n_mapped=seq_lenth-i_total_clip
                if n_mapped < global_values.LRD_MIN_MAP_LEN:  # require at least the mapped region is 1k (by default)
                    continue
                shead = ins_chrm + global_values.SEPERATOR + str(i_regin_start) + global_values.SEPERATOR \
                        + query_name + "N"
                fout_reads.write(">" + shead + "\n")
                fout_reads.write(query_seq + "\n")
        samfile.close()
####
    ####
    def _check_contained_insertion_from_cigar2(self, map_pos, l_cigar, s_seq, ins_pos):
        ref_pos = map_pos
        seq_pos = 0
        l_seqs = []
        i_rclip_pos = -1
        n_mapped = 0
        s_full_seq=""
        b_full=False
        for (opn, lth) in l_cigar:
            if opn == 0:  # alignment match (M)
                ref_pos += lth
                seq_pos += lth
                n_mapped += lth  # accumulate the mapped length
            elif opn == 2 or opn == 3:  # deletion (D) or skipped region (N)
                ref_pos += lth
            elif opn == 4:  # soft-clip (S)
                seq_pos += lth
                if len(l_cigar) > 0 and l_cigar[-1][0] == 4:  # right clip
                    i_rclip_pos = ref_pos
            elif opn == 5 or opn == 6:  # hard-clip (H) or padding (P)
                continue
            elif opn == 7 or opn == 8:  #sequence match (=) or sequence mismatch (X)
                ref_pos += lth
                seq_pos += lth
                n_mapped += lth
            elif opn == 1:  # insertion
                if lth >= global_values.LRD_MIN_INS_LTH:
                    if abs(ref_pos - ins_pos) <= global_values.LRD_MAX_BRKPNT_CHK_REGION:
                        ####here we extend both end (2kbp by default)!!!!!!!!!!
                        if global_values.LRD_EXTND_LEN<=0:
                            #l_seqs.append(s_seq)
                            b_full=True
                            s_full_seq=s_seq
                        else:
                            istart = seq_pos - global_values.LRD_EXTND_LEN
                            if istart < 0:
                                istart = 0
                            iend = seq_pos + lth + global_values.LRD_EXTND_LEN
                            # print istart, iend
                            s_tmp = s_seq[istart:iend]
                            l_seqs.append(s_tmp)
                seq_pos += lth
        if b_full==True:#make sure the whole read is only saved for one time
            l_seqs.append(s_full_seq)
        return l_seqs, i_rclip_pos, n_mapped

####
    def _get_right_clip_pos(self, map_pos, l_cigar):
        ref_pos = map_pos
        seq_pos = 0
        i_rclip_pos = -1
        n_mapped = 0
        for (opn, lth) in l_cigar:
            if opn == 0:  # alignment match (M)
                ref_pos += lth
                seq_pos += lth
                n_mapped += lth  # accumulate the mapped length
            elif opn == 2 or opn == 3:  # deletion (D) or skipped region (N)
                ref_pos += lth
            elif opn == 4:  # soft-clip (S)
                seq_pos += lth
                if len(l_cigar) > 0 and l_cigar[-1][0] == 4:  # right clip
                    i_rclip_pos = ref_pos
            elif opn == 5 or opn == 6:  # hard-clip (H) or padding (P)
                continue
            elif opn == 7 or opn == 8:  # sequence match (=) or sequence mismatch (X)
                ref_pos += lth
                seq_pos += lth
                n_mapped += lth
        return i_rclip_pos
####