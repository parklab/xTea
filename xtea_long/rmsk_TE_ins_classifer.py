##10/17/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

#this module define the full length SVA from repeatmasker output

'''
Rules:
1) The first row of the consensus contain ~60bp (ccctct)n heximer, if the masked sequence should start from this range
2) The length should longer than 700bp
3) The sequence are masked as a unique type
4) Contain polyA
5) New subfamily with Alu at 5' or 3' (or both) side.

Output:
1) How many are masked (>85% of the whole sequence, pay attention to the new type) as SVA (any type);
2) How many are masked as specific type;
3) How many are failed to be masked as SVA;
4) Calculate the length of the masked ones;

SVA subfamily consensus sequence:
SVA-A: VNTR region:[420-1080], total length: 1640
SVA-B: VNTR region:[360-800], total length: 1383
SVA-C: VNTR region:[420-800], total length: 1384
SVA-D: VNTR region:[420-800], total length: 1386
SVA-E: VNTR region:[420-800], total length: 1382
SVA-F: VNTR region:[410-790], total length: 1375
SVA-F1: VNTR region:[TBD]
'''

import os
import sys
import pysam
from x_contig import *
from l_dimorphic_HERV import *
from global_values import *
from x_polyA import *
####

# 7288   1.3  0.1  0.1  AK1_chr10_68475191     32213340 32214177 (15915718) +  SVA_F Retroposon/SVA 5312 6155    (0) 4482983
# 4741   0.6  0.0  0.0  AK1_chr10_68475191     17916705 17917238 (30212657) C  SVA_F Retroposon/SVA (0) 6155   5622 4462908
class SVAClassifier():#
    def __init__(self,sf_hex, sf_wfolder, n_jobs):
        #self.sf_annotation=sf_anno
        self.m_vntr={}
        #each record [hexamer_length, vntr_start, vntr_end, cns_polyA_length, total_length]
        #self.m_vntr['SVA_A']=(65, 420,1080, 20, 1640) #this is from consensus sequence
        self.m_vntr['SVA_A'] = (65, 420, 800, 20, 1387)
        self.m_vntr['SVA_B'] = (70, 360, 800, 20, 1383)
        self.m_vntr['SVA_C'] = (70, 420, 800, 20, 1384)
        self.m_vntr['SVA_D'] = (70, 420, 800, 20, 1386)
        self.m_vntr['SVA_E'] = (70, 420, 800, 20, 1382)
        self.m_vntr['SVA_F'] = (70, 410, 790, 20, 1375)
        self.vntr_start=360
        self.vntr_end=900
        self.sf_wfolder=sf_wfolder
        self.n_jobs=n_jobs
        self.i_hexamer_lenth=60
        self.m_hexamer={}
        self.m_hexamer_rc={}
        self.TYPE_AMBIGOUS="Unknown"
        self.load_in_hexamer(sf_hex)#
        self.MIN_TD_LEN=35#minimal transduction length
        self.TD5=global_values.SEPERATOR+"5TD"
        self.TD3=global_values.SEPERATOR+"3TD"
        self.MAST2=global_values.SEPERATOR+"MAST2"
        self.MAST2_only="MAST2"
        self.FULL_COPY="full_copy"
        self.TRUNCATED="truncated"
        self.NOT_ALU="NULL"
        self.NULL="NULL"
        self.ALU="Alu"
        self.simple_rep="Simple_repeat"
        self.mid_segmt_start=-1
        self.mid_segmt_end = -1
        self.TRANSDCT_UNIQ_MAPQ=30
        self.SEPRATOR="~"#by default, this is step as "~"
        pa=PolyA()
        self.m_polyA=pa.get_pre_defined_polyA_in_rmsk()
        self.m_polyT=pa.get_pre_defined_polyT_in_rmsk()
####

####
    def collect_full_length_ref_SVA(self, sf_rmsk, i_max_dist, sf_ref, sf_out_sites):
        m_non_rc=self.collect_full_length_cannonical_SVA_from_ref_rmsk(sf_rmsk, i_max_dist)
        m_rc=self.collect_full_length_cannonical_rc_SVA_from_ref_rmsk(sf_rmsk, i_max_dist)
        m_non_rc_full, l_non_rc_rg, l_non_rc_sine_r, m_non_rc_trnc=self.slct_ref_full_length_SVA(m_non_rc, self.m_hexamer)
        b_rc=True
        #each rcd in format: chrm~start~end : (chrm, i_start, i_end, sine_r_start, sine_r_end, s_sub_family)
        m_rc_full, l_rc_rg, l_rc_sine_r, m_rc_trnc=self.slct_ref_full_length_SVA(m_rc, self.m_hexamer_rc, b_rc)

        xref = XReference()
        l_non_rc_seqs=xref.get_ref_seqs_of_sites(sf_ref, l_non_rc_sine_r)
        l_rc_seqs=xref.get_ref_seqs_of_sites(sf_ref, l_rc_sine_r)
        ####
        sf_fa=sf_out_sites+".fa"
        with open(sf_out_sites, "w") as fout_sites, open(sf_out_sites+".all_sites", "w") as \
                fout_all_sites, open(sf_fa, "w") as fout_fa:
            for l_one_site, s_seq in zip(l_non_rc_rg, l_non_rc_seqs):
                chrm=l_one_site[0]
                i_start=l_one_site[1]
                i_end=l_one_site[2]
                s_id=chrm+global_values.SEPERATOR+str(i_start)+global_values.SEPERATOR+str(i_end)
                (chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family)=m_non_rc_full[chrm][i_start]
                s_info="%s\t%s\t%s\t%s\t%s\t+\t%s\n" % (chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family)
                fout_sites.write(s_info)
                fout_all_sites.write(s_info)

                if len(s_seq)<400:
                    continue
                s_fa_id=s_id+global_values.SEPERATOR+s_sub_family+global_values.SEPERATOR+"SINE_R"
                fout_fa.write(">"+s_fa_id+"\n")
                fout_fa.write(s_seq+"\n")

            for l_one_site, s_seq in zip(l_rc_rg, l_rc_seqs):
                chrm=l_one_site[0]
                i_start=l_one_site[1]
                i_end=l_one_site[2]
                s_id=chrm+global_values.SEPERATOR+str(i_start)+global_values.SEPERATOR+str(i_end)
                (chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family)=m_rc_full[chrm][i_start]#
                s_info="%s\t%s\t%s\t%s\t%s\tC\t%s\n" % (chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family)
                fout_sites.write(s_info)
                fout_all_sites.write(s_info)

                if len(s_seq)<400:
                    continue
                s_fa_id = s_id + global_values.SEPERATOR + s_sub_family + global_values.SEPERATOR + "SINE_R_RC"
                fout_fa.write(">" + s_fa_id + "\n")
                fout_fa.write(self._gnrt_seq_rc(s_seq) + "\n")

            #save those truncated to "fout_all_sites"
            for chrm in m_non_rc_trnc:
                for i_start in m_non_rc_trnc[chrm]:
                    (chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family) = m_non_rc_trnc[chrm][i_start]
                    s_info = "%s\t%s\t%s\t%s\t%s\t+\t%s\n" % (
                    chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family)
                    fout_all_sites.write(s_info)
            for chrm in m_rc_trnc:
                for i_start in m_rc_trnc[chrm]:
                    (chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family) = m_rc_trnc[chrm][i_start]
                    s_info = "%s\t%s\t%s\t%s\t%s\tC\t%s\n" % (
                        chrm, i_start1, i_end1, sine_r_start, sine_r_end, s_sub_family)
                    fout_all_sites.write(s_info)
####
####
####
    #Check the reference genome annotation, for each "SVA" copy: (SVA fusion cases are not considered)
    #1. If a single SVA, and hit the hexamer and SINE-R region, then view it as a hit;
    #2. If several SVA (and heximer, polyA) records that are close to each other, and hit hexamer, SINE-R, then collect.
    #3. This is only work for the non reverse complementary ones
    def collect_full_length_cannonical_SVA_from_ref_rmsk(self, sf_rmsk, i_max_dist):
        m_candidates={}#
        with open(sf_rmsk) as fin_rmsk:
            pre_chrm=""
            pre_start=0
            pre_end=0
            pre_b_sva=False
            pre_b_hexamer=False
            pre_s_sva=""
            pre_b_polyA=False
            l_tmp=[]
            for line in fin_rmsk:
                fields = line.split()
                if len(fields) < 4:
                    continue
                if ("SW" in fields[0]) or ("score" in fields[0]):
                    continue
                chrm = fields[4]
                start_pos = int(fields[5])
                end_pos = int(fields[6])
                sub_type = fields[9]
                s_super_family=fields[10]

                i_seq_dist_from_end = int(fields[7][1:-1])
                b_rc = False
                if fields[8] == "C":
                    b_rc = True
                csn_start = -1
                csn_end = -1
                if b_rc == True:
                    csn_start = int(fields[13])
                    csn_end = int(fields[12])
                else:
                    csn_start = int(fields[11])
                    csn_end = int(fields[12])

                b_sva=False
                b_hexamer=False
                b_polyA=False
                s_sva="none"
                if self._is_sub_type_hexamer(sub_type, s_super_family)==True:
                    b_hexamer=True
                    b_sva=True
                    s_sva = "sva"
                elif sub_type in self.m_vntr:
                    if b_rc==False:#only consider the non-reverse complementary
                        b_sva=True
                        s_sva = "sva"
                elif self._is_polyA(sub_type, s_super_family)==True:
                    b_polyA=True
                    b_sva = True
                    s_sva = "sva"

                if b_sva==True:#
                    #if previous is also True
                    if pre_b_sva == True:
                        if (chrm == pre_chrm) and ((start_pos - pre_end) < i_max_dist) and (pre_s_sva==s_sva):
                            l_tmp.append((start_pos, end_pos, sub_type, csn_start, csn_end, b_rc, i_seq_dist_from_end))
                        else:
                            if pre_b_hexamer == False:
                                #save the previous one
                                s_tmp_id=chrm+global_values.SEPERATOR+str(pre_start)
                                m_candidates[s_tmp_id]=[]#note, this position is not the correct position
                                for tmp_rcd in l_tmp:
                                    m_candidates[s_tmp_id].append(tmp_rcd)
                                    #m_candidates_rc[s_tmp_id]=b_rc
                            #clean old one
                            del l_tmp[:]
                            #save new one
                            l_tmp.append((start_pos, end_pos, sub_type, csn_start, csn_end, b_rc, i_seq_dist_from_end))
                    else:
                        if b_polyA==True:#polyA should not at the front part!!!!
                            b_sva=False
                            s_sva="none"
                        else:
                            # clean old one
                            del l_tmp[:]
                            # save new one
                            l_tmp.append((start_pos, end_pos, sub_type, csn_start, csn_end, b_rc, i_seq_dist_from_end))
                else:
                    if pre_b_sva == True:#
                        if pre_b_hexamer==False:#if previous just hexamer, then not save, otherwise save
                            # save the previous one
                            s_tmp_id = chrm + global_values.SEPERATOR + str(pre_start)
                            m_candidates[s_tmp_id] = []
                            for tmp_rcd in l_tmp:
                                m_candidates[s_tmp_id].append(tmp_rcd)
                                #m_candidates_rc[s_tmp_id] = b_rc
                        # clean old one
                        del l_tmp[:]
####
                pre_s_sva=s_sva
                pre_b_sva=b_sva
                pre_b_hexamer=b_hexamer
                pre_chrm=chrm
                pre_start=start_pos
                pre_end=end_pos
                pre_b_polyA=b_polyA
        return m_candidates


    def _is_sub_type_hexamer(self, sub_type, s_super_family):
        if (sub_type in self.m_hexamer) or (s_super_family == self.simple_rep and "CCC" in sub_type):
            return True
        else:
            return False

    #
    def _is_polyA(self, sub_type, s_super_family):
        if (sub_type in self.m_polyA) or (s_super_family == self.simple_rep and "AAA" in sub_type):
            return True
        else:
            return False

####
    # Check the reference genome annotation, for each "SVA" copy: (SVA fusion cases are not considered)
    # 1. If a single SVA, and hit the hexamer and SINE-R region, then view it as a hit;
    # 2. If several SVA (and heximer, polyA) records that are close to each other, and hit hexamer, SINE-R, then collect.
    # 3. This is only work for the reverse complementary ones
    def collect_full_length_cannonical_rc_SVA_from_ref_rmsk(self, sf_rmsk, i_max_dist):
        m_candidates = {}
        with open(sf_rmsk) as fin_rmsk:
            pre_chrm = ""
            pre_start = 0
            pre_end = 0
            pre_b_sva = False
            pre_b_hexamer = False
            pre_b_polyT=False
            pre_s_sva = ""
            l_tmp = []
            for line in fin_rmsk:
                fields = line.split()
                if len(fields) < 4:
                    continue
                if ("SW" in fields[0]) or ("score" in fields[0]):
                    continue
                chrm = fields[4]
                start_pos = int(fields[5])
                end_pos = int(fields[6])
                sub_type = fields[9]
                s_super_family=fields[10]

                i_seq_dist_from_end = int(fields[7][1:-1])
                b_rc = False
                if fields[8] == "C":
                    b_rc = True
                csn_start = -1
                csn_end = -1
                if b_rc == True:
                    csn_start = int(fields[13])
                    csn_end = int(fields[12])
                else:
                    csn_start = int(fields[11])
                    csn_end = int(fields[12])

                b_sva = False
                b_hexamer = False
                b_polyT=False
                s_sva = "none"
                if self._is_sub_type_hexamer_rc(sub_type, s_super_family)==True:
                    b_hexamer = True
                    b_sva = True
                    s_sva = "sva"
                elif sub_type in self.m_vntr:
                    if b_rc == True:  # only consider the reverse complementary cases
                        b_sva = True
                        s_sva = "sva"
                elif self._is_polyT(sub_type, s_super_family)==True:
                    b_sva = True
                    s_sva = "sva"
                    b_polyT=True

                if b_sva == True:
                    # if previous is also True
                    if pre_b_sva == True:
                        if b_polyT==True:#polyT should not at the tail
                            b_sva=False
                            s_sva="none"
                        else:
                            if (chrm == pre_chrm) and ((start_pos - pre_end) < i_max_dist) and (pre_s_sva == s_sva):##
                                l_tmp.append((start_pos, end_pos, sub_type, csn_start, csn_end, b_rc, i_seq_dist_from_end))
                            else:
                                if pre_b_hexamer == False or len(l_tmp)>1:
                                    # save the previous one
                                    s_tmp_id = chrm + global_values.SEPERATOR + str(pre_start)
                                    m_candidates[s_tmp_id] = []  # note, this position is not the correct position
                                    for tmp_rcd in l_tmp:
                                        m_candidates[s_tmp_id].append(tmp_rcd)
                                        # m_candidates_rc[s_tmp_id]=b_rc
                                # clean old one
                                del l_tmp[:]
                                # save new one
                                l_tmp.append((start_pos, end_pos, sub_type, csn_start, csn_end, b_rc, i_seq_dist_from_end))
                    else:
                        # clean old one
                        del l_tmp[:]
                        # save new one
                        l_tmp.append((start_pos, end_pos, sub_type, csn_start, csn_end, b_rc, i_seq_dist_from_end))
                else:
                    if pre_b_sva == True:
                        if pre_b_hexamer == False or len(l_tmp)>1:  # if previous just hexamer, then not save, otherwise save
                            # save the previous one
                            s_tmp_id = chrm + global_values.SEPERATOR + str(pre_start)
                            m_candidates[s_tmp_id] = []
                            for tmp_rcd in l_tmp:
                                m_candidates[s_tmp_id].append(tmp_rcd)
                        # clean old one
                        del l_tmp[:]
                pre_s_sva = s_sva
                pre_b_sva = b_sva
                pre_b_hexamer = b_hexamer
                pre_b_polyT=b_polyT
                pre_chrm = chrm
                pre_start = start_pos
                pre_end = end_pos
        return m_candidates

    def _is_sub_type_hexamer_rc(self, sub_type, s_super_family):
        if (sub_type in self.m_hexamer_rc) or (s_super_family == self.simple_rep and "GGG" in sub_type):
            return True
        else:
            return False
    #
    def _is_polyT(self, sub_type, s_super_family):
        if (sub_type in self.m_polyT) or (s_super_family == self.simple_rep and "TTT" in sub_type):
            return True
        else:
            return False

####
    def slct_ref_full_length_SVA(self, m_all_candidates, m_hexamer, b_rc=False):
        m_full_length={}
        m_truncated={}
        l_regions=[]
        l_sine_r=[]
        for s_tmp_id in m_all_candidates:
            b_with_hex = True
            b_with_vntr = True
            s_id_fields=s_tmp_id.split(global_values.SEPERATOR)
            chrm=s_id_fields[0]
            l_segmts=m_all_candidates[s_tmp_id]
            rslt_rcd, pos_rcd = self.classify_cannonical_SVA(l_segmts, m_hexamer, b_rc)
            if rslt_rcd is None:
                continue
            n_hex=pos_rcd[-2]
            if n_hex<=0:#not full length
                b_with_hex=False
                #continue
            n_vntr=pos_rcd[-1]#skip those only with hexamer
            if n_vntr<=0:
                b_with_vntr=False
                #continue
            #find the start and end position on consensus
            i_start=1000000000000
            i_end=0
            for rcd in l_segmts:
                if rcd[0]<i_start:
                    i_start=rcd[0]
                if rcd[1]>i_end:
                    i_end=rcd[1]
            s_sub_family=rslt_rcd[4]
            b_hit_hex=rslt_rcd[5]
            if b_hit_hex==False:
                b_with_hex=False
                #continue

            #s_id=chrm+global_values.SEPERATOR+str(i_start)+global_values.SEPERATOR+str(i_end)
            sine_r_start,sine_r_end=self.parse_SINE_R_region(l_segmts, b_rc)

            if b_with_hex==True and b_with_vntr==True:#full length
                if chrm not in m_full_length:
                    m_full_length[chrm]={}
                m_full_length[chrm][i_start]=(chrm, i_start, i_end, sine_r_start, sine_r_end, s_sub_family)
                l_sine_r.append((chrm, sine_r_start, sine_r_end))
                l_regions.append((chrm, i_start, i_end))
            else:#Truncated
                if chrm not in m_truncated:
                    m_truncated[chrm]={}
                m_truncated[chrm][i_start]=(chrm, i_start, i_end, sine_r_start, sine_r_end, s_sub_family)
        return m_full_length, l_regions, l_sine_r, m_truncated
####
####
    ####
    def classify_SVA_from_rmsk_output(self, sf_annotation, m_contigs, sf_mast2, sf_ref, f_min_sva_cutoff, f_mast2_cutoff, sf_out):
        m_records = {}  # save the records
        with open(sf_annotation) as fin_rmsk:
            for line in fin_rmsk:
                fields = line.split()
                if len(fields) < 4:
                    continue
                if ("SW" in fields[0]) or ("score" in fields[0]):
                    continue

                s_id = fields[4]  # in format: HG02106_chr6_24648990
                tmp_fields = s_id.split("_")
                if len(tmp_fields) > 3:  # skip those duplicated ones (like HG02106_chr6_24648990_1)
                    continue
                start_pos = int(fields[5])
                end_pos = int(fields[6])
                i_seq_dist_from_end = int(fields[7][1:-1])
                b_rc = False
                if fields[8] == "C":
                    b_rc = True
                sub_type = fields[9]#
                csn_start = -1
                csn_end = -1
                if b_rc == True:
                    csn_start = int(fields[13])
                    csn_end = int(fields[12])
                else:
                    csn_start = int(fields[11])
                    csn_end = int(fields[12])
                if s_id not in m_records:
                    m_records[s_id] = []
                m_records[s_id].append((start_pos, end_pos, sub_type, csn_start, csn_end, b_rc, i_seq_dist_from_end))

        n_ambigous = 0
        n_qualified_sva = 0
        sf_td_mast2_fa = sf_out + ".fa"
        sf_stage1_out=sf_out+".stage1_out"
        b_with_seq = True
        m_classified, m_cns_classified, n_ambigous, n_qualified_sva = self.classify_records(m_records, m_contigs, f_min_sva_cutoff,
                                                                            sf_stage1_out, sf_td_mast2_fa, b_with_seq)
        #align the transduction, middle potential MAST2 region to MAST2 gene
        sf_algnmt=self.sf_wfolder+"td_mast2_region_algn_2_MAST2.bam"
        self.algn_seq_to_MAST2(sf_mast2, sf_td_mast2_fa, sf_algnmt)
        sf_not_algned=sf_out + ".not_algn_MAST2.fa"
        #each record of m_hit_mast2: (map_pos, i_map_len)
        m_hit_mast2 = self.parse_MAST2_algnmt(sf_algnmt, sf_mast2, f_mast2_cutoff, sf_not_algned)

        #realign the unmapped to human reference genome
        sf_algnmt2=self.sf_wfolder+"td_mid_region_algn_2_ref.bam"
        self.algn_TD_mid_contig_to_ref(sf_ref, sf_not_algned, sf_algnmt2)
        #parse out those mapped, and those hit other exons
        m_hit_ref=self.parse_second_algnmt_to_ref(sf_algnmt2, sf_ref, f_mast2_cutoff)
        #print m_hit_ref
        #re-classify based on the re-aligned information
        #print m_cns_classified
        m_final_list=self.re_classify_INS_type(m_classified, m_cns_classified, m_hit_mast2, m_hit_ref, sf_out)#
        return m_final_list
####
####
    def get_ins_sine_r_seq(self, s_id, l_segmts):
        pass

    ####
    def classify_records(self, m_records, m_contigs, f_min_sva_cutoff, sf_rslts, sf_target_seq, b_with_seq):##
        m_seq_rslts = {}
        m_cns_pos_rslts={}
        n_ambigous = 0
        n_qualified_sva = 0
        fout_tmp_fa = None
        if b_with_seq == True:
            fout_tmp_fa = open(sf_target_seq, "w")
        with open(sf_rslts, "w") as fout, open(sf_rslts+".sine_r.fa", "w") as fout_sine_r:
            fout.write("sample_pos_id,full_or_truncated,pre_Alu,SVA_subfamily,suffix_Alu,"
                       "5_trnsdct_len,MAST2_pos,Alu_like_pos,Sine_R_pos,3_trnsdct_len\n")
            for s_id in m_records:
                if s_id not in m_contigs:
                    print(("%s not in m_contigs\n" % s_id))
                    continue
                s_contig_seq = m_contigs[s_id]#contig seq
                l_segmts = m_records[s_id]
                seq_rcd, cns_pos_rcd, s_type, i_tmp_len, t_ins_sine_r = self.classify_masked_record(s_id, l_segmts)
                if t_ins_sine_r is not None:
                    i_sine_r_start=int(t_ins_sine_r[1])
                    i_sine_r_end = int(t_ins_sine_r[2])
                    fout_sine_r.write(">"+s_id+"_SINE_R\n")
                    fout_sine_r.write(s_contig_seq[i_sine_r_start:i_sine_r_end+1]+"\n")
####
                if s_type is not None:
                    n_ambigous += 1
                elif seq_rcd is not None:
                    n_qualified_sva += 1

                    idx_extra_hex, s_pre_Alu, i_MAST2_start, i_MAST2_end, s_sub_family, b_hexamer, \
                    i_masked_lenth, s_suf_Alu, idx_sva, b_sva_rc = seq_rcd

                    if float(i_masked_lenth)/float(len(s_contig_seq)) < f_min_sva_cutoff:
                        continue

                    #get the transduction sequence

                    i_5_unmask_len, i_3_unmask_len, i_extra_polyA, s_5_seq, s_3_seq = \
                        self.parse_TD_lenth_seq(l_segmts, s_contig_seq, idx_extra_hex, idx_sva, b_sva_rc)

                    # if "chr9_33130550" in s_id:
                    #     print idx_extra_hex,idx_sva,b_sva_rc,i_5_unmask_len,i_3_unmask_len,"test!!!"

                    s_full = self.FULL_COPY
                    if b_hexamer == False:
                        s_full = self.TRUNCATED
                    sinfo = s_id + "," + s_full + "," + str(
                        i_5_unmask_len) + "," + s_pre_Alu + "," + s_sub_family + "," \
                            + s_suf_Alu + "," + str(i_3_unmask_len)
                    fout.write(sinfo + "\n")
                    m_seq_rslts[s_id] = (s_full, i_5_unmask_len, s_pre_Alu, s_sub_family, s_suf_Alu,
                                         i_3_unmask_len, i_masked_lenth+i_extra_polyA, b_sva_rc)
                    m_cns_pos_rslts[s_id]=cns_pos_rcd
                    if b_with_seq == True:
                        if i_5_unmask_len > self.MIN_TD_LEN:  #save 5-TD
                            s_tmp_5 = s_id + self.TD5
                            fout_tmp_fa.write(">" + s_tmp_5 + "\n" + str(s_5_seq) + "\n")
                        if (i_MAST2_start > 0 and i_MAST2_end > 0):  # save MAST2
                            s_tmp_mast2 = s_id + self.MAST2
                            s_seq_mast2 = s_contig_seq[i_MAST2_start:i_MAST2_end]
                            fout_tmp_fa.write(">" + s_tmp_mast2 + "\n" + s_seq_mast2 + "\n")
                        if i_3_unmask_len > self.MIN_TD_LEN:  # save 3-TD
                            s_tmp_3 = s_id + self.TD3
                            fout_tmp_fa.write(">" + s_tmp_3 + "\n" + str(s_3_seq) + "\n")
        if b_with_seq == True:
            fout_tmp_fa.close()
        return m_seq_rslts, m_cns_pos_rslts, n_ambigous, n_qualified_sva
####

    ####
    def algn_seq_to_MAST2(self, sf_mast2, sf_seq, sf_algnmt):
        xtea_contig = XTEContig(self.sf_wfolder, self.n_jobs)
        xtea_contig.align_short_contigs_2_cns_minimap2(sf_mast2, sf_seq, self.n_jobs, sf_algnmt)
####
####
    # parse out the qualified and save the id, position, length
    # output the "unmapped" sequences, which will be realigned using "minimap2 splice" to the reference
    def parse_MAST2_algnmt(self, sf_algnmt, sf_ref, f_cutoff, sf_unmapped):
        m_hit_mast2 = {}
        with open(sf_unmapped, "w") as fout_unmapped:
            # Here we request unique alignment (as only the sequence of gene MAST2)
            try:
                samfile = pysam.AlignmentFile(sf_algnmt, "rb", reference_filename=sf_ref)  # read in the sam file
                for algnmt in samfile.fetch(until_eof=True):  # check each alignment, and find "left" and "right" flank
                    query_name = algnmt.query_name
                    query_seq = algnmt.query_sequence
                    if query_seq is None:
                        continue
                    if algnmt.is_secondary:  # filter out secondary, but keep the supplimentary
                        continue
                    if algnmt.is_unmapped == True:  ##unmapped
                        fout_unmapped.write(">" + query_name + "\n" + query_seq + "\n")
                        continue
                    if algnmt.is_duplicate == True:  ##duplciate
                        continue
                    if algnmt.mapping_quality < global_values.TRANSDCT_UNIQ_MAPQ:  # require unique map, by default 50
                        fout_unmapped.write(">" + query_name + "\n" + query_seq + "\n")
                        continue

                    # get the map location, and length
                    query_len = len(query_seq)
                    map_pos = algnmt.reference_start  #map position
                    l_cigar=algnmt.cigar
                    i_lclip_len = 0
                    i_rclip_len = 0
                    if l_cigar[0][0] == 4:  # left clipped
                        i_lclip_len = l_cigar[0][1]
                    if l_cigar[-1][0] == 4:  # right clipped
                        i_rclip_len = l_cigar[-1][1]
                    if l_cigar[0][0] == 5:
                        i_lclip_len = l_cigar[0][1]
                        query_len += l_cigar[0][1]
                    if l_cigar[-1][0] == 5:
                        i_rclip_len = l_cigar[-1][1]
                        query_len += l_cigar[-1][1]

                    i_map_len = query_len - i_rclip_len - i_lclip_len
                    if (float(i_map_len) / float(query_len)) < f_cutoff:
                        fout_unmapped.write(">" + query_name + "\n" + query_seq + "\n")
                        continue
                    m_hit_mast2[query_name] = (map_pos, i_map_len)
            except ValueError:#
                print(sf_algnmt, "is empty")
        return m_hit_mast2
####
####
    def algn_TD_mid_contig_to_ref(self, sf_ref, sf_contig, sf_algnmt):
        xtea_contig = XTEContig(self.sf_wfolder, self.n_jobs)
        xtea_contig.align_short_contigs_with_secondary_supplementary2(sf_ref, sf_contig, self.n_jobs, sf_algnmt)

    #major goals:
    # 1) find where the transduction regions are mapped to;
    # 2) New "fusion" cases, where the middle part aligned to?
    # 3) Issue here: if only request
    def parse_second_algnmt_to_ref(self, sf_algnmt, sf_ref, f_cutoff):
        m_hit_ref = {}
        # Here we request unique alignment (as only the sequence of ref)
        try:
            samfile = pysam.AlignmentFile(sf_algnmt, "rb", reference_filename=sf_ref)  # read in the sam file
            for algnmt in samfile.fetch():  # check each alignment, and find "left" and "right" flank
                query_name = algnmt.query_name
                query_seq = algnmt.query_sequence
                if algnmt.is_secondary:  # filter out secondary, but keep the supplimentary
                    continue
                if algnmt.is_unmapped == True:  ##unmapped
                    continue
                if algnmt.is_duplicate == True:  ##duplciate
                    continue
                if algnmt.mapping_quality < self.TRANSDCT_UNIQ_MAPQ:  # require unique map, by default 50
                    continue

                # get the map location, and length
                chrm=algnmt.reference_name #which chrom it aligns to
                query_len = len(query_seq)
                map_pos = algnmt.reference_start  #map position
                l_cigar=algnmt.cigar
                i_lclip_len = 0
                i_rclip_len = 0
                if l_cigar[0][0] == 4:  # left clipped
                    i_lclip_len = l_cigar[0][1]
                if l_cigar[-1][0] == 4:  # right clipped
                    i_rclip_len = l_cigar[-1][1]
                if l_cigar[0][0] == 5:
                    i_lclip_len = l_cigar[0][1]
                    query_len += l_cigar[0][1]
                if l_cigar[-1][0] == 5:
                    i_rclip_len = l_cigar[-1][1]
                    query_len += l_cigar[-1][1]

                i_map_len = query_len - i_rclip_len - i_lclip_len
                if (float(i_map_len) / float(query_len)) < f_cutoff:
                    continue
                m_hit_ref[query_name] = (chrm, map_pos, i_map_len)
        except ValueError:
            print(sf_algnmt, "is empty")
        return m_hit_ref

####

    ####
    def define_major_type(self, s_5TD_type, s_pre_Alu, s_mid_sgmt, s_sub_family, s_suf_Alu):
        s_type=s_sub_family
        if (self.ALU in s_pre_Alu) and (self.ALU in s_suf_Alu):# and (self.MAST2 in s_mid_sgmt)
            s_type="CH10_"+s_sub_family
        elif (self.MAST2_only in s_5TD_type) and (self.ALU in s_suf_Alu):
            s_type="CH10_"+s_sub_family
        elif (self.MAST2_only in s_5TD_type) and ("SVA_F" in s_sub_family) \
                and (self.ALU not in s_pre_Alu) and (self.ALU not in s_suf_Alu):
            s_type="SVA_F1"
        elif (self.ALU in s_pre_Alu) or (self.ALU in s_suf_Alu):
            s_type = "Other_" + s_sub_family
        return s_type
####
####
    #1. re-classify insertion based on:
    ##5'/3' transduction sequence, and the middle sequence (algned to MAST2 or other regions)
    #m_cns_pos save the positions the breakpoint fall in (Alu-like, VNTR, SINE-R region, or so)
    def re_classify_INS_type(self, m_classified, m_cns_pos, m_hit_mast2, m_hit_ref, sf_new_type):
        m_final_list={}
        #with the MAST2 output, re-classify the insertion
        with open(sf_new_type, "w") as fout_new:
            fout_new.write("Id,Type,Full_copy,TD5,Pre_Alu,MAST2,Hex,Alu_like,VNTR,SINE_R,Suf_Alu,TD3,Orientation,Total_len\n")
            for s_id in m_classified:
                (s_full, i_5_unmask_len, s_pre_Alu, s_sub_family, s_suf_Alu, i_3_unmask_len,
                 i_masked_lenth, b_sva_rc) = m_classified[s_id]
                (i_start_pre_Alu, i_end_pre_Alu, i_start_suf_Alu, i_end_suf_Alu, i_pos_hex,
                 i_pos_Alu_like, i_pos_vntr, i_pos_SINE_R, n_hex, n_vntr) = m_cns_pos[s_id]#
                s_tmp_5 = s_id +  self.TD5
                s_tmp_MAST2 = s_id + self.MAST2

                i_start_MAST2=-1
                i_end_MAST2=-1
                s_5TD_type=self.NULL
                if s_tmp_5 in m_hit_mast2:#5TD hit MAST2
                    (tmp_start, tmp_len) = m_hit_mast2[s_tmp_5]
                    s_5TD_type=self.MAST2_only
                    i_start_MAST2 = tmp_start
                    i_end_MAST2 = tmp_start + tmp_len

                s_mid_sgmt=self.NULL
                if s_tmp_MAST2 in m_hit_mast2:
                    s_mid_sgmt=self.MAST2_only
                    (tmp_start, tmp_len)=m_hit_mast2[s_tmp_MAST2]
                    i_start_MAST2=tmp_start
                    i_end_MAST2=tmp_start+tmp_len

                #print s_id, s_mid_sgmt, s_sub_family, s_pre_Alu, s_suf_Alu

                TD5_src="unknown_src_"+str(i_5_unmask_len) #or NULL
                if i_5_unmask_len<self.MIN_TD_LEN:#by default 35
                    TD5_src=self.NULL
                elif s_tmp_5 in m_hit_mast2:
                    TD5_src = self.MAST2_only+"_"+str(m_hit_mast2[s_tmp_5][1])
                elif s_id+self.TD5 in m_hit_ref:
                    (src_chrm,src_pos,src_hit_len)=m_hit_ref[s_id+self.TD5]
                    TD5_src=src_chrm+"_"+str(src_pos)+"_"+str(src_hit_len)
                TD3_src="unknown_src_"+str(i_3_unmask_len)
                if i_3_unmask_len<self.MIN_TD_LEN:
                    TD3_src=self.NULL
                elif s_id+self.TD3 in m_hit_ref:
                    (src_chrm, src_pos, src_hit_len) = m_hit_ref[s_id + self.TD3]
                    TD3_src = src_chrm + "_" + str(src_pos) + "_" + str(src_hit_len)

                ####define major types
                s_major_type = self.define_major_type(s_5TD_type, s_pre_Alu, s_mid_sgmt, s_sub_family, s_suf_Alu)

                s_full_copy = self.TRUNCATED
                if ((n_hex > 0) and (n_vntr > 300)):  # or (i_end_pre_Alu>0) or (s_5TD_type != self.NULL):####
                    s_full_copy = self.FULL_COPY  #
                elif ((self.ALU in s_pre_Alu) and (self.ALU in s_suf_Alu) and (n_vntr > 200)):
                    s_full_copy = self.FULL_COPY  #
                elif (("SVA_F1" in s_major_type) and (n_vntr > 200)):
                    s_full_copy = self.FULL_COPY  #

                pre_Alu=s_pre_Alu+"_"+str(i_start_pre_Alu)+"_"+str(i_end_pre_Alu)
                s_mast2=str(i_start_MAST2)+"_"+str(i_end_MAST2)
                s_hex=str(n_hex)
                s_Alu_like=str(i_pos_Alu_like)
                s_VNTR=str(n_vntr)
                s_SINE_R=str(i_pos_SINE_R)
                suf_Alu=s_suf_Alu+"_"+str(i_start_suf_Alu)+"_"+str(i_end_suf_Alu)
                s_sva_rc="+"
                if b_sva_rc==True:
                    s_sva_rc="-"
                #s_polyA=str(20)+s_sva_rc#for temporary
                #s_id_fields=s_id.split("_")

                sinfo = s_id + "," + s_major_type + "," + s_full_copy+ ","+ TD5_src + "," + pre_Alu + "," + s_mast2 \
                        +","+s_hex + "," + s_Alu_like + "," + s_VNTR +","+s_SINE_R+","+suf_Alu+","+TD3_src+\
                        ","+s_sva_rc+","+str(i_masked_lenth)
                if b_sva_rc==True:
                    sinfo = s_id + "," + s_major_type + "," + s_full_copy + "," + TD5_src + "," + suf_Alu + "," + s_mast2 \
                            + "," + s_hex + "," + s_Alu_like + "," + s_VNTR + "," + s_SINE_R + "," +  pre_Alu + "," + TD3_src + \
                            "," + s_sva_rc + "," + str(i_masked_lenth)
                fout_new.write(sinfo+"\n")
                m_final_list[s_id]=1
        return m_final_list


    def _gnrt_seq_rc(self, s):
        s_rc=""
        for c in reversed(s):
            if c=="A" or c=="a":
                s_rc+="T"
            elif c=="C" or c=="c":
                s_rc+="G"
            elif c == "G" or c == "g":
                s_rc += "C"
            elif c == "T" or c == "t":
                s_rc += "A"
        return s_rc

####
    def _gnrt_reverse_complementary(self, s):
        if len(s)<1:
            return ""
        s_hex_rc = self._gnrt_seq_rc(s)
        s_hex_rc2="("+s_hex_rc+")n"
        return s_hex_rc2
####
    def _rotate_seq(self, s_seq):
        m_new={}
        i_len=len(s_seq)
        s_template=s_seq+s_seq
        for i in range(i_len):
            s_tmp=s_template[i:i+i_len]
            m_new[s_tmp]=1
        return m_new

####
    def load_in_hexamer(self, sf_hex):
        with open(sf_hex) as fin_hex:
            for line in fin_hex:
                s_hex_ori=line.rstrip() #in format: (ccctct)n
                s_hex=s_hex_ori[1:-2]
                self.m_hexamer[s_hex_ori]=1
                m_rotated=self._rotate_seq(s_hex)
                for s_rotate_hex in m_rotated:
                    s_rotate_hex2="("+s_rotate_hex+")n"
                    self.m_hexamer[s_rotate_hex2]=1
                    s_hex_rc=self._gnrt_reverse_complementary(s_rotate_hex)
                    self.m_hexamer_rc[s_hex_rc]=1

    def is_SVA(self, l_segmts):
        b_has_SVA=False
        b_has_Alu=False
        for rcd in l_segmts:
            segmt_id=rcd[2]
            if segmt_id in self.m_vntr:
                b_has_SVA=True
            elif self.ALU in segmt_id:
                b_has_Alu=True
        return b_has_SVA, b_has_Alu
####
####
    ####
    def is_ambigous_sub_type(self, l_segmts):
        '''
        :param l_segmts: each record in format (start_pos, end_pos, sub_type, csn_start, csn_end)
        :return:
        '''
        b_ambigous=False
        m_tmp={}
        i_lth=0
        sva_start=100000000
        sva_end=0
        for rcd in l_segmts:
            segmt_id=rcd[2]
            seq_start=int(rcd[0])
            seq_end=int(rcd[1])
            cns_end=int(rcd[4])
            if (segmt_id in self.m_vntr) and (cns_end>self.i_hexamer_lenth):
                m_tmp[segmt_id]=1
                if sva_start>seq_start:
                    sva_start=seq_start
                if sva_end<seq_end:
                    sva_end=seq_end
        if len(m_tmp)>1:
            b_ambigous=True
        i_lth=sva_end-sva_start
        if i_lth<0:
            i_lth=0
        return b_ambigous, i_lth
####

####
    #this is assume the sequence is always in positive orientation
    def classify_cannonical_SVA(self, l_segmts, m_hexamer, b_sva_rc):
        '''
        Function: heximer, SVA, polyA
        :param l_segmts: each in format (start_pos, end_pos, sub_type, csn_start, csn_end)
        :return:
        '''
        n_hex = 0  # hexamer length
        n_vntr = 0  # vntr length
        # one of the following will be recorded
        i_pos_hex = -1
        i_pos_Alu_like = -1# alu-like region
        i_pos_vntr = -1# this is approximate
        i_pos_SINE_R = -1

        sva_start = 100000000
        sva_end = 0
        sva_cns_start = 100000000
        sva_cns_end = 0
        b_hit_hexamer=False
        s_sub_family=""
        idx_extra_hexamer = -1  #
        idx_SVA = -1 #if there are several SVA records, and this is rc, then save the 1st one, otherwise, save the last
        idx_SVA_hex=-1
        idx_tmp = 0
        for rcd in l_segmts:#
            (start_pos, end_pos, sub_type, cns_start, cns_end, b_rc, i_seq_dist_from_end) = rcd
            segmt_id = sub_type
            if segmt_id in m_hexamer:
                b_hit_hexamer=True
                n_hex+=(end_pos-start_pos)
                i_pos_hex = 1
                idx_extra_hexamer=idx_tmp
            seq_start = int(start_pos)
            seq_end = int(end_pos)
            if segmt_id in self.m_vntr:#for those can be masked as "SVA" regions
                if idx_SVA==-1:
                    idx_SVA=idx_tmp
                    idx_SVA_hex=idx_tmp##
                else:
                    if b_sva_rc==False:#not rc, then always save the last one
                        idx_SVA=idx_tmp
                    else:
                        idx_SVA_hex=idx_tmp
                s_sub_family = segmt_id  # sub-family
                if cns_start<=self.i_hexamer_lenth:
                    b_hit_hexamer=True

                if sva_start > seq_start:
                    sva_start = seq_start
                if sva_end < seq_end:
                    sva_end = seq_end
                if sva_cns_start > cns_start:
                    sva_cns_start=cns_start
                if sva_cns_end < cns_end:
                    sva_cns_end=cns_end
            idx_tmp += 1

#!!!!!Bug for SVA2!!!!!!
        #this is just an VNTR expansion
        if sva_cns_end < self.vntr_end and sva_cns_start>self.vntr_start:
            return None,None

        (hexamer_len, vntr_start, vntr_end, cns_polyA_len, cns_len) = self.m_vntr[s_sub_family]
        if sva_cns_start < hexamer_len:
            if i_pos_hex<0:
                i_pos_hex = sva_cns_start
            n_hex += (hexamer_len - sva_cns_start)
        elif sva_cns_start >= hexamer_len and sva_cns_start < vntr_start:  # hit the Alu-like region
            i_pos_Alu_like = sva_cns_start
        elif sva_cns_start >= vntr_start and sva_cns_start <= vntr_end:  # vntr region
            i_pos_vntr = sva_cns_start - vntr_start
        elif sva_cns_start > vntr_end:
            i_pos_SINE_R = sva_cns_start - vntr_end

        if (sva_cns_start > vntr_end) or (sva_cns_end<vntr_start):#
            n_vntr = 0
        else:
            n_vntr = (sva_end - sva_start) - (sva_cns_end - sva_cns_start) + (vntr_end - vntr_start) - n_hex

        #i_5_unmask_len = self.get_5_transduct_len(l_segmts)
        #i_3_unmask_len, n_extra_polyA = self.get_3_transduct_len(l_segmts)
        i_masked_lenth = sva_end-sva_start
        idx_hex=idx_extra_hexamer
        if idx_hex==-1:
            idx_hex=idx_SVA_hex
        seq_rcd=(idx_hex, self.NOT_ALU, self.mid_segmt_start, self.mid_segmt_end, s_sub_family, b_hit_hexamer,
                i_masked_lenth, self.NOT_ALU , idx_SVA, b_sva_rc)
        pos_rcd = (-1, -1, -1, -1, i_pos_hex,
                   i_pos_Alu_like, i_pos_vntr, i_pos_SINE_R, n_hex, n_vntr)

        #save the SINE-R record
        sine_r_rcd=()
        return seq_rcd, pos_rcd
####
####
####
    #when parse the sequence, we first reverse-complementary the seq
    def parse_SINE_R_region(self, l_segmts, b_rc=False):
        sine_r_start = -1
        sine_r_end = -1

        for rcd in l_segmts:
            (start_pos, end_pos, sub_type, cns_start, cns_end, b_rc_rcd, i_seq_dist_from_end) = rcd
            segmt_id = sub_type

            if segmt_id in self.m_vntr:  # for those can be masked as "SVA" regions
                i_vntr_end=self.m_vntr[segmt_id][2]
                i_sva_len=self.m_vntr[segmt_id][4]
                if cns_end<i_vntr_end:
                    continue
                tail_offset=i_sva_len-cns_end
                sine_r_cns=i_sva_len-i_vntr_end
                seq_sine_r=sine_r_cns-tail_offset

                if seq_sine_r > (sine_r_end-sine_r_start):
                    if b_rc==False:
                        sine_r_end=end_pos
                        sine_r_start=sine_r_end-seq_sine_r
                    else:
                        sine_r_start = start_pos
                        sine_r_end = start_pos + seq_sine_r
        return sine_r_start, sine_r_end

####
    def is_SVA_reverse_complementary(self, l_segmts):
        m_ori={}
        m_ori[0]=0 #count number of records that are NOT reverse complementary
        m_ori[1]=0
        for rcd in l_segmts:
            (start_pos, end_pos, sub_type, cns_start, cns_end, b_rc, i_seq_dist_from_end) = rcd
            # double check masked hexamer is hexamer
            if "SVA" in sub_type:
                if b_rc==True:
                    m_ori[1]+=1
                else:
                    m_ori[0] += 1
        if m_ori[0]>m_ori[1]:
            return False
        return True
#
####
    #this is the subfamily containing Alu.
    # Cannonical structure: 5TD+Alu+MAST2+SVA+Alu+3TD+polyA
    def classify_CH10_subfamily_SVA(self, l_segmts, b_sva_rc):
        i_start_pre_Alu=-1#start position on consensus
        i_end_pre_Alu = -1#end position on consensus
        i_start_suf_Alu=-1
        i_end_suf_Alu = -1

        n_hex = 0  # hexamer length
        n_vntr = 0  # vntr length
        #one of the following will be recorded
        i_pos_hex=-1#hit position on the hexamer region
        i_pos_Alu_like=-1#alu-like region
        i_pos_vntr=-1 #this is approximate
        i_pos_SINE_R = -1

        NO_WHERE=-10000
        idx_SVA=NO_WHERE
        idx_SVA_end=NO_WHERE
        sva_start = 100000000
        sva_end = 0
        sva_cns_start=100000000
        sva_cns_end=0
        #sva_cns_dist_from_end=100000000
        b_hit_hexamer = False
        tmp_idx=0
        s_sub_family=""

        for rcd in l_segmts:
            (start_pos, end_pos, sub_type, cns_start, cns_end, b_rc, i_seq_dist_from_end)=rcd#
#double check masked hexamer is hexamer
            segmt_id = rcd[2]
            if segmt_id in self.m_hexamer:
                n_hex+=(end_pos-start_pos)
                b_hit_hexamer = True#

            seq_start = int(start_pos)
            seq_end = int(end_pos)
            if segmt_id in self.m_vntr:# for those can be masked as "SVA" regions
                s_sub_family = segmt_id
                if idx_SVA==NO_WHERE:#
                    idx_SVA=tmp_idx#save the first SVA record start index
                idx_SVA_end=tmp_idx#save the last SVA record end index
                #cns_start = rcd[3]
                if cns_start <= self.i_hexamer_lenth:
                    b_hit_hexamer = True

                if sva_start > seq_start:
                    sva_start = seq_start
                if sva_end < seq_end:
                    sva_end = seq_end
                if sva_cns_start > cns_start:
                    sva_cns_start=cns_start
                if sva_cns_end < cns_end:
                    sva_cns_end=cns_end
            tmp_idx+=1
        if idx_SVA==NO_WHERE:
            return None, None
        if sva_cns_end<self.vntr_start+60:#No VNTR region, it possible that same segment is masked twice as SVA and Alu.
            return None, None
####
        (hexamer_len, vntr_start, vntr_end, cns_polyA_len, cns_len)=self.m_vntr[s_sub_family]
        if sva_cns_start<hexamer_len:
            i_pos_hex=sva_cns_start
            n_hex+=(hexamer_len-sva_cns_start)
        elif sva_cns_start>=hexamer_len and sva_cns_start<vntr_start:#hit the Alu-like region
            i_pos_Alu_like=sva_cns_start
        elif sva_cns_start>=vntr_start and sva_cns_start<=vntr_end:#vntr region
            i_pos_vntr=sva_cns_start-vntr_start
        elif sva_cns_start>vntr_end:
            i_pos_SINE_R=sva_cns_start-vntr_end

        if sva_cns_start>vntr_end:
            n_vntr=0
        else:
            n_vntr=(sva_end-sva_start)-(sva_cns_end-sva_cns_start)+(vntr_end-vntr_start)-n_hex
####
####
        #check the records before SVA records, find the last one if more than one
        s_pre_Alu=self.NOT_ALU
        i_pre_Alu_end=-1#pre-Alu end position (on the seq)
        idx_pre_Alu=-1
        idx_tmp=0
        for rcd in l_segmts[:idx_SVA]:#save the Alu subfamily
            # (start_pos, end_pos, sub_type, cns_start, cns_end, b_rc, i_seq_dist_from_end)=rcd
            segmt_id = rcd[2] #
            if self.ALU in segmt_id:
                s_pre_Alu=segmt_id
                i_pre_Alu_end = int(rcd[1])
                i_start_pre_Alu=int(rcd[3])
                i_end_pre_Alu=int(rcd[4])
                idx_pre_Alu=idx_tmp
            idx_tmp+=1

        #check the records after SVA records, find the first one if more than one
        idx_suf_Alu=-1
        idx_tmp=idx_SVA_end+1#
        s_suf_Alu=self.NOT_ALU
        for rcd in l_segmts[idx_SVA_end+1:]:
            segmt_id = rcd[2]
            if self.ALU in segmt_id:
                s_suf_Alu=segmt_id
                i_start_suf_Alu=int(rcd[3])
                i_end_suf_Alu = int(rcd[4])
                idx_suf_Alu=idx_tmp
                break
            idx_tmp+=1
        ####
        idx_5_Alu=idx_pre_Alu
        idx_3_Alu=idx_suf_Alu
        i_MAST2_start = i_pre_Alu_end
        i_MAST2_end = sva_start
        b_sva_rc=self.is_SVA_reverse_complementary(l_segmts)
        if b_sva_rc==True:
            i_MAST2_start = sva_end
            i_MAST2_end = i_start_suf_Alu

#check the transduction !!!!
####
        #i_5_unmask_len = self.get_5_transduct_len(l_segmts)
        #i_3_unmask_len, n_extra_polyA = self.get_3_transduct_len(l_segmts)
        i_masked_lenth = sva_end - sva_start #only the SVA part

####
        idx_start=idx_5_Alu
        if idx_start<0:
            idx_start=idx_SVA
        idx_end=idx_3_Alu
        if idx_end<0:
            idx_end=idx_SVA_end
        if b_sva_rc==True:
            idx_start=idx_3_Alu
            if idx_start<0:
                idx_start=idx_SVA_end
            idx_end=idx_5_Alu
            if idx_end<0:
                idx_end=idx_SVA

        seq_rcd=(idx_start, s_pre_Alu, i_MAST2_start, i_MAST2_end, s_sub_family, b_hit_hexamer, i_masked_lenth,
                s_suf_Alu, idx_end, b_sva_rc)
        #template_rcd=()#save the position on the template
        pos_rcd=(i_start_pre_Alu, i_end_pre_Alu, i_start_suf_Alu, i_end_suf_Alu, i_pos_hex,
                 i_pos_Alu_like, i_pos_vntr, i_pos_SINE_R, n_hex, n_vntr)
        return seq_rcd, pos_rcd

    ####
    def parse_TD_lenth_seq(self, l_segmts, s_contig, idx_hex, idx_sva, b_sva_rc):
        i_5_unmask_len = 0
        i_3_unmask_len = 0
        i_extra_polyA = 0
        s_5_seq=""
        s_3_seq=""

        i_seq_len=len(s_contig)
        if b_sva_rc==True:#reverse complementary
            if idx_sva >=0:#check the 3' side (at the tail)
                idx_tmp=idx_sva
                if idx_sva>0:
                    s_sgmt_id = l_segmts[idx_sva - 1][2]
                    if ("(T)" in s_sgmt_id) or ("(TT" in s_sgmt_id) or ("TT)" in s_sgmt_id):
                        i_extra_polyA = l_segmts[idx_sva - 1][1] - l_segmts[idx_sva - 1][0]
                        idx_tmp-=1
                i_3_start=l_segmts[idx_tmp][0]
                s_3_seq=s_contig[:i_3_start]
                i_3_unmask_len=i_3_start

            if idx_hex>=0:#check the 5' side (at the front)
                i_5_end=l_segmts[idx_hex][1]
                s_5_seq=s_contig[i_5_end:]
                i_5_unmask_len=i_seq_len-i_5_end

        else:#non reverse complementary
            if idx_hex>=0:#check the 5' side (at the front)
                i_5_end=l_segmts[idx_hex][0]
                s_5_seq=s_contig[:i_5_end]
                i_5_unmask_len=i_5_end
            if idx_sva>=0:#check the 3' side (at the tail)
                idx_tmp=idx_sva
                if idx_sva<(len(l_segmts)-1):
                    s_sgmt_id=l_segmts[idx_sva+1][2]
                    if ("(A)" in s_sgmt_id) or ("(AA" in s_sgmt_id) or ("AA)" in s_sgmt_id):
                        #SVA extra polyA
                        #if abs(l_segmts[idx_sva][1]-l_segmts[idx_sva+1][0])<100:
                        i_extra_polyA=l_segmts[idx_sva+1][1]-l_segmts[idx_sva+1][0]
                        idx_tmp+=1

                i_3_start=l_segmts[idx_tmp][1]
                s_3_seq=s_contig[i_3_start:]
                i_3_unmask_len=i_seq_len-i_3_start
        return i_5_unmask_len, i_3_unmask_len, i_extra_polyA, s_5_seq, s_3_seq

####
    def classify_masked_record(self, s_id, l_segmts):#
        t_ins_sine_r = None
        i_sva_size=0
        b_has_sva, b_has_alu=self.is_SVA(l_segmts)
        if b_has_sva==False:
            return None, None, None, i_sva_size, t_ins_sine_r
        b_ambigous_type, i_sva_size=self.is_ambigous_sub_type(l_segmts)
        if b_ambigous_type==True:
            print(s_id, l_segmts)
            return None, None, self.TYPE_AMBIGOUS, i_sva_size, t_ins_sine_r

        rslt_rcd=None
        pos_rcd=None
        b_sva_rc = self.is_SVA_reverse_complementary(l_segmts)

        if b_has_alu==True:
            #this is the CH10 Group (or novel)
            rslt_rcd, pos_rcd=self.classify_CH10_subfamily_SVA(l_segmts, b_sva_rc)
        else:#this is cannonical type
####add b_sva_rc
            m_tmp_hexamer=self.m_hexamer
            if b_sva_rc==True:
                m_tmp_hexamer=self.m_hexamer_rc
            rslt_rcd, pos_rcd=self.classify_cannonical_SVA(l_segmts, m_tmp_hexamer, b_sva_rc)

            sine_r_start, sine_r_end=self.parse_SINE_R_region(l_segmts)
            t_ins_sine_r=(s_id, sine_r_start, sine_r_end)
        return rslt_rcd, pos_rcd, None, i_sva_size, t_ins_sine_r
####
####
    def load_in_contigs(self, sf_contig):
        m_contigs={}
        with pysam.FastxFile(sf_contig) as fh:
            for entry in fh:
                m_contigs[str(entry.name)]=str(entry.sequence)
        return m_contigs

####

####
class TEClassifier():
    def __init__(self, sf_wfolder, n_jobs):
        self.sf_wfolder = sf_wfolder
        self.n_jobs = n_jobs
        self.TE_type="None"
        self.SUPER_FAMILY="None"
        self.MIN_CNS_END=0
        self.SEPARATOR="_"
        self.rmsker="RepeatMasker"

    def set_separator(self, s_sep):#
        self.SEPARATOR=s_sep

    def set_TE_type(self, s_type):
        self.TE_type=s_type
        self.SUPER_FAMILY=s_type

    def set_TE_super_family(self, s_type):
        self.SUPER_FAMILY=s_type

    def set_min_cns_end(self, cns_end):
        self.MIN_CNS_END=cns_end

####
    def classify_from_rmsk_output(self, sf_annotation, m_contigs, m_already_called, f_hit_cutoff, i_min_len,
                                  i_max_len, sf_out, b_sample_id, s_pre_define_type=None):
        m_records = {}  # save the records
        with open(sf_annotation) as fin_rmsk:
            for line in fin_rmsk:
                fields = line.split()
                if len(fields) < 4:
                    continue
                if ("SW" in fields[0]) or ("score" in fields[0]):
                    continue

                s_id = fields[4]  # in format: HG02106_chr6_24648990
                if s_id in m_already_called==True:#skip already exist
                    continue
                tmp_fields = s_id.split(self.SEPARATOR)
####change to 4 temporarily
                if len(tmp_fields) > 3:  # skip those duplicated ones (like HG02106_chr6_24648990_1)
                    continue
                start_pos = int(fields[5])
                end_pos = int(fields[6])
                i_seq_dist_from_end = int(fields[7][1:-1])
                b_rc = False
                if fields[8] == "C":
                    b_rc = True
                sub_type = fields[9]
                if (s_pre_define_type is not None) and (s_pre_define_type not in sub_type):
                    continue
                s_tye=fields[10]
                csn_start = -1
                csn_end = -1
                if b_rc == True:
                    csn_start = int(fields[13])
                    csn_end = int(fields[12])
                else:
                    csn_start = int(fields[11])
                    csn_end = int(fields[12])
                if s_id not in m_records:
                    m_records[s_id] = []
                m_records[s_id].append((start_pos, end_pos, sub_type, s_tye, csn_start, csn_end, i_seq_dist_from_end))
                #print s_id, start_pos, end_pos, sub_type, csn_start, csn_end, i_seq_dist_from_end, "test_label"

        #print len(m_records), "number of records checked for Alu!"
        #print m_records
        m_final_list=self.classify_records(m_records, m_contigs, f_hit_cutoff, i_min_len, i_max_len, sf_out, b_sample_id)
        return m_final_list
####

####
    def classify_records(self, m_records, m_contigs, f_cutoff, i_min_len, i_max_len, sf_rslts, b_sample_id=False):
        m_final_list={}
        m_sites={}
        with open(sf_rslts, "w") as fout:
            for s_id in m_records:#
                if s_id not in m_contigs:
                    print(("%s not in m_contigs\n" % s_id))
                    continue
                s_contig_seq = m_contigs[s_id]  # contig seq
                n_seq_len=len(s_contig_seq)
                if n_seq_len<i_min_len:#require the minimal length
                    continue
                if n_seq_len>i_max_len:
                    continue
                l_segmts = m_records[s_id]
                #print s_id, l_segmts
                b_qualified, rcd=self.is_qualified_record(l_segmts, n_seq_len, f_cutoff)
                if b_qualified==False:
                    #print rcd, "not qualified"
                    continue
                (s_sub_family, cns_start, cns_end, seq_start, seq_end) = rcd

                chrm = None
                pos=None
                if b_sample_id==False: #in format: chrm~pos
                    s_id_fields=s_id.split(self.SEPARATOR)
                    if len(s_id_fields)<2:
                        continue
                    chrm=s_id_fields[0]
                    pos=int(s_id_fields[1])
                else:#in format: HG02106_chr6_24648990
                    s_id_fields = s_id.split(self.SEPARATOR)
                    if len(s_id_fields) < 3:
                        continue
                    chrm = s_id_fields[1]
                    pos = int(s_id_fields[2])
                if chrm not in m_sites:
                    m_sites[chrm]={}
                i_start=pos-global_values.NEARBY_REGION
                i_end=pos+global_values.NEARBY_REGION
                b_exist=False
                for i in range(i_start, i_end):
                    if i in m_sites[chrm]:
                        b_exist=True
                        break
                if b_exist==True:
                    continue
                if int(cns_end)<self.MIN_CNS_END:
                    continue

                m_sites[chrm][pos]=1
                #in vcf format:
                fout.write(chrm+"\t"+str(pos)+"\t.\t.\t"+s_contig_seq+
                           "\t60\t.\tSUB_FAMILY="+s_sub_family+";CNS="+str(cns_start)+"-"+str(cns_end)+
                           ";SEQ_REGION="+str(seq_start)+"-"+str(seq_end)+"\n")
                # fout.write(s_id_fields[0]+"\t"+s_id_fields[1]+"\t"+s_sub_family+"\t"+
                #            str(cns_start)+"\t"+str(cns_end)+"\t"+str(seq_start)+"\t"+str(seq_end)+"\n")
                m_final_list[s_id]=(chrm, pos, s_sub_family, cns_start, cns_end, seq_start, seq_end, s_contig_seq)
        return m_final_list
####
####
    def run_RMSK(self, sf_fa):
        cmd_runer=CMD_RUNNER()
        cmd="%s %s" % (self.rmsker, sf_fa)
        cmd_runer.run_cmd_small_output(cmd)

####
    def is_qualified_record(self, l_segmts, i_seq_len, f_cutoff):
        rcd=None
        if self.is_specific_TE(l_segmts)==False:
            return False, rcd
        rcd, n_masked=self.classify_cannonical_TE_ins(l_segmts)
        #(s_sub_family, cns_start, cns_end, seq_start, seq_end)=rcd
        n_extra_polyA=self.get_extra_polyA_len(l_segmts)
        n_cnt=n_masked + n_extra_polyA

        if float(n_cnt)/float(i_seq_len)>f_cutoff:
            return True, rcd
        return False, rcd

####
    def is_specific_TE(self, l_segmts):
        b_is_TE = False
        for rcd in l_segmts:
            segmt_id = rcd[2]
            if self.TE_type in segmt_id:
                b_is_TE = True
        return b_is_TE

    def classify_cannonical_TE_ins(self, l_segmts):
        start = 100000000
        end = 0
        cns_start = 100000000
        cns_end = 0
        s_sub_family = ""
        n_masked=0

        for rcd in l_segmts:#
            (start_pos, end_pos, sub_type, s_type, cns_start_tmp, cns_end_tmp, i_seq_dist_from_end) = rcd
            segmt_id = sub_type
            seq_start = int(start_pos)
            seq_end = int(end_pos)
            if (self.TE_type in segmt_id) and (self.SUPER_FAMILY in s_type):  # for those can be masked as "SVA" regions
                s_sub_family = segmt_id  # sub-family
                if start > seq_start:
                    start = seq_start
                if end < seq_end:
                    end = seq_end
                if cns_start > cns_start_tmp:
                    cns_start = cns_start_tmp
                if cns_end < cns_end_tmp:
                    cns_end = cns_end_tmp
                n_masked+=(seq_end-seq_start)
        #print s_sub_family, cns_start, cns_end, start, end, n_masked
        return (s_sub_family, cns_start, cns_end, start, end), n_masked
####
####
    def get_extra_polyA_len(self, l_segmts):
        s_id=l_segmts[-1]
        n_extra_polyA=0
        if ("(A" in s_id) or ("A-rich" in s_id):  # last one is polyA
            # find the gap between previous and current one
            n_extra_polyA = l_segmts[-1][1] - l_segmts[-1][0]
        return n_extra_polyA

####
class AluClassifier(TEClassifier):
    def __init__(self, sf_wfolder, n_jobs):
        TEClassifier.__init__(self, sf_wfolder, n_jobs)
        self.ALU = "Alu"
        self.set_TE_type(self.ALU)

####
class LINE1Classifier(TEClassifier):
    def __init__(self, sf_wfolder, n_jobs):
        TEClassifier.__init__(self, sf_wfolder, n_jobs)
        self.LINE1 = "L1"
        self.set_TE_type(self.LINE1)
        cns_end=5500#at least 5500bp
        self.set_min_cns_end(cns_end)

class HERVClassifier(TEClassifier):
    def __init__(self, sf_wfolder, n_jobs):
        TEClassifier.__init__(self, sf_wfolder, n_jobs)
        self.HERV = "HERV"
        self.set_TE_type(self.HERV)
        self.ERV = "ERV"
        self.set_TE_super_family(self.ERV)
        cns_end=5500
        self.set_min_cns_end(cns_end)
        self.set_separator(global_values.SEPERATOR)

    # def run_rmsk(self, sf_fa):
    #     cmd="RepeatMasker %s" % sf_fa

    # def slct_hit_cns_end(self, sf_in, sf_out):
    #     with open(sf_in) as fin_rslt:
    #         for line in fin_rslt:
    #             fields=line.split()
####