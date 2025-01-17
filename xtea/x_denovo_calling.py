##03/23/2021
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

#this module is designed for calling de novo TE insertion

from x_somatic_calling import *

####
class DenovoMode(CaseControlMode):
    def __init__(self, sf_ref, s_wfolder, n_jobs):
        CaseControlMode.__init__(self, sf_ref, s_wfolder, n_jobs)
        self.s_wfolder=s_wfolder
        if os.path.exists(s_wfolder)==False:
            self.s_wfolder="./"
        elif s_wfolder[-1]!="/":
            s_wfolder+="/"
        self.sf_ref=sf_ref
        self.n_jobs=n_jobs
####
#####Here assume we have the "case" results from the proband
    #then here we check for both the contrls
    #Input: assume "sf_bam_list" in format: each case bam in one line
    def call_denovo_TE_insertions(self, sf_bam_list, sf_case_cns, sf_case_cns2, extnd, nclip_cutoff, ndisc_cutoff,
                                  npolyA_cutoff, sf_rep_cns, sf_flank, i_flk_len, bin_size, sf_out):

        # Input: candidate_disc_filtered_cns.txt and transduction events
        #First, merge the two inputs
        sf_case_candidates=sf_case_cns+".merged_case_cns_candidates_for_denovo.txt"
        if os.path.isfile(sf_case_cns2)==False:
            sf_case_candidates=sf_case_cns
        else:
            with open(sf_case_cns) as fin_cns, open(sf_case_cns2) as fin_cns2, open(sf_case_candidates,"w") as fout:
                for line in fin_cns:
                    fout.write(line.rstrip()+"\n")
                for line in fin_cns2:
                    if global_values.NOT_TRANSDUCTION not in line:
                        fout.write(line.rstrip() + "\n")

        n_ctrl=1
        l_sf_bam=[]
        with open(sf_bam_list) as fin_bam_list:#assume each line is one ctrl bam
            for line in fin_bam_list:#
                sf_bam_tmp=sf_bam_list+".ctrl_%d.txt" % n_ctrl
                l_sf_bam.append(sf_bam_tmp)
                with open(sf_bam_tmp,"w") as fout_tmp:
                    fout_tmp.write(line.rstrip())
                n_ctrl+=1

        sf_out_tmp=sf_out+".after_parent1.txt"
        if len(l_sf_bam)==1:
            sf_out_tmp = sf_out
        #check against the first parent
        self.call_somatic_TE_insertion_case_ctrl_denovo(l_sf_bam[0], sf_case_candidates, extnd, nclip_cutoff, ndisc_cutoff,
                                  npolyA_cutoff, sf_rep_cns, sf_flank, i_flk_len, bin_size, sf_out_tmp)
        #check against the second parent
        if len(l_sf_bam)>1:##
            self.call_somatic_TE_insertion_case_ctrl_denovo(sf_bam_list, sf_out_tmp, extnd, nclip_cutoff, ndisc_cutoff,
                                           npolyA_cutoff, sf_rep_cns, sf_flank, i_flk_len, bin_size, sf_out)
####
####
    ####nclip, ndisc are the cutoff when calling somatic events
    #sf_bam_list is in format: s_id sf_control1 sf_control2
    def call_somatic_TE_insertion_case_ctrl_denovo(self, sf_bam_list, sf_case_candidates, extnd, nclip_cutoff, ndisc_cutoff,
                                  npolyA_cutoff, sf_rep_cns, sf_flank, i_flk_len, bin_size, sf_out, b_tumor=False):
        #separate the transduction and non-transduction cases
        #for non-transduction cases:
        xclip_disc = XClipDiscFilter(sf_bam_list, self.s_wfolder, self.n_jobs, self.sf_ref)
        sf_non_td = sf_case_candidates + ".non_td.tmp"
        sf_td=sf_case_candidates+".td.tmp"
        sf_td_sibling=sf_case_candidates+".td_sibling.tmp"
        sf_orphan=sf_case_candidates+".orphan.tmp"
        xclip_disc.sprt_TEI_to_td_orphan_non_td(sf_case_candidates, sf_non_td, sf_td, sf_td_sibling, sf_orphan)#

        #parse, realign, cluster clip and discordant reads to consensus, and save to a file
        sf_tmp_cluster=sf_out+".candidate_somatic_cluster_from_ctrl.txt"
        m_clip_cluster_non_td=xclip_disc.call_clip_disc_cluster(sf_non_td, self.iextnd, self.bin_size, sf_rep_cns,
                                          self.bmapped_cutoff,self.i_concord_dist, self.f_concord_ratio, nclip_cutoff,
                                          ndisc_cutoff, self.s_wfolder, sf_tmp_cluster)
        #m_read_depth, m_gntp_info=self.extract_sites_raw_gntp_feature_depth(sf_non_td, sf_bam_list) ####

        #First round filter for TD cases (this is for all TDs, including )
        # parse, realign, cluster clip and discordant reads to consensus, and save to a file
        sf_tmp_td_cluster = sf_out + ".candidate_somatic_td_cluster_from_ctrl.txt"
        m_clip_cluster_td = xclip_disc.call_clip_disc_cluster(sf_td, self.iextnd, self.bin_size, sf_rep_cns,
                                                                  self.bmapped_cutoff, self.i_concord_dist,
                                                                  self.f_concord_ratio, nclip_cutoff,
                                                                  ndisc_cutoff, self.s_wfolder, sf_tmp_td_cluster)#
        #m_read_depth_td, m_gntp_info_td = self.extract_sites_raw_gntp_feature_depth(sf_td, sf_bam_list)

        #Second round filter for TD cases
        #parse, realign to flanking sequences, and save to file
        xtd=XTransduction(self.s_wfolder, self.n_jobs, self.sf_ref)
        sf_tmp_cluster_td=sf_out+".candidate_somatic_cluster_from_ctrl_td.txt"
        m_td_cluster, m_td_polyA=xtd.collect_realign_reads(sf_td, sf_case_candidates, xclip_disc, sf_flank,
                                                           i_flk_len, extnd, bin_size, ndisc_cutoff, sf_rep_cns,
                                                           sf_tmp_cluster_td)
####
        with open(sf_out+".tmp_ctrl_clip_polyA", "w") as fout_clip_polyA:
            for ins_chrm in m_td_polyA:
                for ins_pos in m_td_polyA[ins_chrm]:
                    rcd=m_td_polyA[ins_chrm][ins_pos]
                    sinfo="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t\t{7}\n".format(ins_chrm, ins_pos, rcd[0], rcd[1], rcd[2],
                                                                              rcd[3], rcd[4], rcd[5])
                    fout_clip_polyA.write(sinfo)
####
        # Comment (09-15-20): it is possible a wrongly annotated sibling, but actually a germline SV
        # we need a module to filter out these events

        #collect features for orphan cases from control bam
        #check each sites and filter out germline cases
        xorphan=XOrphanTransduction(self.s_wfolder, self.n_jobs, self.sf_ref)

        m_need_filter_out=xorphan.check_features_for_given_sites(sf_orphan, sf_bam_list, nclip_cutoff, ndisc_cutoff)

        #(this is for sibling cases, we add another checking 09-01-2020)
        m_filter_out_sibling = xorphan.check_features_for_given_sites(sf_td_sibling, sf_bam_list, nclip_cutoff,
                                                                      ndisc_cutoff)#
####
        #print m_clip_cluster
        #parse out the somatic events
        m_cluster_info=xclip_disc.load_TEI_info_from_file(sf_tmp_cluster)
        m_cluster_info_td = xclip_disc.load_TEI_info_from_file(sf_tmp_td_cluster)
        sf_before_final=sf_out+".before_final_filter"
        m_unique_tmp={}
        with open(sf_before_final, "w") as fout_rslt:
            with open(sf_non_td) as fin_sites:
                for line in fin_sites:
                    fields=line.split()
                    chrm = fields[0]
                    pos = int(fields[1])
                    if chrm not in m_unique_tmp:
                        m_unique_tmp[chrm]={}
                    if pos in m_unique_tmp[chrm]:
                        continue
                    else:
                        m_unique_tmp[chrm][pos]=1
                    if self.pass_cns_chk(m_cluster_info, m_clip_cluster_non_td, fields, npolyA_cutoff)==True:
                        fout_rslt.write(line)
            with open(sf_td) as fin_sites:
                for line in fin_sites:
                    l_fields=line.split()
                    chrm=l_fields[0]
                    pos=int(l_fields[1])
                    if chrm not in m_unique_tmp:
                        m_unique_tmp[chrm]={}
                    if pos in m_unique_tmp[chrm]:
                        continue
                    else:
                        m_unique_tmp[chrm][pos]=1
                    if self.pass_td_chk(m_td_cluster, m_td_polyA, l_fields, npolyA_cutoff)==True:
                        if self.pass_cns_chk(m_cluster_info_td, m_clip_cluster_td, l_fields, npolyA_cutoff) == True:
                            if (chrm in m_filter_out_sibling) and (pos in m_filter_out_sibling[chrm]):
                                continue
                            fout_rslt.write(line)
            with open(sf_orphan) as fin_orphan:
                for line in fin_orphan:
                    fields=line.split()
                    chrm=fields[0]
                    pos=int(fields[1])
                    if chrm not in m_unique_tmp:
                        m_unique_tmp[chrm]={}
                    if pos in m_unique_tmp[chrm]:
                        continue
                    else:
                        m_unique_tmp[chrm][pos]=1
                    if (chrm in m_need_filter_out) and (pos in m_need_filter_out[chrm]):
                        continue
                    fout_rslt.write(line)

        f_disc_ef_ratio=0.1 #minimum effective discordant pair ratio
        f_ef_disc_cutoff=0.05 #minimal AF of effective discordant pairs in case
        f_disc_cutoff= 0.15#minimal AF of discordant pairs in case
        f_case_ctrl_cutoff= 0.35#maximum ctrl-ratio/case-ratio
        f_ctrl_clip_ratio=0.75 #maximum clip ratio (higher means lots of clipped reads at the breakpoints)
        f_ctrl_raw_clip_flex_ratio=0.35 #this is a flexible ratio, if case has > ratio, and in case no polyA reads
        f_ctrl_polyA_cutoff=0.3
        f_case_raw_clip_cutoff=0.05
        self.filter_by_case_ctrl_disc_depth_ratio(sf_before_final, sf_bam_list, f_disc_ef_ratio, f_ef_disc_cutoff, f_disc_cutoff, f_case_raw_clip_cutoff,
                                                  f_case_ctrl_cutoff, f_ctrl_clip_ratio, f_ctrl_raw_clip_flex_ratio, f_ctrl_polyA_cutoff, sf_out)
####