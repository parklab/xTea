##09/05/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com

import os
import glob
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
from x_clip_disc_filter import *
from x_sites import *
from bwa_align import *
import global_values
from cmd_runner import *

def unwrap_self_algn_read(arg, **kwarg):
    return XMutation.align_reads_to_cns_one_site(*arg, **kwarg)
def unwrap_self_call_snp(arg, **kwarg):
    return XMutation.call_snp_one_site(*arg, **kwarg)

####
class XMutation():
    def __init__(self, s_working_folder):
        #note, the working_folder should be same as the reads collection step
        self.working_folder = s_working_folder
        if "/" != self.working_folder[-1]:
            self.working_folder[-1] += "/"
        if os.path.exists(self.working_folder)==False:
            os.mkdir(self.working_folder)

        self.cmd_runner = CMD_RUNNER()
        #create the mutation folder
        sf_mutation_folder=self.working_folder+global_values.MUTATION_FOLDER
        if os.path.exists(sf_mutation_folder)==False:
            cmd="mkdir {0}".format(sf_mutation_folder)
            #Popen(cmd, shell=True, stdout=PIPE).communicate()
            self.cmd_runner.run_cmd_small_output(cmd)
        ####
        self.sf_bam_list=""
        self.n_jobs=1
        self.sf_bam_ref=""
        self.sf_candidate_list=""
        self.INTERNAL="internal" #indicates this is for internmal mutations
        self.FLANKING="flank" #indicates this is for flanking mutations
        self.n_batch_vcf=500 #number of vcf files for merging per batch
        self.rep_type_L1="LINE1" #this name should match one of the sequence id in the consensus library
        self.rep_type_Alu = "AluY"
        self.rep_type_SVA="SVA_F"
####
    ####
    def set_configuration(self, sf_bam_list, sf_bam_ref, sf_candidate_list, n_jobs):
        self.sf_bam_list = sf_bam_list
        self.sf_bam_ref = sf_bam_ref
        self.sf_candidate_list = sf_candidate_list #insertion list
        self.n_jobs = n_jobs

    def _select_high_quality_sites(self, sf_sites, i_len_cutoff):
        xsites = XSites(sf_sites)
        m_sites = xsites.load_in_qualified_sites_from_xTEA_output(i_len_cutoff)
        return m_sites

    ####this is for L1 consensus, will only extract the first one
    def _prep_consensus_seq(self, sf_old_cns, i_rep_type, sf_new_cns):##
        s_consensus_rep_name=self.__get_rep_name_in_consensus_by_id(i_rep_type)
        b_find_rep=False
        s_cur_rcd=""
        with pysam.FastxFile(sf_old_cns) as fin_old, open(sf_new_cns,"w") as fout:
            for entry in fin_old:
                s_rname = str(entry.name)
                s_cur_rcd=">" + s_rname + "\n" + str(entry.sequence) + "\n"
                if s_rname==s_consensus_rep_name:
                    fout.write(s_cur_rcd)
                    b_find_rep=True
                    break
            if b_find_rep==False:#didn't find, then just use the last record
                fout.write(s_cur_rcd)

        #bwa index the file
        bwa_align = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, self.n_jobs)
        bwa_align.index_ref_file(sf_new_cns)
        bwa_align.index_fa_file(sf_new_cns)  # samtools faidx

    ####
    def __get_rep_name_in_consensus_by_id(self, i_rep_type):
        if i_rep_type==1:
            return self.rep_type_L1
        elif i_rep_type==2:
            return self.rep_type_Alu
        elif i_rep_type==4:
            return self.rep_type_SVA
        return self.rep_type_L1
####
    ####
    def run_mutation_calling_with_given_list(self, extnd, bin_size, sf_rep_cns, i_rep_type, i_dp_cutoff_internal, i_dp_cutoff_flank,
                                             sf_merged_internal_vcf, sf_merged_flanking_vcf):
        #1. first select the sites
        #m_sites=self._select_high_quality_sites(sf_sites, i_len_cutoff)

        #1.1 prepare the new consensus sequence
        sf_new_cns=self.working_folder+"consensus_seq.fa"
        self._prep_consensus_seq(sf_rep_cns, i_rep_type, sf_new_cns)

        #2. extract the reads ##
        sf_reads_wfolder=self.working_folder+"%s/"%self.INTERNAL
        if os.path.exists(sf_reads_wfolder)==False:
            os.mkdir(sf_reads_wfolder)

        #2.1 gnrt the flanking sequences ##
        xref=XReference() ##
        sf_flank_wfolder = self.working_folder + "%s/"%self.FLANKING
        if os.path.exists(sf_flank_wfolder)==False:
            os.mkdir(sf_flank_wfolder)

        #sf_flank_fa=sf_flank_wfolder+"{0}{1}{2}_flanks.fa" ##
        print("Generate the flanking sequences") ####
        xref.gnrt_target_flank_seq_for_sites(self.sf_candidate_list, extnd, self.n_jobs, self.sf_bam_ref, sf_flank_wfolder)

        ####m_site_fa_folder: s_id:(sf_tmp_fa, s_site_folder)
        print("Begin to extract the reads")
        m_site_fa_folder=self.extract_flanking_and_insertion_reads(extnd, bin_size, sf_reads_wfolder)##
        ####
        l_site_internal_vcfs=[]
        l_site_flanking_vcfs=[]
        for s_site_id in m_site_fa_folder:##
            (sf_site_fa, sf_site_folder)=m_site_fa_folder[s_site_id]
            if len(sf_site_folder)>0 and sf_site_folder[-1]!="/":
                sf_site_folder+="/"
            if os.path.exists(sf_site_folder)==False:
                os.mkdir(sf_site_folder)
            #align the reads to consensus
            sf_site_out_bam=sf_site_folder+"%s_algn_2_cns.bam"%s_site_id
            if os.path.isfile(sf_site_fa)==False:
                print("%s doesn't exist!"%sf_site_fa)##
                continue
            self._algn_collected_reads_one_site(sf_site_fa, sf_new_cns, sf_site_out_bam)
            #call insertion internal mutation from alignments
            sf_site_vcf=sf_site_folder + "{0}_{1}.vcf.gz".format(s_site_id, self.INTERNAL)
            self._call_mutation_one_site(sf_site_out_bam, sf_new_cns, i_dp_cutoff_internal, sf_site_vcf)###
            l_site_internal_vcfs.append(sf_site_vcf)

            #3.2 then extra the subset flanking reads, and align them to the flanking sequences
            sf_site_flanking_clip= sf_reads_wfolder + s_site_id + global_values.S_FLANKING_READS_of_CLIP
            sf_site_flanking_disc=sf_reads_wfolder + s_site_id + global_values.S_FLANKING_READS_of_DISC

            sf_slcted_flanking_fa=sf_reads_wfolder + s_site_id + "_slcted_flanking_reads.fa"
            self._further_slct_flanking_reads(sf_site_out_bam, sf_new_cns, sf_site_flanking_clip, sf_site_flanking_disc, sf_slcted_flanking_fa)
#           ####

            # 4.1 align reads to the flanking sequences
            sf_flank_fa = sf_flank_wfolder + "%s_flanks.fa" % s_site_id  ##
            if os.path.isfile(sf_flank_fa)==False:
                print("%s doesn't exist!"%sf_flank_fa)##
                continue
            sf_flank_bam=sf_flank_wfolder + "%s_reads_2_flanks.bam" % s_site_id
            self._algn_collected_reads_one_site(sf_slcted_flanking_fa, sf_flank_fa, sf_flank_bam)##

            #4.2 call flanking mutations
            sf_flank_vcf = sf_flank_wfolder + "{0}_{1}.vcf.gz".format(s_site_id, self.FLANKING)
            l_site_id=s_site_id.split(global_values.SEPERATOR)
            i_offset=int(l_site_id[1])-extnd
            self._call_mutation_one_site_with_offset(sf_flank_bam, sf_flank_fa, i_offset, i_dp_cutoff_flank, sf_flank_vcf)  ##
            l_site_flanking_vcfs.append(sf_flank_vcf)

#           ####clean the files#######################################################
            self.clean_file_by_path(sf_site_fa)
            self.clean_file_by_path(sf_site_out_bam)
            self.clean_file_by_path(sf_site_flanking_clip)
            self.clean_file_by_path(sf_site_flanking_disc)
            self.clean_file_by_path(sf_slcted_flanking_fa)
            self.clean_files_with_a_prefix(sf_slcted_flanking_fa+".*")
            self.clean_file_by_path(sf_flank_fa)
            self.clean_file_by_path(sf_flank_bam)

        ####merge the internal vcfs####
        sf_merged_internal_vcf_tmp=sf_merged_internal_vcf+".tmp"
        self._merge_vcfs_from_list(l_site_internal_vcfs, sf_merged_internal_vcf_tmp)
        self._rephrase_bcftools_merged_vcf(sf_merged_internal_vcf_tmp, sf_merged_internal_vcf)
        self.clean_file_by_path(sf_merged_internal_vcf_tmp)

        ####merge the flanking vcfs
        sf_merged_flanking_vcf_tmp=sf_merged_flanking_vcf+".tmp"
        self._merge_vcfs_from_list(l_site_flanking_vcfs, sf_merged_flanking_vcf_tmp)
        self._rephrase_bcftools_merged_flanking_vcf(sf_merged_flanking_vcf_tmp, sf_merged_flanking_vcf)
        self.clean_file_by_path(sf_merged_flanking_vcf_tmp)
####

    ####
    def _algn_collected_reads_one_site(self, sf_fa, sf_rep_cns, sf_out_bam):
        # # ##re-align the clipped reads
        bwa_align = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, self.n_jobs)
        #sf_algnmt = working_folder + "temp_clip.sam"
        #self.clean_file_by_path(sf_disc_fa)
        bwa_align.index_ref_file(sf_rep_cns)
        bwa_align.index_fa_file(sf_rep_cns) #samtools faidx
        bwa_align.realign_reads_to_bam(global_values.SAMTOOLS_PATH, sf_rep_cns, sf_fa, sf_out_bam)

    def _call_mutation_one_site(self, sf_bam, sf_cns, i_dp_cutoff, sf_vcf):
        self.call_snp_from_reads_alignment(sf_bam, sf_cns, i_dp_cutoff, sf_vcf)
    ####

    def _call_mutation_one_site_with_offset(self, sf_bam, sf_cns, i_offset, i_dp_cutoff, sf_vcf):
        self.call_snp_from_reads_alignment_with_offset(sf_bam, sf_cns, i_offset, i_dp_cutoff, sf_vcf)##
    ####

    ####
    def _further_slct_flanking_reads(self, sf_cns_algnmt, sf_cns_fa, sf_ori_clip_fa, sf_ori_disc_fa, sf_slcted_fa):
        m_qualified_reads={}
        samfile = pysam.AlignmentFile(sf_cns_algnmt, "rb", reference_filename=sf_cns_fa)
        for algnmt in samfile.fetch():
            if algnmt.is_duplicate == True:  ##duplciate
                continue
            if algnmt.is_unmapped == True:  # unmapped
                continue
            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  # wrong alignment
                continue

            s_id=algnmt.query_name
            l_id_fields=s_id.split(global_values.SEPERATOR)
            s_query_name=l_id_fields[6]
            m_qualified_reads[s_query_name]=1

        samfile.close()
        with open(sf_slcted_fa,"w") as fout:
            with pysam.FastxFile(sf_ori_clip_fa) as fin_clip:
                for entry in fin_clip:
                    s_rname = str(entry.name)
                    # e.g.@chr2~22968144~R~1~22968144~1~0
                    fout.write(">"+s_rname+"\n"+str(entry.sequence)+"\n")

            with pysam.FastxFile(sf_ori_disc_fa) as fin_disc:
                for entry in fin_disc:
                    s_rname = str(entry.name)
                    # e.g. >HISEQ1:48:HAC69ADXX:2:2201:21070:52952~0~1~0~1~0~22968127~chr2~22968144~0
                    fout.write(">" + s_rname + "\n" + str(entry.sequence) + "\n")


    def extract_flanking_and_insertion_reads(self, extnd, bin_size, sf_working_folder):##
        x_cd_filter = XClipDiscFilter(self.sf_bam_list, sf_working_folder, self.n_jobs, self.sf_bam_ref)
        ####
        # collect the clipped and discordant reads
        # each record in format like: @20~41951715~L~1~41952104~0~0,
        # (chrm, map_pos, global_values.FLAG_LEFT_CLIP, is_rc, insertion_pos, n_cnt_clip, sample_id)
        sf_clip_fq = sf_working_folder + "candidate_sites_all_clip.fq"
        sf_disc_fa = sf_working_folder + "candidate_sites_all_disc.fa"
        x_cd_filter.collect_clipped_disc_and_flanking_reads(self.sf_candidate_list, extnd, bin_size, sf_clip_fq, sf_disc_fa)  ##
        ####split the collected reads by site
        m_site_fa_folder=x_cd_filter.split_merged_reads_by_site(sf_clip_fq, sf_disc_fa, sf_working_folder)
        return m_site_fa_folder
####

    def call_mutations_from_reads_algnmt(self, sf_sites, sf_cns, i_len_cutoff, n_jobs, sf_merged_vcf):
        m_sites = self._select_high_quality_sites(sf_sites, i_len_cutoff)

        l_records = []
        for site_chrm in m_sites:
            for pos in m_sites[site_chrm]:
                l_records.append((site_chrm, pos, 1, sf_cns))  # here use 1 core for each job
        #align reads to cns, and collect the qualified reads
        self.align_reads_to_cns_in_parallel(l_records, n_jobs)
        #call the mutations from the alignments
        self.call_snp_in_parallel(l_records, n_jobs, sf_merged_vcf)


    ####
    def align_reads_to_cns_in_parallel(self, l_records, n_jobs):
        pool = Pool(n_jobs)
        pool.map(unwrap_self_algn_read, list(zip([self] * len(l_records), l_records)), 1)
        pool.close()
        pool.join()

    ####
    def call_snp_in_parallel(self, l_records, n_jobs, sf_merged_vcf):
        pool = Pool(n_jobs)
        pool.map(unwrap_self_call_snp, list(zip([self] * len(l_records), l_records)), 1)
        pool.close()
        pool.join()
        #merge all the vcf.gz files, and then index it
        #also remove those sites only happen in one sites
        s_mutation_flder = self.working_folder + global_values.MUTATION_FOLDER + "/"
        cmd="{0} merge ".format(global_values.BCFTOOLS_PATH)
        for rcd in l_records:
            site_chrm = rcd[0]
            pos = rcd[1]
            sf_vcf = s_mutation_flder + "{0}_{1}_{2}.vcf.gz".format(site_chrm, pos, global_values.ALL_HAP_SLCT)
            cmd+=(sf_vcf+" ")
        cmd+="-o {0}".format(sf_merged_vcf)
        self.cmd_runner.run_cmd_small_output(cmd)

    def _merge_vcfs(self, l_vcfs, sf_merged_vcf):
        cmd = "{0} merge ".format(global_values.BCFTOOLS_PATH)
        for sf_vcf in l_vcfs:
            cmd += (sf_vcf + " ")
        cmd += "-o {0}".format(sf_merged_vcf)
        self.cmd_runner.run_cmd_small_output(cmd)

    #this is to merge the files from a list
    #split to batches with 500 each (will fail if the number if big)
    def _merge_vcfs_from_list(self, l_vcfs, sf_merged_vcf):##
        sf_tmp_list=sf_merged_vcf+".list_for_merge" ##
        l_batch_list=[]
        l_batch_files=[]
        n_tmp_cnt=0
        i_batch=1
        n_total_vcf=len(l_vcfs)
        while n_tmp_cnt<n_total_vcf:
            n_tmp_end=n_tmp_cnt+self.n_batch_vcf
            if n_tmp_end>=n_total_vcf:
                n_tmp_end=n_total_vcf
            sf_batch_tmp=sf_tmp_list+"_%d" % i_batch

            with open(sf_batch_tmp, "w") as fout_tmp:
                for sf_vcf in l_vcfs[n_tmp_cnt:n_tmp_end]:
                    fout_tmp.write(sf_vcf + "\n")
            n_tmp_cnt=n_tmp_end

            #merge one batch
            cmd = "%s merge -l %s -Oz " % (global_values.BCFTOOLS_PATH, sf_batch_tmp)
            sf_merged_vcf_batch=sf_merged_vcf+"_%d" % i_batch
            cmd += "-o {0}".format(sf_merged_vcf_batch)
            self.cmd_runner.run_cmd_small_output(cmd)
            l_batch_files.append(sf_batch_tmp)
            #index one batch
            cmd="%s index %s" % (global_values.BCFTOOLS_PATH, sf_merged_vcf_batch)
            self.cmd_runner.run_cmd_small_output(cmd)
            l_batch_list.append(sf_merged_vcf_batch)
            i_batch += 1

        #merge the batches
        with open(sf_tmp_list, "w") as fout_tmp:
            for sf_tmp in l_batch_list:
                fout_tmp.write(sf_tmp+"\n")

        cmd = "%s merge -l %s " % (global_values.BCFTOOLS_PATH, sf_tmp_list)
        cmd += "-o {0}".format(sf_merged_vcf)
        self.cmd_runner.run_cmd_small_output(cmd)

        ####clean the files
        for sf_tmp in l_batch_files:
            with open(sf_tmp) as fin:
                for line in fin:
                    sf_tmp_vcf=line.rstrip()
                    self.clean_file_by_path(sf_tmp_vcf)#
                    self.clean_file_by_path(sf_tmp_vcf+".csi")  #
            self.clean_file_by_path(sf_tmp)
        for sf_tmp in l_batch_list:
            self.clean_file_by_path(sf_tmp)
            self.clean_file_by_path(sf_tmp + ".csi")

    def _rephrase_bcftools_merged_vcf(self, sf_old_vcf, sf_new_vcf):
        #this is mainly to change the format of the header
        with open(sf_old_vcf) as fin_old, open(sf_new_vcf,"w") as fout_new:
            for line in fin_old:
                if "#CHROM" not in line:
                    fout_new.write(line.rstrip()+"\n")
                else:
                    l_fields=line.rstrip().split("\t")
                    s_new_line="\t".join(l_fields[:9])
                    for s_rcd in l_fields[9:]:
                        l_rcd=s_rcd.split("/")
                        s_bam_name=l_rcd[-1]
                        l_bam_name=s_bam_name.split("_")
                        s_id=l_bam_name[0]
                        s_new_line+=("\t"+s_id)
                    fout_new.write(s_new_line+"\n")

    def _rephrase_bcftools_merged_flanking_vcf(self, sf_old_vcf, sf_new_vcf):##
        #this is mainly to change the format of the header
        with open(sf_old_vcf) as fin_old, open(sf_new_vcf,"w") as fout_new:
            for line in fin_old:
                if len(line)<1:
                    continue
                ####
                if "#CHROM" in line:
                    #
                    l_fields=line.rstrip().split("\t")
                    s_new_line="\t".join(l_fields[:9])
                    for s_rcd in l_fields[9:]:
                        l_rcd=s_rcd.split("/")
                        s_bam_name=l_rcd[-1]
                        l_bam_name=s_bam_name.split("_")
                        s_id=l_bam_name[0]
                        s_new_line+=("\t"+s_id)
                    fout_new.write(s_new_line+"\n")
                else:
                    if line[0]=="#":
                        fout_new.write(line.rstrip() + "\n")
                    else:
                        l_fields = line.rstrip().split("\t")
                        l_id=l_fields[0].split("_")
                        s_chrm=l_id[0]
                        s_new_line=s_chrm+"\t"+"\t".join(l_fields[1:])#
                        fout_new.write(s_new_line + "\n")

####
    def align_reads_to_cns_one_site(self, record):##
        site_chrm=record[0]
        pos=record[1]
        n_cores=record[2] #by default, is 1
        sf_cns=record[3]
        s_reads_folder=self.working_folder + global_values.READS_FOLDER + "/"
        bwa_align = BWAlign(global_values.BWA_PATH, global_values.BWA_REALIGN_CUTOFF, n_cores)
        sfa_hap1_slct = s_reads_folder + "{0}_{1}_{2}.fa".format(site_chrm, pos, global_values.HAP1_SLCT)
        sfa_hap2_slct = s_reads_folder + "{0}_{1}_{2}.fa".format(site_chrm, pos, global_values.HAP2_SLCT)
        sfa_hap_unknown_slct = s_reads_folder + "{0}_{1}_{2}.fa".format(site_chrm, pos, global_values.HAP_UNKNOWN_SLCT)
        sfa_hap_all_slct = s_reads_folder + "{0}_{1}_{2}.fa".format(site_chrm, pos, global_values.ALL_HAP_SLCT)

        s_mutation_flder = self.working_folder + global_values.MUTATION_FOLDER + "/"
        ssam_hap1_tmp = s_mutation_flder + "{0}_{1}_{2}.tmp.bam".format(site_chrm, pos, global_values.HAP1_SLCT)
        ssam_hap2_tmp = s_mutation_flder + "{0}_{1}_{2}.tmp.bam".format(site_chrm, pos, global_values.HAP2_SLCT)
        ssam_hap_unknown_tmp = s_mutation_flder + "{0}_{1}_{2}.tmp.bam".format(site_chrm, pos, global_values.HAP_UNKNOWN_SLCT)
        ssam_hap_all_tmp = s_mutation_flder + "{0}_{1}_{2}.tmp.bam".format(site_chrm, pos, global_values.ALL_HAP_SLCT)

        # align reads to the consensus
        bwa_align.realign_reads_to_bam(global_values.SAMTOOLS_PATH, sf_cns, sfa_hap1_slct, ssam_hap1_tmp)
        bwa_align.realign_reads_to_bam(global_values.SAMTOOLS_PATH, sf_cns, sfa_hap2_slct, ssam_hap2_tmp)
        bwa_align.realign_reads_to_bam(global_values.SAMTOOLS_PATH, sf_cns, sfa_hap_unknown_slct, ssam_hap_unknown_tmp)
        bwa_align.realign_reads_to_bam(global_values.SAMTOOLS_PATH, sf_cns, sfa_hap_all_slct, ssam_hap_all_tmp)

        ssam_hap1_slct = s_mutation_flder + "{0}_{1}_{2}.sorted.bam".format(site_chrm, pos, global_values.HAP1_SLCT)
        ssam_hap2_slct = s_mutation_flder + "{0}_{1}_{2}.sorted.bam".format(site_chrm, pos, global_values.HAP2_SLCT)
        ssam_hap_unknown_slct = s_mutation_flder + "{0}_{1}_{2}.sorted.bam".format(site_chrm, pos, global_values.HAP_UNKNOWN_SLCT)
        ssam_hap_all_slct = s_mutation_flder + "{0}_{1}_{2}.sorted.bam".format(site_chrm, pos, global_values.ALL_HAP_SLCT)
        self.select_qualified_alignments(ssam_hap1_tmp, sf_cns, ssam_hap1_slct)
        self.select_qualified_alignments(ssam_hap2_tmp, sf_cns, ssam_hap2_slct)
        self.select_qualified_alignments(ssam_hap_unknown_tmp, sf_cns, ssam_hap_unknown_slct)
        self.select_qualified_alignments(ssam_hap_all_tmp, sf_cns, ssam_hap_all_slct)

    ####
    def call_snp_one_site(self, record):
        site_chrm = record[0]
        pos = record[1]
        sf_cns = record[3]

        s_mutation_flder = self.working_folder + global_values.MUTATION_FOLDER + "/"
        ssam_hap_all_slct = s_mutation_flder + "{0}_{1}_{2}.sorted.bam".format(site_chrm, pos, global_values.ALL_HAP_SLCT)
        sf_vcf=s_mutation_flder+"{0}_{1}_{2}.vcf.gz".format(site_chrm, pos, global_values.ALL_HAP_SLCT)
        i_dp_cutoff=8
        self.call_snp_from_reads_alignment(ssam_hap_all_slct, sf_cns, i_dp_cutoff, sf_vcf)


    #realign_reads_to_bam(self, SAMTOOLS, sf_ref, sf_reads, sf_out_bam)
    #call out the snp and indels from the haplotype-based alignment
    #1. Call out all 1/1 snps for each haplotype, and get the genotype for each snp
    #2. Create a matrix all the insertions.
    #output in vcf.gz format
    def call_snp_from_reads_alignment(self, sf_algnmt, sf_ref, i_dp_cutoff, sf_vcf):
        sf_vcf_tmp=sf_vcf+".tmp"
        cmd="{0} mpileup -IOu -f {1} {2} | {3} call -vmO z -o {4}".format(global_values.BCFTOOLS_PATH, sf_ref, sf_algnmt,
                                                                     global_values.BCFTOOLS_PATH, sf_vcf_tmp)
        self.cmd_runner.run_cmd_small_output(cmd)
        ####index the vcf gz file
        self._index_vcf_gz(sf_vcf_tmp)
        ####filter the results  Hard code here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        cmd="{0} filter -O z -e 'DP<{1} || %QUAL<50' -o {2} {3}".format(global_values.BCFTOOLS_PATH, i_dp_cutoff, sf_vcf, sf_vcf_tmp)
        self.cmd_runner.run_cmd_small_output(cmd)
        ####
        self._index_vcf_gz(sf_vcf)
####

    def call_snp_from_reads_alignment_with_offset(self, sf_algnmt, sf_ref, i_offset, i_dp_cutoff, sf_vcf):
        sf_vcf_tmp=sf_vcf+".tmp"
        cmd="{0} mpileup -IOu -f {1} {2} | {3} call -vmO v -o {4}".format(global_values.BCFTOOLS_PATH, sf_ref, sf_algnmt,
                                                                     global_values.BCFTOOLS_PATH, sf_vcf_tmp)
        self.cmd_runner.run_cmd_small_output(cmd)
        ####index the vcf gz file
        #self._index_vcf_gz(sf_vcf_tmp)

        sf_vcf_tmp2=sf_vcf+".tmp2"
        #add the offset for each called mutation
        self._add_offset_to_vcf(sf_vcf_tmp, i_offset, sf_vcf_tmp2)

        ####filter the results  Hard code here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        cmd="{0} filter -O z -e 'DP<{1} || %QUAL<50' -o {2} {3}".format(global_values.BCFTOOLS_PATH, i_dp_cutoff, sf_vcf, sf_vcf_tmp2)
        self.cmd_runner.run_cmd_small_output(cmd)

        #index the gzipped vcf
        self._index_vcf_gz(sf_vcf)
#

    def _add_offset_to_vcf(self, sf_old_vcf, i_offset, sf_new_vcf):
        with open(sf_new_vcf,"w") as fout, open(sf_old_vcf) as fin:
            for line in fin:
                if len(line)>0:
                    if line[0]=="#":
                        fout.write(line.rstrip()+"\n")
                    else:
                        fields=line.rstrip().split("\t")
                        s_new_line=fields[0]+"\t"+str(int(fields[1])+i_offset)
                        for s_rcd in fields[2:]:
                            s_new_line+=("\t"+s_rcd)
                        fout.write(s_new_line+"\n")

    ####Index xxx.vcf.gz file
    def _index_vcf_gz(self, sf_input):
        #cmd = "{0} -f -p vcf {1}".format(global_values.TABIX_PATH, sf_input)
        cmd = "{0} index {1}".format(global_values.BCFTOOLS_PATH, sf_input)
        #Popen(cmd, shell=True, stdout=PIPE).communicate()
        self.cmd_runner.run_cmd_small_output(cmd)

    ####
    def is_fully_mapped(self, l_cigar, f_cutoff):
        cnt_map = 0
        cnt_total = 0
        for cigar in l_cigar:
            cigar_opr = int(cigar[0])
            cigar_lth = int(cigar[1])
            if cigar_opr == 0:
                cnt_map += cigar_lth
            if cigar_opr != 2:
                cnt_total += cigar_lth
        if cnt_total == 0:
            return False, cnt_map
        if float(cnt_map) / float(cnt_total) >= f_cutoff:
            return True, cnt_map
        return False, cnt_map

    # select reads aligned to consensus (with no or little mutations)
    def select_qualified_alignments(self, sf_algnmt, sf_ref, sf_out_bam):
        bam_in = pysam.AlignmentFile(sf_algnmt, "rb", reference_filename=sf_ref)
        bam_out = pysam.Samfile(sf_out_bam, 'wb', template=bam_in)
        for alignment in bam_in.fetch():
            if alignment.is_unmapped == True:  # unmapped
                continue
            if alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.is_duplicate == True:  ##duplciate
                continue
            rlth = len(alignment.query_sequence)
            l_cigar=alignment.cigar
            ####here select fully mapped and at most 3 mismatch
            #####Hard code here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            b_full_map, cnt_map=self.is_fully_mapped(l_cigar, 0.93)
            if b_full_map==False:
                continue
            if alignment.has_tag("MD") == False:
                continue
            s_tag = alignment.get_tag("MD")
            cnt_mismatch = 0
            for ctmp in s_tag:
                if ctmp >= "A" and ctmp <= "Z":
                    cnt_mismatch += 1
            if rlth>=100:
                if cnt_mismatch>2:
                    continue
            if rlth>=50 and rlth<100:
                if cnt_mismatch>1:
                    continue
            if rlth<50:
                if cnt_mismatch>1:
                    continue
            bam_out.write(alignment)
        bam_out.close()
        bam_in.close()
    ###################################################################################
    ####
    ####

    # 0T38G7G1
    # 49M
    # GCTCTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTCTGTGTGTTT
    # [(0, 803, 't'), (1, 804, 'C'), (2, 805, 'T'), (3, 806, 'C'), (4, 807, 'T'), (5, 808, 'G'), (6, 809, 'T'),
    #  (7, 810, 'G'), (8, 811, 'T'), (9, 812, 'G'), (10, 813, 'T'), (11, 814, 'G'), (12, 815, 'T'), (13, 816, 'G'),
    #  (14, 817, 'T'), (15, 818, 'G'), (16, 819, 'T'), (17, 820, 'G'), (18, 821, 'T'), (19, 822, 'G'), (20, 823, 'T'),
    #  (21, 824, 'G'), (22, 825, 'T'), (23, 826, 'G'), (24, 827, 'T'), (25, 828, 'G'), (26, 829, 'T'), (27, 830, 'G'),
    #  (28, 831, 'T'), (29, 832, 'G'), (30, 833, 'T'), (31, 834, 'G'), (32, 835, 'T'), (33, 836, 'G'), (34, 837, 'T'),
    #  (35, 838, 'G'), (36, 839, 'T'), (37, 840, 'G'), (38, 841, 'T'), (39, 842, 'g'), (40, 843, 'T'), (41, 844, 'G'),
    #  (42, 845, 'T'), (43, 846, 'G'), (44, 847, 'T'), (45, 848, 'G'), (46, 849, 'T'), (47, 850, 'g'), (48, 851, 'T')]
####

    ####Call out the snp from repeat copy alignments: From the MD field
    def call_snp_indel_from_rep_copy_algnmt(self):
        sf_bam = "/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/NA12878_10x_v2/tmp_dbg/cns/picked_disc.sorted.bam"
        sf_reference = "/n/data1/hms/dbmi/park/SOFTWARE/LongRanger/refdata-b37-2.1.0/fasta/genome.fa"
        bamfile = pysam.AlignmentFile(sf_bam, "r", reference_filename=sf_reference)
        ncnt = 0
        for algnmt in bamfile.fetch():
            l_mismatch = algnmt.get_aligned_pairs(False, True)
            if algnmt.has_tag("MD") == False:
                continue
            print(algnmt.get_tag("MD"))
            print(algnmt.cigarstring)
            print(algnmt.query_sequence)
            print(l_mismatch)
            if ncnt > 20:
                break
            ncnt += 1
        bamfile.close()
####

    def clean_file_by_path(self, sf_path):
        if os.path.isfile(sf_path)==True:
            os.remove(sf_path)

    #e.g. prefix path like "./xxx.*"
    def clean_files_with_a_prefix(self, sf_prefix):
        for filename in glob.glob(sf_prefix):
            os.remove(filename)
####