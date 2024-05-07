import pytest

# CS COMMENTED OUT 5/7/24

@pytest.fixture
def TE_Locator_test():
    return TELocator(sf_bam, sf_barcode_bam, s_working_folder, n_jobs, sf_ref)
    
    # def __init__(self, sf_bam, sf_barcode_bam, s_working_folder, n_jobs, sf_ref):
    #     self.sf_bam = sf_bam
    #     self.sf_barcode_bam = sf_barcode_bam
    #     self.working_folder = s_working_folder
    #     self.n_jobs = int(n_jobs)
    #     self.sf_reference = sf_ref  ##reference genome
    #     self.cmd_runner = CMD_RUNNER()


# call_TEI_candidate_sites_from_clip_reads_v2(self, sf_annotation, sf_rep_cns, sf_ref, b_se, cutoff_left_clip,
#                                                     cutoff_right_clip, b_cutoff, sf_pub_folder, idx_bam,
#                                                     b_force, max_cov_cutoff, sf_out):

# def call_TEI_candidate_sites_from_clip_reads_v2_mosaic(self, sf_annotation, sf_rep_cns, sf_ref, b_se,
#                                                     cutoff_left_clip,
#                                                     cutoff_right_clip, b_cutoff, sf_pub_folder, idx_bam,
#                                                     b_force, max_cov_cutoff, sf_out):
    
# collect_all_clipped_reads_only(self, sf_annotation, b_se, s_working_folder, wfolder_pub_clip):        

# run_filter_by_barcode_coverage(self, record):

# filter_candidate_sites_by_barcode_coverage(self, m_candidate_sites, iextend, i_cov_cutoff):

# run_filter_by_discordant_pair_by_chrom_non_barcode(self, record):
        # _get_chrm_id_name(self, samfile): (ONLY USED ONCE)
# filter_candidate_sites_by_discordant_pairs_non_barcode(self, m_candidate_sites, iextend, i_is, f_dev,
#                                                                sf_annotation, n_discordant_cutoff):
       
# output_candidate_sites_by_chrm(self, m_candidate_list, sf_folder, s_suffix):
      
# output_candidate_sites(self, m_candidate_list, sf_out):



#######
def f():
    raise SystemExit(1)


def test_mytest():
    with pytest.raises(SystemExit):
        f()