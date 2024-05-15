##05/08/2024
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com
##@@maintainer: Corinne Sexton, DBMI, Harvard Medical School
##@@contact: corinne_sexton@hms.harvard.edu
##

####
import os

import xtea.global_values
from xtea.x_TEI_locator import TE_Multi_Locator
from xtea.x_intermediate_sites import XIntermediateSites
from xtea.x_basic_info import X_BasicInfo
from xtea.x_parameter import Parameters


def automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs, b_force=False, b_tumor=False, f_purity=0.45):
    ####1. collect the basic information
    search_win = 500
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    rcd=x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win, b_force)
    f_cov=rcd[0]
    # rlth=rcd[1]
    # mean_is=rcd[2]
    # std_var=rcd[3]

    ####2. based on the coverage, set the parameters
    xpar=Parameters()
    if b_tumor is True:
        f_cov=f_cov*f_purity
    par_rcd=xpar.get_par_by_cov(f_cov) #in format (iclip, idisc, i_clip-disc)
    print("\t\tAve coverage is {0}: automatic parameters (clip, disc, clip-disc) with value ({1}, {2} ,{3})\n"
          .format(f_cov, par_rcd[0], par_rcd[1], par_rcd[2]))
    return par_rcd, rcd


def get_candidate_sites(options,annot_path_dict,output_dir,tmp_dir):

# OPTIONS STILL MISSING
    # if options.mit: #if this to call mitochondrial insertion, then will not filter out chrM in "x_intermediate_sites.py"
    #     xtea.global_values.turn_on_mit()
    # if options.dna:
    #     xtea.global_values.turn_off_rna_mediated()
    # if options.cbs:
    #     xtea.global_values.turn_on_check_by_sample()
    # if options.sva:
    #     xtea.global_values.turn_on_sva()

    if options.int_files:
        xtea.global_values.keep_intermediate_files()

    b_tumor=options.tumor #whether this is tumor sample
    f_purity=options.purity #tumor purity, by default 0.45
    b_se = options.single  ##single end reads or not, default is not
    b_mosaic=options.mosaic #this is for mosaic calling from normal tissue

    # NOT SURE THIS IS USED
    # site_clip_cutoff=options.siteclip #this is the cutoff for the exact position, use larger value for 10X
    # xtea.global_values.set_initial_min_clip_cutoff(site_clip_cutoff)

    ###take in the normal illumina reads (10x will be viewed as normal illumina)
    print("Working on \"clipped reads\" step!")

    sf_bam_list = options.input_bams
    n_jobs = int(options.cores)
    b_force = False # removed command line option
    b_resume=options.resume #resume the running, which will skip the step if output file already exists
    
    sf_rep_cns = annot_path_dict['sf_rep_cns']
    sf_rep = annot_path_dict['sf_rep'] #repeat copies "-r"
    sf_annotation = annot_path_dict['sf_annotation']
    sf_ref = annot_path_dict['sf_ref'] #reference genome "-ref"

    s_working_folder = output_dir
    sf_out = f"{s_working_folder}/candidate_list_from_clip.txt"
    
    wfolder_pub_clip = tmp_dir #public clip folder

    # downstream NOT USED (filler values for now, remove downstream later)
    cutoff_left_clip = 10
    cutoff_right_clip = 10
    
    # true clipping cutoff
    cutoff_clip_mate_in_rep = options.cr

    rcd = None
    basic_rcd = None
    if b_resume is False or os.path.isfile(sf_out) is False:
        print("\tGenerating cutoff parameters based on coverage.")
        if cutoff_clip_mate_in_rep is None:
            rcd, basic_rcd=automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs,
                                                        b_force, b_tumor, f_purity)
            cutoff_clip_mate_in_rep=rcd[2]

        tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)



        #by default, if number of clipped reads is larger than this value, then discard
        # max_cov_cutoff=int(15*basic_rcd[0])   # ORIGINALLY WAS A CL PARAMETER (cov) set to 40 by default
        max_cov_cutoff=600

        ##Hard code inside:
        # 1. call_TEI_candidate_sites_from_clip_reads_v2 --> run_cnt_clip_part_aligned_to_rep_by_chrm_sort_version
        # here if half of the seq is mapped, then consider it as aligned work.
        ##2. require >=2 clip reads, whose clipped part is aligned to repeat copies
        print("\tCalling insertion sites.")
        tem_locator.call_TEI_candidate_sites_from_multiple_alignmts(sf_annotation, sf_rep_cns, sf_rep, b_se,
                                                                    cutoff_left_clip, cutoff_right_clip,
                                                                    cutoff_clip_mate_in_rep, b_mosaic,
                                                                    wfolder_pub_clip, b_force, max_cov_cutoff, sf_out)

    print("Working on \"discordant reads\" step!")
    sf_candidate_list = sf_out

    sf_disc_out = f"{s_working_folder}/candidate_list_from_disc.txt"
    PEAK_WINDOW = 100
    if b_mosaic: #for mosaic/somatic events # CS EDIT!
        PEAK_WINDOW = 30

    if b_resume == False or os.path.isfile(sf_disc_out) == False:
        xfilter = XIntermediateSites()
        m_original_sites = xfilter.load_in_candidate_list(sf_candidate_list)
        sf_peak_sites = f"{s_working_folder}/clip_peak_candidate.list"
        # get the peak sites
        m_sites_clip_peak = xfilter.call_peak_candidate_sites_with_std_derivation(m_original_sites, PEAK_WINDOW)
        xfilter.output_candidate_sites(m_sites_clip_peak, sf_peak_sites)  # output the sites
        m_original_sites.clear()  #release the memory

        if rcd is None or basic_rcd is None:
            rcd, basic_rcd = automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs,
                                                        b_force, b_tumor, f_purity)
        rlth = basic_rcd[1]  # read length
        mean_is = basic_rcd[2]  # mean insert size
        std_var = basic_rcd[3]  # standard derivation
            
        max_is = int(mean_is + 3 * std_var) + int(rlth)
        iextend = max_is
        i_is = 100000  ###set the insert size a large value, by default 100k
        f_dev = std_var

        # this is the cutoff for  "left discordant" and "right discordant"
        # Either of them is larger than this cutoff, the site will be reported
        n_disc_cutoff = options.nd
        if n_disc_cutoff is None:
            n_disc_cutoff=rcd[1]
            
        sf_tmp = s_working_folder + "/disc_tmp.list"
        sf_raw_disc=sf_disc_out + xtea.global_values.RAW_DISC_TMP_SUFFIX #save the left and right raw disc for each site
        tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)
        print(f"\tFiltering insertion sites based on discordant cutoff: {n_disc_cutoff}.")
        tem_locator.filter_candidate_sites_by_discordant_pairs_multi_alignmts(m_sites_clip_peak, iextend, i_is,
                                                                                f_dev, n_disc_cutoff, sf_annotation,
                                                                                sf_tmp, sf_raw_disc, b_tumor)
        xfilter.merge_clip_disc(sf_tmp, sf_candidate_list, sf_disc_out)