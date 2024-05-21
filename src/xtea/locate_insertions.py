##05/08/2024
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong.simon.chu@gmail.com
##@@maintainer: Corinne Sexton, DBMI, Harvard Medical School
##@@contact: corinne_sexton@hms.harvard.edu
##

####
import os
from shutil import copyfile

import xtea.global_values
from xtea.x_TEI_locator import TE_Multi_Locator
import xtea.x_annotation
from xtea.x_intermediate_sites import XIntermediateSites
from xtea.x_basic_info import X_BasicInfo
from xtea.x_parameter import Parameters
from xtea.x_clip_disc_filter import XClipDiscFilter
from xtea.x_transduction import XTransduction
from xtea.x_orphan_transduction import XOrphanTransduction
from xtea.x_mosaic_calling import MosaicCaller
from xtea.x_post_filter import XPostFilter
from xtea.x_gvcf import gVCF
from xtea.x_gene_annotation import GFF3



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


def get_clip_sites(options,annot_path_dict,output_dir, wfolder_pub_clip):

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
    
    # downstream USED only in (filler values for now, remove downstream later)
    cutoff_left_clip = options.clip_cutoff
    cutoff_right_clip = options.clip_cutoff
    
    # true clipping cutoff
    cutoff_clip_mate_in_rep = options.cr

    print("\tGenerating cutoff parameters based on coverage.")
    # ALWAYS ATTEMPT TO CALCULATE BECAUSE NEEDED DOWNSTREAM
    rcd, basic_rcd=automatic_gnrt_parameters(sf_bam_list, sf_ref, s_working_folder, n_jobs,
                                                    b_force, b_tumor, f_purity)
    if cutoff_clip_mate_in_rep is None: 
        cutoff_clip_mate_in_rep=rcd[2]
    else:
        rcd[2] = cutoff_clip_mate_in_rep # set to user preset
    if cutoff_left_clip is None:
        cutoff_left_clip=rcd[0]
        cutoff_right_clip=rcd[0]
    else:
        rcd[0] = cutoff_left_clip # set to user preset


    if b_resume is False or os.path.isfile(sf_out) is False:
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
        
    return (rcd,basic_rcd)

    
def get_disc_sites(options,annot_path_dict,output_dir,rcd,basic_rcd):

    print("Working on \"discordant reads\" step!")

    b_tumor=options.tumor #whether this is tumor sample
    b_mosaic=options.mosaic #this is for mosaic calling from normal tissue

    sf_bam_list = options.input_bams
    n_jobs = int(options.cores)
    b_resume=options.resume #resume the running, which will skip the step if output file already exists
    
    sf_annotation = annot_path_dict['sf_annotation']
    sf_ref = annot_path_dict['sf_ref'] #reference genome "-ref"

    s_working_folder = output_dir

    sf_candidate_list = f"{s_working_folder}/candidate_list_from_clip.txt"

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
            n_disc_cutoff = rcd[1]
            rcd[1] = n_disc_cutoff
            
        sf_tmp = s_working_folder + "/disc_tmp.list"
        sf_raw_disc=sf_disc_out + xtea.global_values.RAW_DISC_TMP_SUFFIX #save the left and right raw disc for each site
        tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)
        print(f"\tFiltering insertion sites based on discordant cutoff: {n_disc_cutoff}.")
        tem_locator.filter_candidate_sites_by_discordant_pairs_multi_alignmts(m_sites_clip_peak, iextend, i_is,
                                                                                f_dev, n_disc_cutoff, sf_annotation,
                                                                                sf_tmp, sf_raw_disc, b_tumor)
        xfilter.merge_clip_disc(sf_tmp, sf_candidate_list, sf_disc_out)


def filter_csn(options,annot_path_dict,output_dir,rcd,basic_rcd):
    print("Working on \"clip-disc-filtering\" step!")

    sf_bam_list = options.input_bams
    n_jobs = int(options.cores)
    b_resume=options.resume #resume the running, which will skip the step if output file already exists
    
    sf_ref = annot_path_dict['sf_ref'] # reference genome "-ref"
    sf_cns = annot_path_dict['sf_rep'] # cns ref "-r"

    s_working_folder = output_dir
    print("Current working folder is: {0}\n".format(s_working_folder))

    sf_candidate_list = f"{s_working_folder}/candidate_list_from_disc.txt"
    sf_output = f"{s_working_folder}/candidate_disc_filtered_cns.txt"

    iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
    bin_size = 50000000  # block size for parallelization
    i_concord_dist = 550  # this should be the 3*std_derivation, used to cluster disc reads on the consensus
    f_concord_ratio = 0.45

    bmapped_cutoff = xtea.global_values.MIN_CLIP_MAPPED_RATIO
    
    if b_resume == False or os.path.isfile(sf_output) == False:
        ave_cov = basic_rcd[0]  # ave coverage
        rlth = basic_rcd[1]  # read length
        mean_is = basic_rcd[2]  # mean insert size
        std_var = basic_rcd[3]  # standard derivation

        print("Mean insert size is: {0}\n".format(mean_is))
        print("Standard derivation is: {0}\n".format(std_var))

        max_is = int(mean_is + 3 * std_var)
        if iextnd < max_is: #correct the bias
            iextnd = max_is
        if i_concord_dist < max_is: #correct the bias
            i_concord_dist = max_is

        xtea.global_values.set_read_length(rlth)
        xtea.global_values.set_insert_size(max_is)
        xtea.global_values.set_average_cov(ave_cov)

        print("Read length is: {0}\n".format(rlth))
        print("Maximum insert size is: {0}\n".format(max_is))
        print("Average coverage is: {0}\n".format(ave_cov))

        n_clip_cutoff=rcd[0]
        n_disc_cutoff=rcd[1]

        print("Filter (on cns) cutoff: {0} and {1} are used!!!".format(n_clip_cutoff, n_disc_cutoff))

        x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_ref)
        x_cd_filter.call_MEIs_consensus(sf_candidate_list, iextnd, bin_size, sf_cns, bmapped_cutoff, i_concord_dist, f_concord_ratio,
                                        n_clip_cutoff, n_disc_cutoff, sf_output)
        ####        

def get_transduction(r,options,annot_path_dict,output_dir,rcd,basic_rcd):

    b_tumor=options.tumor #whether this is tumor sample

    sf_bam_list = options.input_bams
    n_jobs = int(options.cores)
    b_resume=options.resume #resume the running, which will skip the step if output file already exists

   
    sf_flank = annot_path_dict['sf_flank']
    sf_reference = annot_path_dict['sf_ref'] #reference genome "-ref"
    sf_rmsk = annot_path_dict["sf_anno1"]
    sf_cns = annot_path_dict['sf_rep']  ####repeat copies/cns here

    s_working_folder = output_dir
    i_rep_type = r

    iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
    bin_size = 50000000  # block size for parallelization
    i_flank_lenth = 3000

    #need to re-collect all the clip, disc reads
    sf_candidate_list = f"{s_working_folder}/candidate_disc_filtered_cns.txt" #this is the output from the "cns" step.
    sf_raw_disc = f"{s_working_folder}/candidate_list_from_disc.txt.clip_sites_raw_disc.txt"  # this is the raw disc file
    sf_output = f"{s_working_folder}/candidate_disc_filtered_cns2.txt"

    print("Current working folder is: {0}\n".format(s_working_folder))

    if b_resume == False or os.path.isfile(sf_output) == False:
        if os.path.isfile(sf_flank)==True:#for Alu and many others, there is no transduction
            ave_cov = basic_rcd[0]  # ave coverage
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard derivation
            print("Mean insert size is: {0}\n".format(mean_is))
            print("Standard derivation is: {0}\n".format(std_var))

            max_is = int(mean_is + 3 * std_var)
            if iextnd < max_is:  # correct the bias
                iextnd = max_is
            i_concord_dist=550
            f_concord_ratio = 0.25
            if i_concord_dist < max_is:  # correct the bias
                i_concord_dist = max_is
                
            xtea.global_values.set_read_length(rlth)
            xtea.global_values.set_insert_size(max_is)
            xtea.global_values.set_average_cov(ave_cov)

            n_clip_cutoff = rcd[0]
            n_disc_cutoff = rcd[1]

            xtransduction = XTransduction(s_working_folder, n_jobs, sf_reference)

            i_win=150 #if a site is close to an existing site, then will not be considered again
            sf_tmp_slct=sf_raw_disc+".slct"
            i_min_copy_len = 225
            xannotation = xtransduction.prep_annotation_interval_tree(sf_rmsk, i_min_copy_len)
            xtransduction.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_candidate_list, n_disc_cutoff, xannotation,
                                            i_rep_type, b_tumor, sf_tmp_slct)
            #now for the selected sites, re-evaluate each one
            x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_reference)
            i_max_cov=ave_cov*(xtea.global_values.MAX_COV_TIMES+1)
            sf_output_tmp=sf_output + xtea.global_values.TD_NON_SIBLING_SUFFIX
            xtransduction.call_candidate_transduction_v3(sf_tmp_slct, sf_candidate_list, x_cd_filter,
                                                            sf_flank, sf_cns, i_flank_lenth, iextnd, bin_size, n_clip_cutoff,
                                                            n_disc_cutoff, i_concord_dist, f_concord_ratio, xannotation,
                                                            sf_bam_list, i_rep_type, i_max_cov, ave_cov, sf_output_tmp)
####
            xorphan=XOrphanTransduction(s_working_folder, n_jobs, sf_reference)
            n_half_disc_cutoff=n_disc_cutoff/2
            i_search_win=2000
            sf_updated_cns=sf_output #this is the final updated

            #1.Call out the sibling transduction events from the current list
            sf_sibling_TD=sf_output+".sibling_transduction_from_existing_list"#
            xorphan.call_sibling_TD_from_existing_list(sf_output_tmp, sf_bam_list, iextnd, n_half_disc_cutoff,
                                                        i_search_win, xannotation, i_rep_type, i_max_cov,
                                                        sf_updated_cns, sf_sibling_TD)

            #update high confident ones (in "cns" filter step, two results are generated)
            sf_ori_hc=sf_candidate_list+xtea.global_values.HIGH_CONFIDENT_SUFFIX
            sf_new_hc=sf_output+xtea.global_values.HIGH_CONFIDENT_SUFFIX
            xorphan.update_high_confident_callset(sf_ori_hc, sf_updated_cns, sf_new_hc)
####
        else:#rename the two files generated in previous step
            copyfile(sf_candidate_list, sf_output)
            sf_ori_hc = sf_candidate_list + xtea.global_values.HIGH_CONFIDENT_SUFFIX
            sf_new_hc=sf_output+xtea.global_values.HIGH_CONFIDENT_SUFFIX
            copyfile(sf_ori_hc, sf_new_hc)

def get_sibling(r,options,annot_path_dict,output_dir,rcd,basic_rcd):
    '''
    Todo: 09-29-2019: Add filtering modules:
    1. using background low mapq reads (multiple mapped reads) for filtering
    2. Set a upper-bound cutoff for discordant reads
    3. using blacklist for filtering
    '''

    b_tumor=options.tumor #whether this is tumor sample

    sf_bam_list = options.input_bams
    n_jobs = int(options.cores)

   
    sf_flank = annot_path_dict['sf_flank']
    sf_reference = annot_path_dict['sf_ref'] #reference genome "-ref"
    sf_rmsk = annot_path_dict["sf_anno1"]

    s_working_folder = output_dir
    i_rep_type = r

    #need to re-collect all the clip, disc reads
    sf_raw_disc = f"{s_working_folder}/candidate_list_from_disc.txt.clip_sites_raw_disc.txt"  # this is the raw disc file
    sf_pre_step_out = f"{s_working_folder}/candidate_disc_filtered_cns2.txt"
    sf_output = f"{s_working_folder}/candidate_sibling_transduction2.txt"

    # sf_black_list = options.blacklist
    sf_black_list = 'NULL' # NOT USED YET??

    if os.path.isfile(sf_flank) == True:  #for Alu and many others, there is no transduction
        ave_cov = basic_rcd[0]  # ave coverage
        rlth = basic_rcd[1]  # read length
        mean_is = basic_rcd[2]  # mean insert size
        std_var = basic_rcd[3]  # standard derivation
        print("Mean insert size is: {0}\n".format(mean_is))
        print("Standard derivation is: {0}\n".format(std_var))
        max_is = int(mean_is + 3 * std_var)

        i_concord_dist = 550
        if i_concord_dist < max_is:  # correct the bias
            i_concord_dist = max_is
        xtea.global_values.set_read_length(rlth)
        xtea.global_values.set_insert_size(mean_is) #here set mean inset size
        xtea.global_values.set_average_cov(ave_cov)

        n_clip_cutoff = rcd[0]
        n_disc_cutoff = rcd[1]

        xorphan = XOrphanTransduction(s_working_folder, n_jobs, sf_reference)
        i_min_copy_len = 225
        xorphan.set_boundary_extend(mean_is)
        xannotation = xorphan.prep_annotation_interval_tree(sf_rmsk, i_min_copy_len)
        # 2. Call orphan "sibling" transdcution from non_existing list
        # sf_sibling_TD2 = sf_output + ".novel_sibling_transduction"
        b_with_original = False
        sf_tmp_slct2 = sf_raw_disc + ".slct2"
        sf_output_tmp = sf_pre_step_out + xtea.global_values.TD_NON_SIBLING_SUFFIX
        # select the sites to exclude the already called out sites, and filter out sites fall in black_list
        xorphan.re_slct_with_clip_raw_disc_sites(sf_raw_disc, sf_output_tmp, n_disc_cutoff, xannotation,
                                                i_rep_type, b_tumor, sf_tmp_slct2, b_with_original)

        #re-select transduction candidates based on disc-clip consistency (clip position encompass disc ones?)
        m_failed_ori_td=xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list, i_concord_dist,
                                                                        n_clip_cutoff, n_disc_cutoff, sf_black_list,
                                                                        sf_rmsk)
        sf_sibling_TD=sf_output
        if os.path.isfile(sf_black_list)==False:
            print("Blacklist file {0} does not exist!".format(sf_black_list))
        xorphan.call_novel_sibling_TD_from_raw_list(sf_tmp_slct2, sf_bam_list, i_concord_dist, n_clip_cutoff,
                                                    n_disc_cutoff, sf_black_list, sf_rmsk, sf_sibling_TD)

        ####append the newly called events to existing list
        xorphan.append_to_existing_list(sf_sibling_TD, sf_pre_step_out, m_failed_ori_td)
        xorphan.append_to_existing_list(sf_sibling_TD, sf_pre_step_out+xtea.global_values.HIGH_CONFIDENT_SUFFIX,
                                        m_failed_ori_td)
    else:
        if os.path.isfile(sf_pre_step_out)==True:#do td filtering only
            ave_cov = basic_rcd[0]  # ave coverage
            rlth = basic_rcd[1]  # read length
            mean_is = basic_rcd[2]  # mean insert size
            std_var = basic_rcd[3]  # standard derivation
            print("Mean insert size is: {0}\n".format(mean_is))
            print("Standard derivation is: {0}\n".format(std_var))
            max_is = int(mean_is + 3 * std_var)

            i_concord_dist = 550
            if i_concord_dist < max_is:  # correct the bias
                i_concord_dist = max_is
            xtea.global_values.set_read_length(rlth)
            xtea.global_values.set_insert_size(mean_is)  # here set mean inset size
            xtea.global_values.set_average_cov(ave_cov)
            
            n_clip_cutoff = rcd[0]
            n_disc_cutoff = rcd[1]
            xorphan = XOrphanTransduction(s_working_folder, n_jobs, sf_reference)
            # re-select transduction candidates based on disc-clip consistency (clip position encompass disc ones?)
            m_failed_ori_td = xorphan.distinguish_source_from_insertion_for_td(sf_pre_step_out, sf_bam_list,
                                                                                i_concord_dist, n_clip_cutoff,
                                                                                n_disc_cutoff, sf_black_list, sf_rmsk)
            xorphan.update_existing_list_only(sf_pre_step_out, m_failed_ori_td)

def filter_sites_post(r,options,annot_path_dict,output_dir,basic_rcd):  
    
    b_tumor=options.tumor #whether this is tumor sample
    n_jobs = int(options.cores)

    sf_rmsk = annot_path_dict["sf_anno1"]

    s_working_folder = output_dir
    i_rep_type = r

    sf_new_out = f"{s_working_folder}/candidate_disc_filtered_cns_post_filtering.txt"
    sf_xtea_rslt = f"{s_working_folder}/candidate_disc_filtered_cns2.txt"

    # sf_black_list = options.blacklist
    sf_black_list = 'NULL' # NOT USED YET?? TODO

    b_pf_mosaic=options.mosaic #this is for mosaic calling from normal tissue
    i_min_copy_len=225 #when check whether fall in repeat region, require the minimum copy length
    
    if b_pf_mosaic is True:#for mosaic events
        xpf_mosic = MosaicCaller(s_working_folder, n_jobs)
        xpf_mosic.run_call_mosaic(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, sf_black_list, sf_new_out)
    else:
        f_cov = basic_rcd[0] # CS EDIT
        xpost_filter = XPostFilter(s_working_folder, n_jobs)
        #here sf_black_list is the centromere + duplication region
        xpost_filter.run_post_filtering(sf_xtea_rslt, sf_rmsk, i_min_copy_len, i_rep_type, f_cov, sf_black_list,
                                        sf_new_out, b_tumor)


def annotate_genes(options,output_dir):

    s_working_folder = output_dir
    sf_input=f"{s_working_folder}/candidate_disc_filtered_cns.txt.high_confident.post_filtering.txt"
    sf_output=f"{s_working_folder}/candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene.txt"

    sf_gene_annotation=options.genome_gff3

    gff=GFF3(sf_gene_annotation)
    iextnd=xtea.global_values.UP_DOWN_GENE
    i_user_extnd = options.extend
    if i_user_extnd>iextnd:
        iextnd=i_user_extnd
    gff.load_gene_annotation_with_extnd(iextnd)
    gff.index_gene_annotation_interval_tree()
    gff.annotate_results(sf_input, sf_output)


def call_genotypes(r,options,annot_path_dict,output_dir,rcd,basic_rcd):
    return

def generate_VCF(r,options,annot_path_dict,output_dir):

    sf_prefix = output_dir
    s_rep_type = r
    s_sample_id=options.sample_name

    sf_ref = annot_path_dict['sf_ref'] #reference genome "-ref"

    sf_bam = options.input_bams[0] # used below for vcf header chromosomes
    sf_raw_rslt = f"{sf_prefix}/candidate_disc_filtered_cns.txt.high_confident.post_filtering_with_gene_gntp.txt"

    if os.path.isfile(sf_raw_rslt)==True:
        gvcf=gVCF()
        sf_vcf=f"{sf_prefix}/{s_sample_id}_{s_rep_type}.vcf"
        gvcf.cvt_raw_rslt_to_gvcf(s_sample_id, sf_bam, sf_raw_rslt, s_rep_type, sf_ref, sf_vcf)
    else:
        print(f"Missing file: {sf_raw_rslt}")
