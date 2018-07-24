##11/27/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

'''
####@@CC 01/03/2018
Todo List: 1. clip step requires lots of memory. 150G for 10 cores. This should be improved.
Update: 01/09/2018: The issue is caused at the "ClipReadInfo::cnt_clip_part_aligned_to_rep()" function.
Solution:
        1)Skip the load all into memory step, instead only load those has "left(right) realigned" ones. This is
            the current version used.
        2)For the re-aligned sam file, parse out the chrm, positions, and write to a file, and sort the file,
          then, count the number of supported reads. This version should be more memory efficient, but hasn't been
          fully tested.
####

'''

import os
from x_TEI_locator import *
#from x_TEI_source_tracer import *
from x_local_assembly import *
from x_filter import *
from x_reference import *
from x_clip_disc_filter import *
from x_mosaic_calling import *
from x_analysis import *
from optparse import OptionParser
from x_reads_collection import *
from x_mutation import *
from x_sv import *
####
##parse the options
def parse_option():
    parser = OptionParser()
    parser.add_option("-P", "--preprocess",
                      action="store_true", dest="preprocess", default=False,
                      help="Preprocessing stpes")
    parser.add_option("-C", "--clip",
                      action="store_true", dest="clip", default=False,
                      help="Call candidate TEI sites from clipped reads")
    parser.add_option("-S", "--single",
                      action="store_true", dest="single", default=False,
                      help="Call clip positions from single-end reads")
    parser.add_option("-D", "--discordant",
                      action="store_true", dest="discordant", default=False,
                      help="Filter with discordant paired end reads")
    parser.add_option("-N", "--filter_csn",
                      action="store_true", dest="filter_csn", default=False,
                      help="Filter out candidate sites from map position on consensus")
    parser.add_option("-B", "--barcode",
                      action="store_true", dest="barcode", default=False,
                      help="Indicate the input is 10X bam")
    parser.add_option("-E", "--collect",
                      action="store_true", dest="collect", default=False,
                      help="Collect reads for candidate sites")
    parser.add_option("-I", "--mutation",
                      action="store_true", dest="mutation", default=False,
                      help="Call internal mutation")
    parser.add_option("-U", "--collect_Illumina",
                      action="store_true", dest="collect_illumina", default=False,
                      help="Collect reads for candidate sites from normal illumina alignment")
    parser.add_option("-F", "--filter_asm",
                      action="store_true", dest="filter_asm", default=False,
                      help="Filter out candidate sites from assembly")
    parser.add_option("-G", "--contig_realign",
                      action="store_true", dest="contig_realign", default=False,
                      help="Filter out candidate sites from assembly")
    parser.add_option("-T", "--trace",
                      action="store_true", dest="trace", default=False,
                      help="Trace the sources of TEIs")
    parser.add_option("-A", "--assembly",
                      action="store_true", dest="assembly", default=False,
                      help="Do local assembly for collected reads")
    parser.add_option("-L", "--local",
                      action="store_true", dest="local", default=False,
                      help="Assemble the TEIs on local machine")
    parser.add_option("-M", "--map",
                      action="store_true", dest="map", default=False,
                      help="map flank regions to the assembled contigs")
    parser.add_option("-V", "--visualization",
                      action="store_true", dest="visualization", default=False,
                      help="Show the heatmap figure of the selected regions")
    parser.add_option("-K", "--withflank",
                      action="store_true", dest="withflank", default=False,
                      help="Keep the flank regions with the repeat copies")
    parser.add_option("--bed",
                      action="store_true", dest="bed", default=False,
                      help="Input annotation in bed format")
    parser.add_option("--mosaic",
                      action="store_true", dest="mosaic", default=False,
                      help="Call mosaic events")
    parser.add_option("--flk_map",
                      action="store_true", dest="flk_map", default=False,
                      help="Map flanks to contigs")
    parser.add_option("--analysis",
                      action="store_true", dest="analysis", default=False,
                      help="Result analysis")
    parser.add_option("--flank", dest="flank", default=False,
                      help="flank regions")
    parser.add_option("--sv",
                      action="store_true", dest="sv", default=False,
                      help="Call promoted SVs")

    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("-r", "--reference", dest="reference",
                      help="The reference file ", metavar="FILE")
    parser.add_option("-a", "--annotation", dest="annotation",
                      help="The annotation file ", metavar="FILE")
    # parser.add_option("-c", "--copies", dest="copies",
    #                   help="Repeat copies ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-d", "--barcode_bam", dest="barcode_bam",
                      help="Input barcode indexed bam file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")
    parser.add_option("-e", "--extend", dest="extend", type="int",
                      help="extend length")
    parser.add_option("-u", "--dup", dest="duplication",
                      help="duplication files", metavar="FILE")
    parser.add_option("--fflank", dest="fflank",
                      help="flank region file", metavar="FILE")
    parser.add_option("--flklen", dest="flklen", type="int",
                      help="flank region file")
    parser.add_option("--ref", dest="ref",
                      help="repeat consensus/reference", metavar="FILE")
    parser.add_option("--lc", dest="lclip", type="int",
                      help="cutoff of minimum # of left clipped reads")
    parser.add_option("--rc", dest="rclip", type="int",
                      help="cutoff of minimum # of rightt clipped reads")
    parser.add_option("--cr", dest="cliprep", type="int",
                      help="cutoff of minimum # of clipped parts fall in repeats")
    parser.add_option("--nd", dest="ndisc", type="int",
                      help="cutoff of minimum # of discordant pair")
    parser.add_option("--nb", dest="nbarcode", type="int",
                      help="cutoff of maximum # of molecure coverage")
    parser.add_option("--teilen", dest="teilen", type="int",
                      help="minimum length of the insertion for future analysis")

    (options, args) = parser.parse_args()
    return (options, args)

####
##main function
if __name__ == '__main__':
    (options, args) = parse_option()

    if options.preprocess:  # preprocess steps
        s_working_folder = options.wfolder
        sf_ref = options.reference
        sf_annotation = options.annotation
        sf_out_fa = options.output
        flank_lth = options.extend
        b_with_flank = options.withflank  # if not set then no flank region

        b_bed_fmt=options.bed
        x_annotation = XAnnotation(sf_annotation)
        if b_bed_fmt==True:
            x_annotation.collect_seqs_of_TE_from_ref_bed_fmt(sf_ref, sf_out_fa)
        else:# this is for repeatmasker output
            b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
            x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
            x_annotation.load_rmsk_annotation()
            # x_annotation.collect_flank_regions_of_TE_from_ref(sf_ref, flank_lth, sf_out_fa) #only get the flank regions
            x_annotation.collect_seqs_of_TE_from_ref(sf_ref, flank_lth, b_with_flank, sf_out_fa)
            x_annotation.bwa_index_TE_seqs(sf_out_fa)
    elif options.flank:  # preprocess the flank regions steps
        s_working_folder = options.wfolder
        sf_ref = options.reference
        sf_annotation = options.annotation
        sf_out_fa = options.output
        flank_lth = options.extend

        x_annotation = XAnnotation(sf_annotation)
        b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
        x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
        x_annotation.load_rmsk_annotation_no_extnd()

        x_annotation.collect_flank_regions_of_TE_from_ref(sf_ref, flank_lth, sf_out_fa)  # only get the flank regions
        x_annotation.bwa_index_TE_seqs(sf_out_fa)

    elif options.clip:  ###take in the normal illumina reads (10x will be viewed as normal illumina)
        sf_bam_list = options.input
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_rep = options.reference  ####repeat copies "-r"
        sf_annotation = options.annotation
        sf_out = options.output
        b_se = options.single  ##single end reads or not, default is not
        sf_ref=options.ref ###reference genome "-ref"

        # merge the list from different bams of the same individual
        # Here when do the filtering, nearby regions are already considered!
        cutoff_left_clip = options.lclip
        cutoff_right_clip = options.rclip
        cutoff_clip_mate_in_rep = options.cliprep
        tem_locator = TE_Multi_Locator(sf_bam_list, s_working_folder, n_jobs, sf_ref)

        ##Hard code inside:
        # 1. call_TEI_candidate_sites_from_clip_reads_v2 --> run_cnt_clip_part_aligned_to_rep_by_chrm_sort_version
        # here if half of the seq is mapped, then consider it as aligned work.
        tem_locator.call_TEI_candidate_sites_from_multiple_alignmts(sf_annotation, sf_rep, b_se, cutoff_left_clip,
                                                                    cutoff_right_clip, cutoff_clip_mate_in_rep, sf_out)

    elif options.discordant:  # this views all the alignments as normal illumina reads
        sf_list = options.bam  ###read in a bam list file
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_annotation = options.annotation
        sf_candidate_list = options.input
        sf_out = options.output
        sf_ref = options.ref  ###reference genome, some cram file require this file to open
        PEAK_WINDOW = 100

        xfilter = XFilter()
        m_original_sites = xfilter.load_in_candidate_list(sf_candidate_list)
        sf_peak_sites = s_working_folder + "clip_peak_candidate.list"
        m_sites_clip_peak = xfilter.call_peak_candidate_sites(m_original_sites, PEAK_WINDOW)  # get the peak sites
        xfilter.output_candidate_sites(m_sites_clip_peak, sf_peak_sites)  # output the sites
        m_original_sites.clear()  # release the memory

        iextend = 800
        i_is = 100000  ###set the insert size a large value
        f_dev = 50

        # this is the cutoff for  "left discordant" and "right discordant"
        # Either of them is larger than this cutoff, the site will be reported
        n_disc_cutoff = options.ndisc

        sf_tmp = s_working_folder + "disc_tmp.list"
        tem_locator = TE_Multi_Locator(sf_list, s_working_folder, n_jobs, sf_ref)
        tem_locator.filter_candidate_sites_by_discordant_pairs_multi_alignmts(m_sites_clip_peak, iextend, i_is, f_dev,
                                                                              n_disc_cutoff, sf_annotation, sf_tmp)
        xfilter.merge_clip_disc(sf_tmp, sf_candidate_list, sf_out)
    ####
    elif options.filter_csn:  # filter out the FP by the pattern in the consensus repeat.
        sf_bam_list = options.bam  ###read in a bam list file
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_ref = options.ref  ###reference genome, some cram file require this file to open

        sf_candidate_list = options.input
        iextnd = 400  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
        bin_size = 50000000  # block size for parallelization
        sf_cns = options.reference  ####repeat copies/cns here
        bmapped_cutoff = 0.5
        sf_annotation = options.annotation
        i_concord_dist = 550  # this should be the 3*std_derivation, used to cluster disc reads on the consensus
        f_concord_ratio = 0.45
        n_clip_cutoff = options.cliprep# this is the sum of left and right clipped reads
        n_disc_cutoff = options.ndisc #each sample should have at least this number of discordant reads
        sf_output = options.output
        sf_flank=options.fflank
        i_flank_lenth = options.flklen

        x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_ref)
        x_cd_filter.call_MEIs_consensus(sf_candidate_list, iextnd, bin_size, sf_cns, sf_flank, i_flank_lenth,
                                        bmapped_cutoff, i_concord_dist, f_concord_ratio, n_clip_cutoff, n_disc_cutoff,
                                        sf_output)

    elif options.mosaic:  # this is only for normal illumina data
        sf_bam_list = options.bam
        sf_candidate_list = options.input
        n_jobs = options.cores
        sf_black_list = options.annotation
        sf_1KG = options.reference
        sf_dup_list = options.duplication
        s_working_folder = options.wfolder
        sf_output = options.output
        sf_ref = options.ref  ###reference genome, some cram file require this file to open
        extnd = 150
        clip_slack = 25
        af_cutoff = 0.1

        mMEI = MosaicMEICaller(sf_ref)
        sf_af_filtered = sf_output + ".af_filtered"
        mMEI.filter_by_AF(sf_bam_list, sf_candidate_list, extnd, clip_slack, s_working_folder, n_jobs, af_cutoff,
                          sf_af_filtered)

        sf_blacklist_filtered = sf_output + ".black_list_filtered"
        mMEI.filter_by_black_list(sf_af_filtered, sf_black_list, sf_blacklist_filtered)

        sf_duplication_filtered = sf_output + ".dup_filtered"
        mMEI.filter_by_black_list(sf_blacklist_filtered, sf_dup_list, sf_duplication_filtered)

        nearby_slack = 300
        mMEI.filter_by_germline_list(sf_1KG, sf_blacklist_filtered, nearby_slack, sf_output)

    ####
    elif options.barcode:  # this is only for 10X alignmt
        sf_ori_bam = options.bam  # pos indexted bam
        b_barcode = options.barcode
        sf_barcode_bam = options.barcode_bam  # barcode indexed bam
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_annotation = options.annotation
        sf_candidate_list = options.input
        sf_out = options.output
        sf_ref = options.ref  ###reference genome, some cram file require this file to open

        xfilter = XFilter()
        m_input_sites = xfilter.load_in_candidate_list(sf_candidate_list)

        iextend = 2500
        i_cov_cutoff = options.nbarcode

        if s_working_folder[-1] != "/":
            s_working_folder += "/"
        sf_tmp = s_working_folder + "barcode_tmp.list"
        caller = TELocator(sf_ori_bam, sf_barcode_bam, s_working_folder, n_jobs, sf_ref)
        m_sites_barcode = caller.filter_candidate_sites_by_barcode_coverage(m_input_sites, iextend, i_cov_cutoff)
        caller.output_candidate_sites(m_sites_barcode, sf_tmp)  # clip_peak_discord_candidate_barcode.list
        sf_tmp2 = s_working_folder + "barcode_tmp2.list"
        xfilter.merge_clip_disc_barcode(sf_tmp, sf_candidate_list, sf_tmp2)

        ###combine the ones are close to each other
        window_size = 500  # if closer than 500, then the sites are merged as one
        xfilter.combine_closing_sites(sf_tmp2, window_size, sf_out)

    elif options.collect:  # collect the reads for each candidate site
        sf_ori_bam = options.bam
        sf_barcode_bam = options.barcode_bam
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_candidate_list = options.input
        i_cov_cutoff = options.nbarcode
        sf_annotation = options.annotation
        i_flank_lenth = options.flklen #extend the boundary of each annotated repeat
        sf_ref = options.ref  ###reference genome, some cram file require this file to open

        # collect reads for all the sites
        i_extend = 1500 #by default, collect barcode in [-1500, 1500] region
        xread_collection=XReadsCollection(s_working_folder, sf_ref)
        xread_collection.collect_phased_reads_all_TEIs(sf_ori_bam, sf_barcode_bam, sf_candidate_list, i_extend,
                                                       i_cov_cutoff, sf_annotation, i_flank_lenth, n_jobs)

    elif options.mutation:#call out the internal mutations by aligning the reads
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_cns = options.reference  ####repeat copies/cns here
        sf_sites=options.input
        sf_merged_vcf=options.output
        n_len_cutoff=options.teilen
        xmutation=XMutation(s_working_folder)
        xmutation.call_mutations_from_reads_algnmt(sf_sites, sf_cns, n_len_cutoff, n_jobs, sf_merged_vcf)
####

    # this step needs to retrieve the seqs of mate reads
    elif options.collect_illumina:  # this is only for normal illumina data
        sf_ori_bam = options.bam
        sf_candidate_list = options.input
        n_jobs = options.cores
        s_working_folder = options.wfolder
    ####

    elif options.assembly:  # assemble the reads for all the sites
        b_local = options.local
        s_working_folder = options.wfolder
        n_jobs = options.cores
        sf_candidate_list = options.input
        sf_ref = options.ref  ###reference genome, some cram file require this file to open

        xlasm = XLocalAssembly(s_working_folder, sf_ref)
        if b_local == True:
            # xlasm.assemble_all_TEIs_locally(sf_candidate_list, n_jobs)
            xlasm.assemble_all_phased_TEIs_locally(sf_candidate_list, n_jobs)
        else:
            # xlasm.assemble_all_TEIs_slurm_cluster(sf_candidate_list)
            xlasm.assemble_all_phased_TEIs_slurm_cluster(sf_candidate_list)

    elif options.map:####
        sf_ref = options.reference
        s_working_folder = options.wfolder
        n_jobs = int(options.cores)
        sf_sites = options.input
        sf_cns=options.ref####repeat consensus
        sf_final_list=options.output
        sf_final_seqs=sf_final_list+".fa"

        xctg = XTEContig(s_working_folder, n_jobs)
        #xctg.align_asm_contigs_to_reference(sf_sites, sf_ref, s_working_folder)
        xctg.call_MEIs_from_all_group_contigs(sf_sites, sf_ref, sf_cns, s_working_folder, sf_final_list, sf_final_seqs)

    elif options.flk_map:  #gnrt the flank regions and align the regions to the contigs
        sf_ref = options.reference
        s_working_folder = options.wfolder
        n_jobs = int(options.cores)
        sf_sites = options.input

        i_extend = 500
        xref = XReference()
        xref.gnrt_flank_region_for_sites(sf_sites, i_extend, n_jobs, sf_ref, s_working_folder)
        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.align_flanks_to_phased_contig(sf_sites)

    elif options.filter_asm:
        s_working_folder = options.wfolder
        n_jobs = int(options.cores)
        sf_sites = options.input
        sf_keep_sites = options.output
        sf_repeat_copies = options.reference  ####repeat copies here

        # xctg.filter_out_non_TE_from_asm(sf_sites, flank_length, f_map_cutoff, i_slack, sf_keep_sites)
        flank_length = 500
        f_map_cutoff = 0.45 #25% of the bases are required to be mapped! This is a small ratio
        i_slack = 35
        xctg = XTEContig(s_working_folder, n_jobs)
        xctg.validate_TEI_from_phased_asm_algnmt(sf_sites, flank_length, f_map_cutoff, i_slack, sf_repeat_copies,
                                                 sf_keep_sites)
    elif options.sv:
        sf_sites = options.input
        sf_dels=options.output
        te_del=TEDeletion()
        m_del=te_del.call_del_matched_brkpnts(sf_sites)
        te_del.output_del(m_del, sf_dels)

####
    # elif options.trace:
    #     sf_vcf = options.input
    #     sf_bam_coordinate = options.bam
    #     sf_barcode_bam = options.barcode_bam
    #     s_working_folder = options.wfolder
    #     sf_annotation = options.annotation
    #     n_jobs = options.cores
    #     sf_ref = options.ref  ###reference genome, some cram file require this file to open
    #
    #     st = SourceTracer(sf_vcf, sf_bam_coordinate, sf_barcode_bam, sf_annotation, n_jobs, s_working_folder, sf_ref)
    #     st.load_files()
    #     g_network = st.trace_source_for_one_site("chr5", 89450781, 500)

    elif options.visualization:  ##show the shared barcode heatmap between the TEI site and source region
        sf_matrix = options.input

        ####classify the TE insertions to normal ones and complex ones:
        # for insertion with deletion: the other clipped part will not be aligned to a repeat region, but the other breakpoint
        # while for normal ones, the other clipped part will be aligned to a repeat region

        ####1. python

        ####2. python x_TEA_main.py -T -i ${VCF} -r ${ANNOTATION} -b ${BAM} -d ${BARCODE_BAM} -p ${TMP} -n 1
        ####3. python

        ####4.
    elif options.analysis:
        sf_candidate_list=options.input
        sf_annotation = options.annotation
        sf_csv_prefix=options.output
        x_analysis=XAnalysis()
        #cnt # of TE insertions fall in repetitive regions
        x_analysis.cnt_fall_in_rep(sf_annotation, sf_candidate_list, sf_csv_prefix)
####
####