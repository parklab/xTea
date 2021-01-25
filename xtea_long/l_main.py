##11/27/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

'''
#TD-list:
1) wtdbg2 asm step will hang if run in parallel.
Reason possibly because of multiprocess.Popen(), but exact reason is still unknown!
2) Pay attention to pbmm2 output, the cigar is different from normal bam/cram; (solved)
3) need to clean the intermediate files; (solved)
4) Some functions are similar for  "brkpnt collect part" and "clip seq collection part", need to be merged
5) In l_rep_classification.py there is hard code for "b_with_chrm==True"
6) current design by aligning flank regions to assembled contig: will miss those TE insertion promoted deletion cases!!!!!
'''

####
import os
from l_MEI_caller import *
from x_basic_info import *
from x_parameter import *
from l_ghost_TE import *
from l_rep_classification import *
from l_haplotype import *
from l_pseudogene_masker import *
from x_gene_annotation import *
from l_ref_SVA_extractor import *

####
# wfolder="/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/NA12878_10x_v2/promoted_SV/"
# construct_cns_minimap2(ins_chrm, ins_pos, wfolder)
####res
###for short clipped reads and low coverage regions:
# wtdbg2 -i chr6_156041354.fa -l 256 -e 1 -S 1 --rescue-low-cov-edges --node-len 256 --ctg-min-length 256 --ctg-min-nodes 1  -fo chr6_156041354_wtdbg2
####let collect all the whole reads, and then do local assembly, and then align to flank region
####For some cases, it doens't work well, for example 13:86340152 (sometimes because the alignment is not reliable!!!

##parse the options options:
def parse_option():
    parser = OptionParser()

    parser.add_option("-U", "--user",
                      action="store_true", dest="user_specific", default=False,
                      help="User specific cutoff")
    parser.add_option("-C", "--lrd_clip",
                      action="store_true", dest="lrd_clip", default=False,
                      help="Collect the potential clip positions")
    parser.add_option("-A", "--assembly",
                      action="store_true", dest="assembly", default=False,
                      help="Do local assembly for collected long reads")
    parser.add_option("-P", "--Pinpoint",
                      action="store_true", dest="pinpoint", default=False,
                      help="Pinpoint the exact postion and call the insertion seq from assembled seqs")
    parser.add_option("-X", "--complex",
                      action="store_true", dest="complex", default=False,
                      help="For TE insertion promoted SVs")
    parser.add_option("-Y", "--classify",
                      action="store_true", dest="classify", default=False,
                      help="Classify the insertion by alignment to templates")
    parser.add_option("-N", "--polymorphic",
                      action="store_true", dest="polymorphic", default=False,
                      help="Call candidate non reference polymorphic rep copies")
    parser.add_option("-M", "--merge",
                      action="store_true", dest="merge", default=False,
                      help="Merge haplotype results")
    parser.add_option("--pseudo",
                      action="store_true", dest="pseudo", default=False,
                      help="Pseudogene insertion calling")
    parser.add_option("--cmd",
                      action="store_true", dest="asm_cmd", default=False,
                      help="Generate asm command script (for cluster)")
    parser.add_option("--mei_no_asm",
                      action="store_true", dest="mei_no_asm", default=False,
                      help="Call MEI only without asm")
    parser.add_option("--skip_exists",
                      action="store_true", dest="skip_exists", default=False,
                      help="Skip the already assembled ones")
    parser.add_option("--line",
                      action="store_true", dest="call_LINE", default=False,
                      help="Call LINE insertions from rmsk output")
    parser.add_option("--sva",
                      action="store_true", dest="call_SVA", default=False,
                      help="Call SVA insertions from rmsk output")
    parser.add_option("--rsva",
                      action="store_true", dest="rsva", default=False,
                      help="Generate reference SVA copies for samples with given sites")
    parser.add_option("--clean",
                      action="store_true", dest="clean", default=False,
                      help="Remove the intermediate files")
    parser.add_option("-K", "--keep",
                      action="store_true", dest="keep", default=False,
                      help="Keep the intermediate files")
    parser.add_option("--collect_asm",
                      action="store_true", dest="collect_asm", default=False,
                      help="Only collect and asm sequences")
    parser.add_option("--pacbio",
                      action="store_true", dest="pacbio", default=True,
                      help="Pacbio reads")
    parser.add_option("--hg19",#by default it's hg38
                      action="store_true", dest="hg19", default=False,
                      help="Working on reference genome hg19")
    parser.add_option("--cns", dest="consensus",
                      help="repeat consensus file", metavar="FILE")

    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-r", "--ref", dest="reference",
                      help="Reference genome", metavar="FILE") #/or reference consensus
    parser.add_option("--rep", dest="rep_lib",
                      help="repeat folder", metavar="FILE")
    parser.add_option("--gene", dest="gene",
                      help="Gene Annotation file", metavar="FILE")
    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("--i2", dest="input2", default="null",
                      help="input file ", metavar="FILE")
    parser.add_option("--iline", dest="iline",
                      help="LINE-1 results", metavar="FILE")
    parser.add_option("--isva", dest="isva", default="null",
                      help="SVA results", metavar="FILE")
    parser.add_option("-m", "--rmsk", dest="rmsk",
                      help="input file ", metavar="FILE")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")
    parser.add_option("--min", dest="min_copy_len", type="int", default=1000,
                      help="Minimum copy length for collecting polymorphic reads")
    parser.add_option("-s", "--slack", dest="slack", type="int", default=250,
                      help="slack value")
    parser.add_option("-w", "--win", dest="win", type="int", default=75,
                      help="peak window size")
    parser.add_option("-d", "--std", dest="std", type="int", default=60,
                      help="Maximum standard derivation of breakpoints")
    parser.add_option("-c", dest="clip", type="int", default=7,
                      help="cutoff of minimum # of clipped/contained reads")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-y", "--type", dest="type", type="int", default=15,
                      help="Type of repeats working on")
    (options, args) = parser.parse_args()
    return (options, args)

####
def automatic_gnrt_parameters_lrd(sf_bam_list, sf_ref, s_working_folder, n_jobs):
    ####1. collect the basic information
    search_win = 500
    x_basic_info = X_BasicInfo(s_working_folder, n_jobs, sf_ref)
    rcd=x_basic_info.get_cov_is_rlth(sf_bam_list, sf_ref, search_win)
    f_cov=rcd[0]
    rlth=rcd[1]
    mean_is=rcd[2]
    std_var=rcd[3]

    ####2. based on the coverage, set the parameters
    xpar=LRDParameters()
    nclip=xpar.get_par_by_cov(f_cov) #in format (iclip, idisc, i_clip-disc)
    print("Ave coverage is {0}: using parameters clip with value {1}\n".format(f_cov, nclip))
    return nclip, f_cov
####
####
if __name__ == '__main__':
    (options, args) = parse_option()
    l_extd_len = [3500, 1500, 500, -1000]#assembly extended flanking length
    SUFFIX_LEFT_BRKPNT=".left_breakponts"#left-clip reads formed breakpoints
    SUFFIX_RIGHT_BRKPNT = ".right_breakponts"#right-clip reads formed breakpoints
    if options.lrd_clip:#call the candidate sites from clip reads (or internal large insertion positions)
        sf_bam_list = options.bam
        n_jobs = options.cores
        sf_merged_out = options.output
        sf_merged_out_left=sf_merged_out+SUFFIX_LEFT_BRKPNT
        sf_merged_out_right = sf_merged_out + SUFFIX_RIGHT_BRKPNT
        sf_ref = options.reference
        swfolder = options.wfolder
        if swfolder[-1]!="/":
            swfolder+="/"
        #calc coverage
        nclip, fcov=automatic_gnrt_parameters_lrd(sf_bam_list, sf_ref, swfolder, n_jobs)
        if options.user_specific==True:
            nclip=options.clip

        lcaller = L_MEI_Caller(swfolder, n_jobs, sf_ref)
        i_peak_win=options.win
        i_cutoff=int(nclip)
        i_max_cutoff=int(3.2*fcov) #filter out those with large number of clipped reads
        std_dev_cutoff=options.std #standard derivation threshold for clustering clip positions
        lcaller.collect_MEIs_breakpoints(sf_bam_list, sf_ref, i_peak_win, i_cutoff, i_max_cutoff, std_dev_cutoff,
                                         sf_merged_out, sf_merged_out_left, sf_merged_out_right)
####
####
    elif options.assembly:#for given list, assemble the clipped reads and call insertion
        sf_bam_list=options.bam
        sf_sites=options.input
        sf_sites2=options.input2 #another list from a different source
        n_cores=options.cores
        sf_out_fa=options.output
        sf_out_sites=sf_out_fa+".sites"
        sf_ref=options.reference
        swfolder = options.wfolder
        if swfolder[-1]!="/":
            swfolder+="/"
        sf_rep_folder = options.rep_lib  # repeat folder
        i_type = options.type
####
        lcaller=L_MEI_Caller(swfolder, n_cores, sf_ref)
        i_slack=150# will be merged if distance is smaller than this value
        sf_merged_sites=lcaller.merge_sites(sf_sites, sf_sites2, i_slack)

        # negative means collect whole reads, and align abs(size) flanking region
        #l_extd_len = [3500, 2500, 1500, 500, 350]

        #align flank regions to the contig
        #b_asm_cmd=True
        #b_call_seq = False
        #True,False: for asm command only
        #False,True: for call MEI only
        #False, False: for assembly in non-parallel with MEI call (by default)
        b_collect_asm=options.collect_asm #only collect seqs and assemble the sequence
        b_asm_cmd=options.asm_cmd
        b_call_seq=options.mei_no_asm
        b_skip_exist=options.skip_exists
####
        if b_collect_asm==True:
            l_extd_len=[-1000]

        # calc coverage
        nclip, fcov = automatic_gnrt_parameters_lrd(sf_bam_list, sf_ref, swfolder, n_cores)
        if options.user_specific == True:
            nclip = options.clip
        n_cutoff=nclip#int(nclip*0.75)#set a cutoff for number of qualified reads (has kmer support)
        n_min_non_polyA_kmer=21
        n_min_polyA_kmer=3#
        print("clip cutoff is: "+str(nclip))

        #get the consensus file list, for the repeat kmer library
        l_sf_cns=lcaller.set_rep_lib_path_kmer_cnt(i_type, sf_rep_folder)
        sf_script=sf_out_fa+".run_asm_cmd"
        # lcaller.call_MEIs_for_sites_with_cmd_version(sf_bam_list, sf_merged_sites, sf_ref, l_extd_len, sf_out_fa,
        #                                              sf_out_sites, b_collect_asm, b_asm_cmd, sf_script,
        #                                              b_call_seq, b_skip_exist)
        lcaller.call_MEIs_for_sites_with_rep_kmer_filtering(sf_bam_list, sf_merged_sites, sf_ref, l_extd_len, sf_out_fa,
                                                            sf_out_sites, l_sf_cns, n_cutoff, n_min_non_polyA_kmer, n_min_polyA_kmer,
                                                            b_collect_asm, b_asm_cmd, sf_script, b_call_seq, b_skip_exist)
        with open(sf_script+"_all_asm.sh", "w") as fmerged:
            for i_extnd in l_extd_len:
                sf_tmp_cmd=sf_script+"_{0}.sh".format(i_extnd)
                if os.path.isfile(sf_tmp_cmd)==True:
                    with open(sf_tmp_cmd) as fin_tmp:
                        for line in fin_tmp:
                            fmerged.write(line)
        lcaller.clean_intermediate_files2(l_extd_len, sf_merged_sites)##
        #align the contig to the target region (this is not used), this version not used for now
        #####lcaller.call_MEIs_for_sites_with_contig_2_ref(sf_bam_list, sf_sites, sf_ref, l_extd_len, sf_out_fa, sf_out_sites)
####
####for complex events
    elif options.complex:#call TE insertion promoted SV from assembled left and right contigs
        #1. for each site, align the flanking region of itself to the contig, and save the seq of the other side to file
        #2.1 align the left flanks to the newly saved right contigs
        #2.2 align the right flanks to the newly saved left contigs
        #3. If two breakpoints support each other, then call it out
        #4. Based on the orientation and chrms
        sf_lsites = options.input#left breakpoint sites
        sf_rsites = options.input2  # right breakpoint sites
        n_cores = options.cores
        sf_out_prefix = options.output
        sf_out_sites = sf_out_prefix + ".sites"
        sf_ref = options.reference
        swfolder = options.wfolder
        sf_rep_folder=options.rep_lib#repeat library folder
        b_hg19 = options.hg19
        sf_l1_rslt=options.iline
        sf_sva_rslt=options.isva
        if swfolder[-1] != "/":
            swfolder += "/"
        f_flk_map_ratio=0.75
        i_extd_len=-1000
        flk_lenth=3000
        global_values.set_lrd_extnd_len(i_extd_len)
        lcaller = L_MEI_Caller(swfolder, n_cores, sf_ref)
        lcaller.call_complex_MEIs_from_asm_in_parallel(sf_ref, sf_lsites, sf_rsites, swfolder, i_extd_len,
                                                       f_flk_map_ratio, sf_l1_rslt, sf_sva_rslt, flk_lenth,
                                                       sf_rep_folder, b_hg19, sf_out_prefix)
####
    elif options.classify:#classify the insertions to different types
        sf_rep_ins=options.input #this is the input seq in fasta format
        sf_rslt=options.output #labeled output
        sf_rep_folder=options.rep_lib#repeat folder
        sf_ref=options.reference
        n_jobs = options.cores
        swfolder = options.wfolder
        i_type=options.type
        b_hg19=options.hg19
        flk_lenth=3000
        if swfolder[-1] != "/":
            swfolder += "/"
        lrc=LRepClassification(swfolder, n_jobs)
        lrc.set_rep_configuration(i_type, sf_rep_folder, b_hg19)
        lrc.classify_ins_seqs(sf_rep_ins, sf_ref, flk_lenth, sf_rslt)
####
    elif options.merge:#classify the insertions to different types
        sf_prefix1=options.input #this is the input seq in fasta format
        sf_prefix2=options.input2
        sf_out_prefix=options.output #labeled output
        n_jobs = options.cores
        swfolder = options.wfolder
        i_type=options.type
        i_slack=50
        if swfolder[-1] != "/":
            swfolder += "/"
        lhap=IHaplotype()
        lhap.merge_xTEA_output_of_two_hap(sf_prefix1, sf_prefix2, i_type, swfolder, n_jobs, i_slack, sf_out_prefix)
####
    #find ghost full length copies
    elif options.polymorphic:#call polymorphic copies
        sf_bam_list = options.bam
        sf_rmsk_copy = options.rmsk
        i_min_rep_len=options.min_copy_len
        sf_rep_cns = options.consensus  #repeat consensus
        n_jobs = options.cores
        sf_out_fa = options.output
        swfolder = options.wfolder

        if swfolder[-1] != "/":
            swfolder += "/"
        sf_ref = options.reference
        i_slack=options.slack
        global_values.set_polymorphic_brk_chk_win(i_slack)
        lnrp=LNonRefPolymorphic(n_jobs, swfolder)
        lnrp.collect_ghost_polymorphic_rep_reads(sf_rmsk_copy, i_min_rep_len, sf_bam_list, sf_ref, sf_out_fa)

        s_cluster_folder=swfolder+"cluster/"
        b_pacbio=options.pacbio
        sf_flank=sf_out_fa + ".separate_flanking.fa"
        i_max_clip=200 #clip part smaller than 200bp
        iset_cutoff=5 #at least 5 reads in the set
        sf_algnmt=sf_flank+".algn_2_itself.sorted.bam"
        i_min_overlap=3500 #at least 5k overlap
        l_cluster_sfa_list=lnrp.cluster_reads_by_flank_region(sf_algnmt, sf_out_fa, sf_flank, b_pacbio,
                                  i_max_clip, i_min_overlap, iset_cutoff, s_cluster_folder)
        ####
        #asm the cluster to construct the copy and flanking regions
        sf_merged_copies=sf_out_fa+".assembled_final_out.fa"
        lnrp.asm_collect_cluster_reads(l_cluster_sfa_list, sf_merged_copies)
        ####
        #align the consensus to the asm to call out the full length copy, and seperate the flanking regions
        s_sample_id="ghost_polymorphic"
        sf_prefix=os.path.dirname(sf_out_fa)
        lnrp.seprt_cns_flank_of_asm_contig(s_sample_id, sf_merged_copies, sf_rep_cns, sf_ref, sf_prefix)

        #align the flank regions to reference genome
        #collect those unmapped or mapped to unknown contigs
    elif options.rsva:
        sf_bam_list = options.bam
        sf_sites = options.input #here each copy in format: chrm start end
        n_cores = options.cores
        sf_out_fa = options.output
        sf_ref = options.reference
        swfolder = options.wfolder
        if swfolder[-1] != "/":
            swfolder += "/"

        ref_sva=RefSVAExtractor(swfolder, n_cores)
        b_asm_cmd=False
        ref_sva.collect_asm_extract_ref_sva_seq(sf_sites, sf_ref, sf_bam_list, b_asm_cmd, sf_out_fa)

    elif options.pseudo:
        n_jobs = options.cores
        swfolder = options.wfolder
        sf_out = options.output
        sf_exon=options.reference
        sf_contigs=options.input

        if swfolder[-1] != "/":
            swfolder += "/"
        pseudo=PsudogeneMasker(swfolder, n_jobs)
        sf_algnmt=swfolder+"exon_2_contigs.sorted.bam"
        pseudo.align_exon_2_contigs(sf_contigs, sf_exon, sf_algnmt)
        f_max_clip_cutoff=0.15
        f_min_cov=0.9#at least 90% of the contig is covered
        sf_out_tmp=sf_out+".tmp"
        pseudo.parse_exon_2_contig_algnmt(sf_algnmt, sf_contigs, f_max_clip_cutoff, f_min_cov, sf_out_tmp)

        sf_gene_annotation = options.gene
        gff = GFF3(sf_gene_annotation)
        iextnd = global_values.UP_DOWN_GENE
        gff.load_gene_annotation_with_extnd(iextnd)
        gff.index_gene_annotation_interval_tree()
        gff.annotate_results(sf_out_tmp, sf_out)

####
####call from rmsk output
    elif options.call_LINE:
        sf_rmsk = options.rmsk
        sf_sites=options.input
        sf_out = options.output
        cfrmsk=CallFromRMSK()
        cfrmsk.call_TEI_from_rmsk_L1(sf_rmsk, sf_sites, sf_out)
####
    elif options.call_SVA:
        sf_rmsk = options.rmsk
        sf_sites = options.input
        sf_out = options.output
        cfrmsk = CallFromRMSK()
        cfrmsk.call_TEI_from_rmsk_SVA(sf_rmsk, sf_sites, sf_out)

    elif options.clean:#clean the un-necessary files
        sf_sites = options.input
        swfolder = options.wfolder #this is set same as "options.assembly"
        n_jobs = options.cores
        sf_ref=options.reference
        lcaller = L_MEI_Caller(swfolder, n_jobs, sf_ref)
        lcaller.clean_intermediate_files2(l_extd_len, sf_sites)

####
####
#sf_bam = "/n/data1/hms/dbmi/park/simon_chu/projects/data/NA12878_pacbio/NA12878_pacbio_40x_hg19.sorted.bam"

# sf_bam="/n/data1/hms/dbmi/park/simon_chu/projects/data/NA12878_nanopore/ngmlr_alignment/ngm_Nanopore_human_ngmlr-0.2.3_mapped.bam"
# sf_bam="/n/data1/hms/dbmi/park/simon_chu/projects/data/NA12878_pacbio/NA12878_pacbio_40x.sorted.bam"
# sf_bam="/n/data1/hms/dbmi/park/simon_chu/projects/data/1000G_trio/Pacbio/WGS/alignment_hg38/NA19240_2_ref.sorted.bam"

####
####