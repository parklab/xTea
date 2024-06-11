##11/27/2017
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu
##@@maintainer: Corinne Sexton, DBMI, Harvard Medical School
##@@contact: corinne_sexton@hms.harvard.edu
##

####
import global_values
from x_TEI_locator import *
from x_local_assembly import *
from x_intermediate_sites import *
from x_reference import *
from x_clip_disc_filter import *
from x_somatic_calling import *
from x_reads_collection import *
from x_mutation import *
from x_gene_annotation import *
from x_genotype_feature import *
from x_basic_info import *
from x_parameter import *

import toml


def parse_toml(toml_file):
    return toml.load(toml_file)

if __name__ == '__main__':
    options = parse_toml("../params/params_prepreocess.toml")

    if options.preprocess:  #preprocess steps
        s_working_folder = options.wfolder
        sf_ref = options.reference
        sf_annotation = options.annotation
        sf_out_fa = options.output
        flank_lth = options.extend
        b_with_flank = options.withflank  # if not set then no flank region

        b_bed_fmt=options.bed
        x_annotation = XAnnotation(sf_annotation)
        i_lextnd=options.rmsk_extnd
        global_values.set_load_rmsk_left_extnd(i_lextnd)
        if b_bed_fmt==True:
            x_annotation.collect_seqs_of_TE_from_ref_bed_fmt(sf_ref, sf_out_fa, flank_lth)
        else:# this is for repeatmasker output
            b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
            x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
            x_annotation.load_rmsk_annotation()
            # x_annotation.collect_flank_regions_of_TE_from_ref(sf_ref, flank_lth, sf_out_fa) #only get the flank regions
            x_annotation.collect_seqs_of_TE_from_ref(sf_ref, flank_lth, b_with_flank, sf_out_fa)
            x_annotation.bwa_index_TE_seqs(sf_out_fa)

    elif options.flank:#preprocess the flank regions steps
        s_working_folder = options.wfolder
        sf_ref = options.reference
        sf_annotation = options.annotation
        sf_out_fa = options.output
        flank_lth = options.extend

        x_annotation = XAnnotation(sf_annotation)
        b_with_chr = x_annotation.is_ref_chrm_with_chr(sf_ref)
        x_annotation.set_with_chr(b_with_chr)  # if chrm in reference has "chr", then True, otherwise False
        b_bed_fmt = options.bed
        if b_bed_fmt==True:
            print("load from bed")
            x_annotation.load_annotation_no_extnd_from_bed()
        else:
            print("load from rmsk")
            x_annotation.load_rmsk_annotation_no_extnd()

        x_annotation.collect_flank_regions_of_TE_from_ref(sf_ref, flank_lth, sf_out_fa)  # only get the flank regions
        x_annotation.bwa_index_TE_seqs(sf_out_fa)