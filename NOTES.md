# 23 main functions in x_TEA_main.py

## global options set on cmdline
--mit 
--dna
--cbs
--sva
--user_specific
--tumor
--purity
--resume


## Main params
- -P **preprocess**
    - preprocess steps
    wfolder, reference, annotation, output
    extend, withflank, bed, rmsk_extnd
- **flank**
    - preprocess the flank regions steps
    wfolder, reference, annotation, output,
    extend, bed
- -C **clip**
    - take in the normal illumina reads (10x will be viewed as normal illumina)
    input, wfolder, cores, cns, reference,
    annotation, output, single, ref, force
    mosaic, siteclip, lclip, rclip, cliprep,
    cov, cwfolder
    - functions:
    TE_Multi_Locator::call_TEI_candidate_sites_from_multiple_alignmts
        TELocator::call_TEI_candidate_sites_from_clip_reads_v2
- -D **discordant**
    - this views all the alignments as normal illumina reads
- -N **filter_csn** 
    - filter out the FP by the pattern in the consensus repeat
- --transduction **transduction**
    - need to re-collect all the clip, disc reads
- --sibling **sibling**
    - sibling orphan transduction
- --case_control **case_control**
    - case-control mode to call somatic events
- **mosaic**
    - this is only for normal illumina data 
- -B **barcode**
    - this is only for 10X alignmt
- --postF **postF**
    - post filtering step
- --gntp_feature **gntp_feature**
    - generate the genotype features
- --gntp_classify **gntp_classify**
    - predict the genotype
    --train_gntp
    --model
    --input
    --output
- **gVCF**
    - 
- -J **joint**
    - MosaicJointCalling
    - call_mosaic_from_multi_samples
- --igv **igv**
    - prepare the igv screenshot script for multiple bams and sites
- -E **collect**
    - collect the reads for each candidate site
- **mutation**
    - call out the internal mutations by aligning the reads
- **gene**
- **map**
    - align the asm to reference genome
- **flk_map**
    - gnrt the flank regions and align the regions to the contigs
- -F **filter_asm**
- -Q **collect_clip**
    - collect the clipped reads for the sample


## extra params
--resume
--spectrum
-I, --mutation
-U, --collect_Illumina
-G, --contig_realign
-T, --trace
-A, --assembly
-L, --local
-M, --map
-V, --visualization
-K, --withflank
--mit
--dna
--cbs
--sva
