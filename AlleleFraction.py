import pysam
import pathlib
import pandas as pd
import numpy as np

import utils
import MeerkatPhaser as phaser

def AF_from_rslt_list(rslt_list, f_bam, map_min=75, clip_min=13):
    bam = pysam.AlignmentFile(f_bam)

    AF_list = []
    for f_rslt in rslt_list:
        l = AF_from_SVmate(f_rslt, bam, map_min=map_min, clip_min=clip_min)
        AF_list.extend(l)

    df = pd.DataFrame(AF_list, columns=['fusion', 'chr_fuse', 'breakpoint', 'orient', 'n_ref', 'n_alt', 'n_clip', 'n_disc', 'AF'])
    df2 = df[~df[['chr_fuse', 'breakpoint']].duplicated()]

    return df2

def AF_from_SVmate(f_rslt, bam, clip_min=13, map_min=75):
    """ Calculate allele fraction for an SVmate rslt file
        Should provide an AF for upstream and downstream side of each binary fusion
    """
    df = pd.read_table(f_rslt, names=range(10))

    # for key, group in df.groupby(3):
    #     if len(group[group[1] == 

    AF_list = []
    for key, group in df.groupby([4, 5]):
        fusion = group.iloc[0, 3]
        orient = group.iloc[0, 8]

        g_clip = group.loc[(group[1] == 'clip') & (group[2] == orient), :]
        g_disc = group.loc[(group[1] == 'disc') & (group[2] == orient), :]

        bp = None
        if len(g_clip):
            discord_only = False

            if len(set(g_clip[2])) > 1:
                print("WARNING: there is more than one clip direction present for this fusion! Using only reads with {} orientation.".format(orient))

            if len(g_clip) <= 2:
                bp = g_clip[6].tolist()[0]
            else:
                try:
                    bp = g_clip.loc[g_clip[2] == orient, 6].mode()[0]
                except:
                    print("WARNING: no clipped read consensus for {} fusion {}. Using discordant reads instead.".format(orient, fusion))
                    # raise ValueError(print(g_clip))

        # elif len(g_disc):
        if not bp:
            discord_only = True

            if len(set(g_disc[2])) > 1:
                print("WARNING: there is more than one discordant direction present for this fusion! Using only reads with {} orientation.".format(orient))

            if orient == 'downstream':
                bp = g_disc[6].max() + 150  # Hard-coded read length for now 

            else:
                bp = g_disc[6].min()

        # else:
        if not bp:
            print("WARNING: no {} reads support the fusions {}".format(orient, fusion))
            # print(df)
            continue

        reads_alt = df.loc[df[3] == fusion, 0]
        chrom, _, _ = utils.parse_loc_str(group.iloc[0, 4])

        if np.isnan(bp):
            # print(df)
            continue

        print(bp)

        n_disc = len(set(df.loc[(df[3] == fusion) & (df[1]=='disc'), 0]))
        n_clip = len(set(df.loc[(df[3] == fusion) & (df[1]=='clip'), 0]))

        n_ref, n_alt, AF, rl = AF_by_bp(chrom, bp, orient, reads_alt, bam, map_min=map_min, clip_min=clip_min, discord_only=discord_only)
        # n_ref, n_alt, AF, rl = AF_by_bp(chrom, bp, orient, reads_alt, bam, overlap_min, discord_only)

        AF_list.append([fusion, chrom, bp, orient, n_ref, n_alt, n_clip, n_disc, AF])
        # AF_list.append([fusion, chrom, bp, n_ref, n_alt, AF, rl])

    return AF_list

# def AF_by_bp(chrom, bp, orientation, reads_alt, bam, overlap_min=13, discord_only=False):
def AF_by_bp(chrom, bp, orientation, reads_alt, bam, map_min=13, clip_min=13, discord_only=False):
    """ Calculate allele fraction for an event with a known breakpoint and alternate reads

        Inputs:
            chrom: chromosome of event (as string)
            bp: positon of breakpoint
            orientation: "downstream" or "upstream" to indicate whether reads are left-clipped / forward or
                         right clipped / reverse
            reads_alt: list of alternate read names
            bam: pysam bam file opened for reading
    """
    if orientation == 'downstream':
        start = bp - 1000
        end   = bp + 50
        # end   = bp - clip_min
    else:
        # start = bp + clip_min
        start = bp - 50
        end = bp + 1000

    read_list = []
    for read in bam.fetch(chrom, start, end):
        if read.query_name in read_list:
            continue

        if not phaser.read_PASS(read):
            continue

        # Adjust start / end for of read and mate for any softclipped bases
        d_read = phaser._read_bounds(read, as_dict=True)
        d_mate = phaser._mate_bounds(read, as_dict=True)
        pair_dir = _pair_direction(read)
        # is_consistent = _check_read_direction(read, orientation)
        # orient = _pair_orientation(read)
        # insert = _insert_size(qstart, qend, mate_qstart, mate_qend)

        # Mapping qualities
        mapq = read.mapq
        mate_mapq = phaser._mate_mapq(read)

        # Is read purely on the reference side of the bp --> look at mate
        if _read_ref_to_bp(d_read['qstart'], d_read['qend'], bp, orientation):
            # Is mate beyond bp?
            if _mate_beyond_bp(d_read, d_mate, pair_dir, bp, orientation, discord_only):
                read_list.append(read.query_name)
                continue

            # Does mate overlap bp.
            # if phaser._read_overlap_bp2(d_mate['qstart'], d_mate['qend'], bp, orientation, map_min=map_min, clip_min=clip_min):
            # if phaser._read_overlap_bp2(d_mate['start'], d_mate['end'], bp, orientation, map_min=map_min, clip_min=clip_min):
            if phaser._read_overlap_bp2(d_mate['start'], d_mate['end'], bp, orientation, map_min=1, clip_min=clip_min):
            # if phaser._read_overlap_bp(d_mate['qstart'], d_mate['qend'], bp, clip_min):
                if not discord_only:
                    read_list.append(read.query_name)

                continue

        # Does read overlap bp
        # if phaser._read_overlap_bp2(d_read['qstart'], d_read['qend'], bp, orientation, map_min=map_min, clip_min=clip_min):
        if phaser._read_overlap_bp2(d_read['start'], d_read['end'], bp, orientation, map_min=map_min, clip_min=clip_min):
            if not discord_only:
                read_list.append(read.query_name)
            # elif is_consistent:
            #     read_list.append(read.query_name)

            continue

    set_alt = set(reads_alt)
    set_tot = set(read_list) | set_alt

    n_alt = len(set_alt)
    n_ref = len(set_tot) - n_alt
    AF = n_alt / len(set_tot)

    return n_ref, n_alt, AF, read_list

    # 1. Fetch reads in region around bp
    # 2. iterate over reads
    # 3. determine read bounds and mate bounds
    # 4. determine if mate is beyond bp --> count it
    # 5. if mate overlaps bp continue parsing until the mate is reached
    # 6. for reads overlapping bp, check if overlapping by 13 --> count it

def _pair_direction(read):
    """ Determine direction of pair (normal=FR, invers_f=FF, invers_r=RR, reversed=RF)
    """
    if read.is_reverse and read.mate_is_reverse:
        direction = 'inverse_r'
    elif not read.is_reverse and not read.mate_is_reverse:
        direction = 'inverse_f'
    elif read.is_reverse and (read.reference_start < read.next_reference_start):
        direction = 'reversed'
    elif read.mate_is_reverse and (read.next_reference_start < read.reference_start):
        direction = 'reversed'
    else:
        direction = 'normal'

    return direction

def _mate_beyond_bp(d_read, d_mate, pair_dir, bp, orientation, discord_only):
    """ Determine if mate lies beyond the breakpoint (or theoretically should be beyond if pair inverted)
    """
    beyond = False

    # Different chroms --> beyond
    if d_read['chrom'] != d_mate['chrom']:
        return False
        # return True

    # if pair_dir == 'inverse_f' and (d_read['start'] > d_mate['start']):
    #     start_m = d_read['end'] + (d_read['start'] - d_mate['qend'])
    #     end_m = start_m + (d_mate['qend'] - d_mate['qstart'])

    # elif pair_dir == 'inverse_r' and (d_read['start'] < d_mate['start']):
    #     end_m = d_read['qstart'] - (d_mate['qstart'] - d_read['qend'])
    #     start_m = end_m - (d_mate['qend'] - d_mate['qstart'])

    # elif pair_dir == 'reversed' and (d_read['start'] < d_mate['start']):
    #     end_m = d_read['qstart'] - (d_mate['qstart'] - d_read['qend'])
    #     start_m = end_m - (d_mate['qend'] - d_mate['qstart'])

    # elif pair_dir == 'reversed' and (d_read['start'] > d_mate['start']):
    #     start_m = d_read['end'] + (d_read['start'] - d_mate['qend'])
    #     end_m = start_m + (d_mate['qend'] - d_mate['qstart'])
    #     
    # else:
    #     start_m = d_mate['qstart']
    #     end_m = d_mate['qend']
    start_m = d_mate['qstart']
    end_m = d_mate['qend']

    buf = 0
    if discord_only:
        buf=75

    return phaser._read_beyond_bp(start_m, end_m, bp, orientation, buf=buf)

def _check_read_direction(read, orientation):
    is_consistent = False

    if orientation == 'downstream':
        if not read.is_reverse:
            is_consistent = True
    else:
        if read.is_reverse:
            is_consistent = True

    return is_consistent

def _read_ref_to_bp(start, end, bp, orientation):
    """ Is read on the reference side of a breakpoint
    """
    before = False

    if orientation == 'downstream':
        if end <= bp:
            before = True

    if orientation == 'upstream':
        if start >= bp:
            before = True            

    return before

def combine_meerkatBED_AF(df_bed, df_AF):
    """ Unify raw meerkat results (in BED format) with AF calculations of associated breakpoints 

        Inputs:
            df_bed: dataframe of meerkat results in BED format
            df_AF: dataframe of allele fraction results per fusion

        Output:
            dataframe in which each fusion has been matched with the raw event from which it arises
    """
    row_list = []
    idx_list = []
    for i, SV in df_bed.iterrows():
        bounds, _, _ = utils.bounds_by_SV_type(SV)
        fuse_start_str = "{}:{}-{}:{}".format(*bounds['fuse_start'])
        fuse_end_str = "{}:{}-{}:{}".format(*bounds['fuse_end'])

        for i, af in df_AF[df_AF.fusion == fuse_start_str].iterrows():
            if i in idx_list:
                continue

            row = SV[:-3].tolist() + af.tolist()
            row_list.append(row)
            idx_list.append(i)

        for i, af in df_AF[df_AF.fusion == fuse_end_str].iterrows():
            if i in idx_list:
                continue

            row = SV[:-3].tolist() + af.tolist()
            row_list.append(row)
            idx_list.append(i)

    columns = df_bed.columns[:-3].tolist() + df_AF.columns.tolist()
    df = pd.DataFrame(row_list, columns=columns)

    return df
