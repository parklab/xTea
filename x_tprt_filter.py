import global_values
##04/16/2019
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

####1. for tprt both (with polyA and TSD) events, we filter out those:
        ####1.1 polyA dominant cases, this module will mainly used on Alu
                # (as for L1 and SVA, some middle region or transduction will contain the polyA part)
####2. false positive ones comes from low reads quality part caused clip
####
####
class XTEARsltParser():
    def load_in_xTEA_rslt(self, sf_rslt):
        l_rcd=[]
        with open(sf_rslt) as fin_in:
            for line in fin_in:
                fields = line.split()

                n_lpolyA = int(fields[9])
                n_rpolyA = int(fields[10])
                n_ef_lclip=int(fields[5])
                n_ef_rclip=int(fields[6])
                #n_ef_clip = n_ef_lclip + n_ef_rclip
                n_ef_ldisc=int(fields[7])
                n_ef_rdisc=int(fields[8])
                #n_ef_disc = int(fields[7]) + int(fields[8])
                n_clip = int(fields[35])
                n_full_map = int(fields[36])
                n_disc = int(fields[39])
                n_concod = int(fields[40])
                l_rcd.append((n_lpolyA, n_rpolyA, n_ef_lclip, n_ef_rclip, n_ef_ldisc, n_ef_rdisc, n_clip, n_full_map,
                              n_disc, n_concod, line))
        return l_rcd
####
####
class XTPRTFilter():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1] != "/":
            self.swfolder += "/"
        self.n_jobs = n_jobs
########Hard code here!!!!!!!!!!!!!!!!!!!
        self.f_side_polyA_cutoff=global_values.ONE_SIDE_POLYA_CUTOFF
        ####


    def filter_by_polyA_AF(self, sf_in, sf_out):
        xtea_rslt=XTEARsltParser()
        l_rcd=xtea_rslt.load_in_xTEA_rslt(sf_in)
        af_filter=AFConflictFilter(self.swfolder, self.n_jobs)
        l_types = af_filter.get_rep_type()
        m_cutoff = af_filter.get_cutoff_by_type(l_types)
        with open(sf_out, "w") as fout_new:
            for rcd in l_rcd:
                s_info = rcd[-1]
                if self.is_polyA_dominant_two_side(rcd, self.f_side_polyA_cutoff)==True:
                    continue
                if af_filter.is_qualified_rcd(s_info, m_cutoff)==False:
                    continue
                fout_new.write(s_info)

####n_clip is clip reads cutoff
    def is_polyA_dominant_two_side(self, rcd, f_cutoff):
        n_lpolyA = rcd[0]
        n_rpolyA = rcd[1]
        n_ef_lclip = rcd[2]
        n_ef_rclip = rcd[3]

        #n_ef_clip = n_ef_lclip + n_ef_rclip
        if n_ef_lclip == 0 or n_ef_rclip==0:
            return False

        b_lpolyA= ((float(n_lpolyA) / float(n_ef_lclip)) > f_cutoff)
        b_rpolyA = ((float(n_rpolyA) / float(n_ef_rclip)) > f_cutoff)
        if b_lpolyA and b_rpolyA:
            return True

        return False
####

    ####this is count the overall
    def is_polyA_dominant(self, rcd, f_cutoff):
        n_lpolyA=rcd[0]
        n_rpolyA=rcd[1]
        n_ef_lclip=rcd[2]
        n_ef_rclip=rcd[3]
        n_ef_clip=n_ef_lclip+n_ef_rclip
        if n_ef_clip==0:
            return False
        n_polyA=n_lpolyA+n_rpolyA
        if float(n_polyA)/float(n_ef_clip) > f_cutoff:
            return True
        return False
####

####
class AFConflictFilter():
    def __init__(self, swfolder, n_jobs):
        self.swfolder = swfolder
        if self.swfolder[-1] != "/":
            self.swfolder += "/"
        self.n_jobs = n_jobs

    ####
    def get_rep_type(self):
        l_types = []
        l_types.append(global_values.ONE_SIDE_FLANKING)
        l_types.append(global_values.TWO_SIDE)
        l_types.append(global_values.TWO_SIDE_TPRT_BOTH)
        l_types.append(global_values.TWO_SIDE_TPRT)
        l_types.append(global_values.ONE_HALF_SIDE)
        l_types.append(global_values.ONE_HALF_SIDE_TRPT_BOTH)
        l_types.append(global_values.ONE_HALF_SIDE_TRPT)
        l_types.append(global_values.ONE_HALF_SIDE_POLYA_DOMINANT)
        l_types.append(global_values.ONE_SIDE)
        l_types.append(global_values.ONE_SIDE_COVERAGE_CONFLICT)
        l_types.append(global_values.ONE_SIDE_TRSDCT)
        l_types.append(global_values.ONE_SIDE_WEAK)
        l_types.append(global_values.ONE_SIDE_OTHER)
        l_types.append(global_values.ONE_SIDE_SV)
        l_types.append(global_values.ONE_SIDE_POLYA_DOMINANT)
        l_types.append(global_values.TWO_SIDE_POLYA_DOMINANT)
        l_types.append(global_values.HIGH_COV_ISD)
        l_types.append(global_values.OTHER_TYPE)
        return l_types

####
    ####
    def get_cutoff_by_type(self, l_types):
        m_cutoff = {}
        for s_type in l_types:
            if ("two_side" in s_type) or ("both-side" in s_type) or ("one_side_and_half_transduction" is s_type):
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
            elif "one_side" in s_type:  #
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
            elif ("one_half" in s_type) or ("one-half" in s_type):
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
            else:
                m_cutoff[s_type] = (0.075, 0.075, 0.075, 0.075)
        return m_cutoff

    def is_qualified_rcd(self, s_line, m_cutoff):
        fields = s_line.split()
        n_lpolyA = int(fields[9])
        n_rpolyA = int(fields[10])

        n_ef_clip = int(fields[5]) + int(fields[6])
        n_ef_disc = int(fields[7]) + int(fields[8])
        n_clip = int(fields[35])
        n_full_map = int(fields[36])
        n_disc = int(fields[39])
        n_concod = int(fields[40])

        f_ef_clip = 0.0
        if n_clip != 0:
            f_ef_clip = float(n_ef_clip) / float(n_clip)
        f_ef_disc = 0.0
        if n_disc != 0:
            f_ef_disc = float(n_ef_disc) / float(n_disc)
        f_clip_full_map = 0.0
        if (n_clip + n_full_map) != 0:
            f_clip_full_map = float(n_clip) / float(n_clip + n_full_map)
        f_disc_concod = 0.0
        if (n_disc + n_concod) != 0:
            f_disc_concod = float(n_disc) / float(n_disc + n_concod)

        s_type_ins = fields[32]
        b_pass = self.is_ins_pass_cutoff(m_cutoff, s_type_ins, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod)
        return b_pass
####
    ####
    def is_ins_pass_cutoff(self, m_cutoff, s_type, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod):
        (f_ef_clip_cutoff, f_ef_disc_cutoff, f_clip_full_cutoff, f_disc_concod_cutoff) = m_cutoff[s_type]
        b_ef_clip = (f_ef_clip > f_ef_clip_cutoff)
        b_ef_disc = f_ef_disc > f_ef_disc_cutoff
        b_clip_full = f_clip_full_map > f_clip_full_cutoff
        b_disc_concd = f_disc_concod > f_disc_concod_cutoff
        b_pass = b_ef_clip and b_ef_disc and b_clip_full and b_disc_concd
        return b_pass

    ####
    def filter_by_af_conflict(self, sf_hc, n_clip_cutoff, n_disc_cutoff):
        l_types = self.get_rep_type()
        m_cutoff = self.get_cutoff_by_type(l_types)
        self.calc_ratio(m_cutoff, n_clip_cutoff, n_disc_cutoff, sf_hc)
    ####
    def calc_ratio(self, m_cutoff, n_clip_cutoff, n_disc_cutoff, sf_in):
        sf_out = sf_in + ".after_filter"
        with open(sf_in) as fin_in, open(sf_out, "w") as fout_af_filter:
            n_total = 0
            n_hard_pass = 0
            n_pass = 0
            for line in fin_in:
                fields = line.split()
                n_lpolyA = int(fields[9])
                n_rpolyA = int(fields[10])

                n_ef_clip = int(fields[5]) + int(fields[6])
                n_ef_disc = int(fields[7]) + int(fields[8])
                n_clip = int(fields[35])
                n_full_map = int(fields[36])
                n_disc = int(fields[39])
                n_concod = int(fields[40])
                n_total += 1

                if n_ef_clip < n_clip_cutoff:
                    print "ef_clip", line
                    continue
                if n_ef_disc < n_disc_cutoff:
                    print "ef_disc", line
                    continue
                if n_lpolyA + n_rpolyA < 1:
                    print "no polyA", line
                    continue

                n_hard_pass += 1

                f_ef_clip = 0.0
                if n_clip != 0:
                    f_ef_clip = float(n_ef_clip) / float(n_clip)
                f_ef_disc = 0.0
                if n_disc != 0:
                    f_ef_disc = float(n_ef_disc) / float(n_disc)
                f_clip_full_map = 0.0
                if (n_clip + n_full_map) != 0:
                    f_clip_full_map = float(n_clip) / float(n_clip + n_full_map)
                f_disc_concod = 0.0
                if (n_disc + n_concod) != 0:
                    f_disc_concod = float(n_disc) / float(n_disc + n_concod)

                s_type_ins = fields[32]
                b_pass = self.is_ins_pass_cutoff(m_cutoff, s_type_ins, f_ef_clip, f_ef_disc, f_clip_full_map, f_disc_concod)

                if b_pass is True:
                    n_pass += 1
                    fout_af_filter.write(line)
                    #                 if "two" in s_type_ins:
                    #                     print line.rstrip()
                    #                     print fields[5], fields[6], fields[7],fields[8], fields[35],fields[36], fields[39], fields[40]
                    #                     print fields[0],fields[1],f_ef_clip,f_ef_disc,f_clip_full_map,f_disc_concod
                else:
                    s = 1
                    print line.rstrip()
                    #                 if "two" in s_type_ins:
                    #                     print fields[5], fields[6], fields[7],fields[8], fields[35],fields[36], fields[39], fields[40]
                    #                     print fields[0],fields[1],f_ef_clip,f_ef_disc,f_clip_full_map,f_disc_concod

            print n_hard_pass, n_pass, n_total
####