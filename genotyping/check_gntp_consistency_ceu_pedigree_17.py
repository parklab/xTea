import os

INCONSIST_TYPE1=0
INCONSIST_TYPE2=1
INCONSIST_TYPE3=2
INCONSIST_TYPE4=3
INCONSIST_TYPE5=4

EXTND = 200
BRKPNT_SLACK=100
LOG_ON=True

def process_chrm_name(chrm, b_with_chr):
    b_chrm_with_chr = False
    if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
        b_chrm_with_chr = True

    # print chrm, self.b_with_chr, b_chrm_with_chr ############################################################################
    if b_with_chr == True and b_chrm_with_chr == True:
        return chrm
    elif b_with_chr == True and b_chrm_with_chr == False:
        return "chr" + chrm
    elif b_with_chr == False and b_chrm_with_chr == True:
        return chrm[3:]
    else:
        return chrm

####
#m_hit save the position of in m_dict1
def cmp_dict(m_dict1, m_bcmk):
    n_overlap = 0
    m_unhit = {}
    m_hit={}

    # print m_dict1.keys(), "dict keys"
    # print m_bcmk.keys(), "bcmk keys"

    b_contain_chr = False
    if "chr1" in m_bcmk:
        b_contain_chr = True

    for chrm_1 in m_bcmk:
        chrm = process_chrm_name(chrm_1, b_contain_chr)
        if chrm not in m_dict1:
            #print "chromosome {0} is not included in database\n".format(chrm)
            continue

        for pos in m_bcmk[chrm]:
            b_hit = False
            for i in range(-1 * EXTND, EXTND):
                if (pos + i) in m_dict1[chrm]:
                    n_overlap += 1
                    b_hit = True
                    if chrm not in m_hit:
                        m_hit[chrm]={}
                    m_hit[chrm][pos+i]=m_dict1[chrm][pos+i]
                    break

            if b_hit == False:
                if chrm not in m_unhit:
                    m_unhit[chrm] = {}
                m_unhit[chrm][pos] = m_bcmk[chrm][pos]
    return n_overlap, m_unhit, m_hit

def load_in_bcmk(l_bcmk):
    m_bcmk={}
    for rcd in l_bcmk:
        s_sample=rcd[0]
        sf_sample=rcd[1]
        m_bcmk[s_sample]={}
        with open(sf_sample) as fin_sample:
            for line in fin_sample:
                fields=line.split()
                chrm=fields[0]
                pos=int(fields[1])
                if chrm not in m_bcmk[s_sample]:
                    m_bcmk[s_sample][chrm]={}
                m_bcmk[s_sample][chrm][pos]=1
    return m_bcmk

def load_rslts_from_vcf(sf_vcf):
    m_rslts = {}
    with open(sf_vcf) as fin_vcf:
        for line in fin_vcf:
            if line[0]=="#":
                continue
            if "orphan" in line:
                continue
        
            fields=line.split()
            ins_chrm=fields[0]
            ins_pos=int(fields[1])
            gntp_fields=fields[-1].split(":")
            gntp=gntp_fields[0]
            if gntp=="0/0":
                continue
            if ins_chrm not in m_rslts:
                m_rslts[ins_chrm]={}
            m_rslts[ins_chrm][ins_pos]=gntp
    return m_rslts


########
def load_in_relationship(sf_relationship):
    m_relationship={}
    m_children={}
    m_parents={}#
    m_samples={}
    with open(sf_relationship) as fin_relationship:
        for line in fin_relationship:
            fields=line.split()
            s_child=fields[0]
            s_father=fields[1]
            s_mother=fields[2]
            m_samples[s_child]=1
            m_samples[s_father]=1
            m_samples[s_mother]=1

            m_relationship[s_child]=(s_father, s_mother)
            m_parents[s_father]=s_child
            m_parents[s_mother]=s_child
    return m_relationship, m_samples

########
def prepare_xTEA_rslt_list(sf_rslt_list):
    m_xTEA_rslt_list={}
    with open(sf_rslt_list) as fin_xTEA:
        for line in fin_xTEA:
            fields=line.split("/")
            s_id=fields[-3]
            if s_id not in m_xTEA_rslt_list:
                m_xTEA_rslt_list[s_id]={}
            s_type=fields[-2] #L1, Alu, SVA
            m_xTEA_rslt_list[s_id][s_type]=line.rstrip()
    return m_xTEA_rslt_list

####
def prepare_Melt_rslt_list(sf_melt_prefix, m_samples, l_rep_type2, l_rep_type):
    m_melt_rslt_list={}
    for s_sample in m_samples:
        for s_type, s_type_m in zip(l_rep_type2, l_rep_type):
            sf_rslt=sf_melt_prefix+s_sample+"/{0}.final_comp.vcf".format(s_type_m)#
            if s_sample not in m_melt_rslt_list:
                m_melt_rslt_list[s_sample]={}
            m_melt_rslt_list[s_sample][s_type]=sf_rslt
    return m_melt_rslt_list

def find_hit_parents(m_rslt, ins_chrm, ins_pos):
    if ins_chrm not in m_rslt:
        return False
    i_slack=100
    
    for tmp_pos in range(ins_pos-i_slack, ins_pos+i_slack):
        if tmp_pos in m_rslt[ins_chrm]:
            return True
    return False

def is_gntp_consistent(s_child, s_father, s_mother):
    if s_father=="1/1" and s_mother=="0/0" and s_child!="0/1":
        return False, INCONSIST_TYPE1
    if s_father=="0/0" and s_mother=="1/1" and s_child!="0/1":
        return False, INCONSIST_TYPE1
    
    if s_father=="1/1" and s_mother=="0/1" and s_child=="0/0":
        return False, INCONSIST_TYPE2
    if s_father=="0/1" and s_mother=="1/1" and s_child=="0/0":
        return False, INCONSIST_TYPE2
    
    if s_father=="1/1" and s_mother=="1/1" and s_child!="1/1":
        return False, INCONSIST_TYPE3
    
    if s_father=="0/0" and s_mother=="0/0" and s_child!="0/0":
        return False, INCONSIST_TYPE4
    
    if s_father=="0/0" and s_mother=="0/1" and s_child=="1/1":
        return False, INCONSIST_TYPE5
    if s_father=="0/1" and s_mother=="0/0" and s_child=="1/1":
        return False, INCONSIST_TYPE5
    return True, -1

####
def check_overlap(m_xTEA, m_melt):
    #n_bcmk, m_unhit1, m_hit1 = self.cmp_dict(m_xTEA, m_xTEA)
    n_hit, m_melt_specific, m_overlap = cmp_dict(m_xTEA, m_melt)
    n_hit_fp, m_xTEA_specific, m_hit3 = cmp_dict(m_melt, m_xTEA)
    return m_overlap, m_xTEA_specific, m_melt_specific

def check_gntp_for_site(chrm, pos, m_sites):
    if chrm not in m_sites:
        return '0/0'
    for tmp_pos in m_sites[chrm]:
        if abs(tmp_pos-pos)<BRKPNT_SLACK:
            return m_sites[chrm][tmp_pos]
    return '0/0'

def get_overlapped_sites(m_sites, m_bcmk):
    b_contain_chr=False
    for chrm in m_sites:
        if len(chrm)>3 and chrm[:3]=="chr":
            b_contain_chr=True

    m_overlap={}
    for chrm_1 in m_bcmk:
        chrm = process_chrm_name(chrm_1, b_contain_chr)
        if chrm not in m_sites:
            continue

        for tmp_pos in m_bcmk[chrm_1]:
            for pos in range(tmp_pos-EXTND, tmp_pos+EXTND):
                if pos in m_sites[chrm]:
                    if chrm not in m_overlap:
                        m_overlap[chrm]={}
                    m_overlap[chrm][pos]=1
                    break
    return m_overlap

####
def check_gntp_consistent(m_sites, s_rep_type, s_father, s_mother, s_p_grandpa, s_p_grandma, s_m_grandpa, s_m_grandma, m_rslt_list, l_bcmk=None, s_child="not_set"):
    global LOG_ON
    n_consist_level1=0#number of consistent sites
    m_level1_inconsist={}
    n_consist_level2=0#number of consistent sites
    m_level2_inconsist={}
    sf_father=m_rslt_list[s_father][s_rep_type]
    sf_mother=m_rslt_list[s_mother][s_rep_type]
    sf_p_grandpa=m_rslt_list[s_p_grandpa][s_rep_type]
    sf_p_grandma=m_rslt_list[s_p_grandma][s_rep_type]
    sf_m_grandpa=m_rslt_list[s_m_grandpa][s_rep_type]
    sf_m_grandma=m_rslt_list[s_m_grandma][s_rep_type]

    m_bcmk = None#this is the results of third-party benchmark data
    if l_bcmk is not None:
        m_bcmk = load_in_bcmk(l_bcmk)

    m_father1=load_rslts_from_vcf(sf_father)
    m_father=m_father1
    if s_father in m_bcmk:
        m_father=get_overlapped_sites(m_father1, m_bcmk[s_father])
    m_mother1 = load_rslts_from_vcf(sf_mother)
    m_mother=m_mother1
    if s_mother in m_bcmk:
        m_mother=get_overlapped_sites(m_mother1, m_bcmk[s_mother])

    m_p_grandpa=load_rslts_from_vcf(sf_p_grandpa)
    m_p_grandma=load_rslts_from_vcf(sf_p_grandma)
    m_m_grandpa=load_rslts_from_vcf(sf_m_grandpa)
    m_m_grandma=load_rslts_from_vcf(sf_m_grandma)
    ####
    for chrm in m_sites:
        for pos in m_sites[chrm]:
            s_gntp_child=m_sites[chrm][pos]
            #level 1 checking: check father and mother
            s_gntp_fa=check_gntp_for_site(chrm, pos, m_father)
            s_gntp_ma=check_gntp_for_site(chrm, pos, m_mother)
            b_lev1, in_consist_type = is_gntp_consistent(s_gntp_child, s_gntp_fa, s_gntp_ma)
            if b_lev1==True:
                n_consist_level1 += 1
                if LOG_ON==True:
                    print "level1 passed"
            else:
                if in_consist_type not in m_level1_inconsist:
                    m_level1_inconsist[in_consist_type]=0
                m_level1_inconsist[in_consist_type]+=1
                continue
            #level 2 checking, check grandparents
            s_gntp_p_gpa=check_gntp_for_site(chrm, pos, m_p_grandpa)
            s_gntp_p_gma=check_gntp_for_site(chrm, pos, m_p_grandma)
            s_gntp_m_gpa=check_gntp_for_site(chrm, pos, m_m_grandpa)
            s_gntp_m_gma=check_gntp_for_site(chrm, pos, m_m_grandma)
            b_p_lev2, incsist_type_p=is_gntp_consistent(s_gntp_fa, s_gntp_p_gpa, s_gntp_p_gma)
            b_m_lev2, incsist_type_m=is_gntp_consistent(s_gntp_ma, s_gntp_m_gpa, s_gntp_m_gma)
            
            if LOG_ON==True:
                print s_gntp_child, s_gntp_fa, s_gntp_ma, s_child, chrm, pos, s_rep_type, "test1"
                print s_gntp_fa, s_gntp_p_gpa, s_gntp_p_gma, s_father, s_p_grandpa, s_p_grandma, "test2"
                print s_gntp_ma, s_gntp_m_gpa, s_gntp_m_gma, s_mother, s_m_grandpa, s_m_grandma, "test3"
            
            
            if (b_p_lev2==True) and (b_m_lev2==True):
                if LOG_ON==True:
                    print "level 2 passed"
                n_consist_level2+=1
            if b_p_lev2==False:
                if LOG_ON==True:
                    print "level 2 parternal failed"
                if incsist_type_p not in m_level2_inconsist:
                    m_level2_inconsist[incsist_type_p]=0
                m_level2_inconsist[incsist_type_p]+=1
            elif (b_m_lev2==False):
                if LOG_ON==True:
                    print "level 2 maternal failed"
                if incsist_type_m not in m_level2_inconsist:
                    m_level2_inconsist[incsist_type_m]=0
                m_level2_inconsist[incsist_type_m]+=1
    LOG_ON=False
    return n_consist_level1, n_consist_level2, m_level1_inconsist, m_level2_inconsist
########

####
def count_inconsistent_type(m_inconsistent):#
    l_cnt=[0,0,0,0,0]
    for i_type in m_inconsistent:
        l_cnt[i_type]=m_inconsistent[i_type]
    return l_cnt

####
def cnt_sites(m_sites):
    n_site=0
    for chrm in m_sites:
        n_site+=len(m_sites[chrm])
    return n_site

####
def count_gntp_consistent_by_category(m_relationship, m_xTEA_rslt_list, m_melt_rslt_list, l_bcmk=None):
    #print m_melt_rslt_list
    m_sample_cnt={}
    for s_child in m_xTEA_rslt_list:
        if s_child not in m_relationship:
            continue
        s_father=m_relationship[s_child][0]
        s_mother=m_relationship[s_child][1]
        if (s_father not in m_relationship) or (s_mother not in m_relationship):
            continue

        s_p_grandpa=m_relationship[s_father][0]
        s_p_grandma=m_relationship[s_father][1]
        s_m_grandpa=m_relationship[s_mother][0]
        s_m_grandma=m_relationship[s_mother][1]
        for s_type in m_xTEA_rslt_list[s_child]:
            if s_child not in m_melt_rslt_list:
                print s_child, "not found in Melt list"
            if s_type not in m_melt_rslt_list[s_child]:
                print s_type, "not found in Melt list of sample", s_child
            #compare the two results and find overlap, xTEA_specific, Melt_specific
            sf_xTEA_child=m_xTEA_rslt_list[s_child][s_type]
            sf_melt_child=m_melt_rslt_list[s_child][s_type]
            m_xTEA_child=load_rslts_from_vcf(sf_xTEA_child)
            m_melt_child=load_rslts_from_vcf(sf_melt_child)
            m_overlap, m_xTEA_specific, m_melt_specific=check_overlap(m_xTEA_child, m_melt_child)

            #check consistent for overlap, xTEA_specific, and melt_specific
            #for overlap ones, we check both xTEA and Melt results
            #first check overlap for xTEA 
            n_xTEA_ovlp_lv1, n_xTEA_ovlp_lv2, m_ovlp_lv1_xTEA, m_ovlp_lv2_xTEA = \
                check_gntp_consistent(m_overlap, s_type, s_father, s_mother, s_p_grandpa, s_p_grandma, s_m_grandpa, s_m_grandma, m_xTEA_rslt_list, l_bcmk, s_child)
            l_err_ovlp_lv1_xTEA=count_inconsistent_type(m_ovlp_lv1_xTEA)
            l_err_ovlp_lv2_xTEA=count_inconsistent_type(m_ovlp_lv2_xTEA)
            n_melt_ovlp_lv1, n_melt_ovlp_lv2, m_ovlp_lv1_melt, m_ovlp_lv2_melt = \
                check_gntp_consistent(m_overlap, s_type, s_father, s_mother, s_p_grandpa, s_p_grandma, s_m_grandpa, s_m_grandma, m_melt_rslt_list, l_bcmk)
            
            n_xTEA_spec_lv1, n_xTEA_spec_lv2, m_spec_lv1_xTEA, m_spec_lv2_xTEA = \
                check_gntp_consistent(m_xTEA_specific, s_type, s_father, s_mother, s_p_grandpa, s_p_grandma, s_m_grandpa, s_m_grandma, m_xTEA_rslt_list, l_bcmk)
            
            n_melt_spec_lv1, n_melt_spec_lv2, m_spec_lv1_melt, m_spec_lv2_melt = \
                check_gntp_consistent(m_melt_specific, s_type, s_father, s_mother, s_p_grandpa, s_p_grandma, s_m_grandpa, s_m_grandma, m_melt_rslt_list, l_bcmk)
            
            #by sample:rep-type:(xxx,xxx,xxx)
            if s_child not in m_sample_cnt:
                m_sample_cnt[s_child]={}
            n_overlap=cnt_sites(m_overlap)
            n_xTEA_spec=cnt_sites(m_xTEA_specific)
            n_melt_spec=cnt_sites(m_melt_specific)
            l_ovlp_rslt=(n_overlap, n_xTEA_spec, n_melt_spec)
            l_consistnt_rslt=(n_xTEA_ovlp_lv1, n_xTEA_ovlp_lv2, n_melt_ovlp_lv1, n_melt_ovlp_lv2, n_xTEA_spec_lv1, n_xTEA_spec_lv2, n_melt_spec_lv1, n_melt_spec_lv2)
            m_sample_cnt[s_child][s_type]=(l_ovlp_rslt, l_consistnt_rslt)
    return m_sample_cnt

####
def dump_statistic_rslts(m_cnt, sf_rslts):
    with open(sf_rslts,"w") as fout_rslts:#
        fout_rslts.write("Sample,Repeat_type,Num_overlap,Num_xTEA_spec,Num_Melt_spec,Ovlp_xTEA_level1,Ovlp_xTEA_level2,Ovlp_Melt_level1,Ovlp_Melt_level2,Spec_xTEA_level1,Spec_xTEA_level2,Spec_Melt_level1,Spec_Melt_level2\n")
        for s_sample in m_cnt:
            for s_type in m_cnt[s_sample]:
                (n_overlap, n_xTEA_spec, n_melt_spec)=m_cnt[s_sample][s_type][0]
                (n_xTEA_ovlp_lv1, n_xTEA_ovlp_lv2, n_melt_ovlp_lv1, n_melt_ovlp_lv2, n_xTEA_spec_lv1, n_xTEA_spec_lv2, n_melt_spec_lv1, n_melt_spec_lv2)=m_cnt[s_sample][s_type][1]
                
                f_ovlp_lv1_xTEA=float(n_xTEA_ovlp_lv1)/float(n_overlap)
                f_ovlp_lv2_xTEA=float(n_xTEA_ovlp_lv2)/float(n_overlap)
                f_ovlp_lv1_melt=float(n_melt_ovlp_lv1)/float(n_overlap)
                f_ovlp_lv2_melt=float(n_melt_ovlp_lv2)/float(n_overlap)
                f_spec_lv1_xTEA=float(n_xTEA_spec_lv1)/float(n_xTEA_spec)
                f_spec_lv2_xTEA=float(n_xTEA_spec_lv2)/float(n_xTEA_spec)
                f_spec_lv1_melt=float(n_melt_spec_lv1)/float(n_melt_spec)
                f_spec_lv2_melt=float(n_melt_spec_lv2)/float(n_melt_spec)
                
                sinfo="%f,%f,%f,%f,%f,%f,%f,%f" % (f_ovlp_lv1_xTEA, f_ovlp_lv2_xTEA, f_ovlp_lv1_melt, f_ovlp_lv2_melt, f_spec_lv1_xTEA, f_spec_lv2_xTEA, f_spec_lv1_melt, f_spec_lv2_melt)
                fout_rslts.write(s_sample+","+s_type+","+str(n_overlap)+","+str(n_xTEA_spec)+","+str(n_melt_spec)+","+sinfo+"\n")
####
def dump_statistic_rslts2(m_cnt, sf_rslts):
    with open(sf_rslts,"w") as fout_rslts:#
        fout_rslts.write("Sample,Repeat_type,Level,Num_overlap,Num_xTEA_spec,Num_Melt_spec,Ovlp_xTEA,Ovlp_Melt,Spec_xTEA,Spec_Melt\n")
        for s_sample in m_cnt:
            for s_type in m_cnt[s_sample]:
                for s_level in ["Level1","Level2"]:
                    (n_overlap, n_xTEA_spec, n_melt_spec)=m_cnt[s_sample][s_type][0]
                    (n_xTEA_ovlp_lv1, n_xTEA_ovlp_lv2, n_melt_ovlp_lv1, n_melt_ovlp_lv2, n_xTEA_spec_lv1, n_xTEA_spec_lv2, n_melt_spec_lv1, n_melt_spec_lv2)=m_cnt[s_sample][s_type][1]

                    f_ovlp_lv1_xTEA=float(n_xTEA_ovlp_lv1)/float(n_overlap)
                    f_ovlp_lv2_xTEA=float(n_xTEA_ovlp_lv2)/float(n_overlap)
                    f_ovlp_lv1_melt=float(n_melt_ovlp_lv1)/float(n_overlap)
                    f_ovlp_lv2_melt=float(n_melt_ovlp_lv2)/float(n_overlap)
                    f_spec_lv1_xTEA=float(n_xTEA_spec_lv1)/float(n_xTEA_spec)
                    f_spec_lv2_xTEA=float(n_xTEA_spec_lv2)/float(n_xTEA_spec)
                    f_spec_lv1_melt=float(n_melt_spec_lv1)/float(n_melt_spec)
                    f_spec_lv2_melt=float(n_melt_spec_lv2)/float(n_melt_spec)

                    #sinfo="%f,%f,%f,%f,%f,%f,%f,%f" % (f_ovlp_lv1_xTEA, f_ovlp_lv2_xTEA, f_ovlp_lv1_melt, f_ovlp_lv2_melt, f_spec_lv1_xTEA, f_spec_lv2_xTEA, f_spec_lv1_melt, f_spec_lv2_melt)
                    sinfo="%f,%f,%f,%f" % (f_ovlp_lv1_xTEA, f_ovlp_lv1_melt, f_spec_lv1_xTEA, f_spec_lv1_melt)
                    if s_level=="Level2":
                        sinfo="%f,%f,%f,%f" % (f_ovlp_lv2_xTEA, f_ovlp_lv2_melt, f_spec_lv2_xTEA, f_spec_lv2_melt)
                    fout_rslts.write(s_sample+","+s_type+","+s_level+","+str(n_overlap)+","+str(n_xTEA_spec)+","+str(n_melt_spec)+","+sinfo+"\n")

####
def dump_statistic_rslts3(m_cnt, sf_rslts):
    with open(sf_rslts,"w") as fout_rslts:#
        fout_rslts.write("Sample,Repeat_type,Level,Overlap_type,Method,Num,Ratio\n")
        for s_sample in m_cnt:
            for s_type in m_cnt[s_sample]:
                for s_level in ["Level1","Level2"]:
                    #for s_ovlp_type in ["xTEA_Overlap", "Melt_Overlap", "xTEA_Specific", "Melt_Specific"]
                    (n_overlap, n_xTEA_spec, n_melt_spec)=m_cnt[s_sample][s_type][0]
                    (n_xTEA_ovlp_lv1, n_xTEA_ovlp_lv2, n_melt_ovlp_lv1, n_melt_ovlp_lv2, n_xTEA_spec_lv1, n_xTEA_spec_lv2, n_melt_spec_lv1, n_melt_spec_lv2)=m_cnt[s_sample][s_type][1]

                    f_ovlp_lv1_xTEA=float(n_xTEA_ovlp_lv1)/float(n_overlap)
                    f_ovlp_lv2_xTEA=float(n_xTEA_ovlp_lv2)/float(n_overlap)
                    f_ovlp_lv1_melt=float(n_melt_ovlp_lv1)/float(n_overlap)
                    f_ovlp_lv2_melt=float(n_melt_ovlp_lv2)/float(n_overlap)
                    f_spec_lv1_xTEA=float(n_xTEA_spec_lv1)/float(n_xTEA_spec)
                    f_spec_lv2_xTEA=float(n_xTEA_spec_lv2)/float(n_xTEA_spec)
                    f_spec_lv1_melt=float(n_melt_spec_lv1)/float(n_melt_spec)
                    f_spec_lv2_melt=float(n_melt_spec_lv2)/float(n_melt_spec)
                    
                    #save xTEA overlap
                    s_ovlp_type="xTEA_overlap"
                    sinfo="%d,%f" % (n_overlap, f_ovlp_lv1_xTEA)
                    if s_level=="Level2":
                        sinfo="%d,%f" % (n_overlap,f_ovlp_lv2_xTEA)
                    fout_rslts.write(s_sample+","+s_type+","+s_level+","+s_ovlp_type+",xTEA,"+sinfo+"\n")
                    #save Melt overlap
                    s_ovlp_type="Melt_overlap"
                    sinfo="%d,%f" % (n_overlap, f_ovlp_lv1_melt)
                    if s_level=="Level2":
                        sinfo="%d,%f" % (n_overlap, f_ovlp_lv2_melt)
                    fout_rslts.write(s_sample+","+s_type+","+s_level+","+s_ovlp_type+",Melt,"+sinfo+"\n")
                    #save xTEA specific
                    s_ovlp_type="xTEA_specific"
                    sinfo="%d,%f" % (n_xTEA_spec, f_spec_lv1_xTEA)
                    if s_level=="Level2":
                        sinfo="%d,%f" % (n_xTEA_spec, f_spec_lv2_xTEA)
                    fout_rslts.write(s_sample+","+s_type+","+s_level+","+s_ovlp_type+",xTEA,"+sinfo+"\n")
                    #save Melt specific
                    s_ovlp_type="Melt_specific"
                    sinfo="%d,%f" % (n_melt_spec, f_spec_lv1_melt)
                    if s_level=="Level2":
                        sinfo="%d,%f" % (n_melt_spec, f_spec_lv2_melt)
                    fout_rslts.write(s_sample+","+s_type+","+s_level+","+s_ovlp_type+",Melt,"+sinfo+"\n")
####
                    
##main function
if __name__ == '__main__':
    ####str(n_overlap)+","+str(n_xTEA_spec)+","+str(n_melt_spec)
    ####
    sf_relationship="/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/CEU_Pedigree_17/sample_relationship.txt"
    m_relationship, m_samples=load_in_relationship(sf_relationship)
    ####
    sf_prefix="/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/CEU_Pedigree_17/"
    #sf_rep_type="Alu"
    #sf_melt_rep_type="ALU"
    l_rep_type2=["Alu", "L1", "SVA"]
    l_rep_type=["ALU", "LINE1", "SVA"]
    #l_rep_type=["Alu", "L1", "SVA"]

    #prepare xTEA result list
    sf_xTEA_rslt_list="xTEA_vcf.list"
    m_xTEA_rslt_list=prepare_xTEA_rslt_list(sf_xTEA_rslt_list)
    ####
    #prepare Melt results list
    sf_melt_prefix="/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/CEU_Pedigree_17/melt_results/"
    m_melt_rslt_list=prepare_Melt_rslt_list(sf_melt_prefix, m_samples, l_rep_type2, l_rep_type)

    ####
    sf_NA12878_pbsv_ccs="/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/GIAB/released_rslts/NA12878/HG001_hs37d5.pbsv_INS.vcf"#NA12878 results
    l_bcmk=[]
    l_bcmk.append(("NA12878", sf_NA12878_pbsv_ccs))
    m_cnted_rslts=count_gntp_consistent_by_category(m_relationship, m_xTEA_rslt_list, m_melt_rslt_list, l_bcmk)
    sf_out="Melt_xTEA_genotype_precision_ratio2.csv"
    dump_statistic_rslts2(m_cnted_rslts, sf_out)
    sf_out3 = "Melt_xTEA_genotype_precision_ratio3.csv"
    dump_statistic_rslts3(m_cnted_rslts, sf_out3)
    ####
    ####
####