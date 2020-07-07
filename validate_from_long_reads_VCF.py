import sys
import os
from Bio import SeqIO
from subprocess import *
from multiprocessing import Pool

class LongReadValidation:
    def __init__(self, chrm, sf_folder):
        self.chrm=chrm
        self.sf_work_folder=sf_folder

    def read_pos(self, sf_f): ##read in the candidate positions
        m_temp={}
        with open(sf_f) as fin_f:
            for line in fin_f:
                fields=line.split()
                if len(fields)<2:
                    continue
                pos=int(fields[0])
                m_temp[pos]=1
        return m_temp

    def check_trio_consistent(self, sf_p1, sf_p2, sf_c):
        m_p1=self.read_pos(sf_p1)
        m_p2=self.read_pos(sf_p2)
        m_c=self.read_pos(sf_c)

        cnt=0
        for pos in m_c:
            for i in range(pos-10,pos+11):
                if m_p1.has_key(i) or m_p2.has_key(i):
                    cnt=cnt+1
                    break
        print cnt

    def gnrt_flank(self, sf_pos, sf_ref, chrm, flank_lenth, sf_out):
        chrm_seq=""
        for record in SeqIO.parse(sf_ref, "fasta"):
            if str(record.id)==chrm:
                chrm_seq=str(record.seq)
                break

        with open(sf_out,"w") as fout_flank:
            cnt=0
            with open(sf_pos) as fin_pos:
                for line in fin_pos:
                    fields=line.split()
                    pos=int(fields[0])

                    fout_flank.write(">"+str(cnt)+"_left\n")
                    fout_flank.write(chrm_seq[pos-flank_lenth-5:pos-5]+"\n")

                    fout_flank.write(">"+str(cnt)+"_right\n")
                    fout_flank.write(chrm_seq[pos+5:pos+flank_lenth+5]+"\n")
                    cnt=cnt+1


    def check_alignment(self, sf_algnmt, dist, lenth_cutoff):
        m_algn={}
        with open(sf_algnmt) as fin_algn:
            for line in fin_algn:
                fields=line.split()
                qname=fields[0]
                qlenth=int(fields[1])
                qstart=int(fields[2])
                qend=int(fields[3])
                strand=fields[4]
                read_id=fields[5] ####for mini ##map to which long read
                target_start=int(fields[7])
                target_end=int(fields[8])

                if float(abs(qend-qstart))/float(qlenth) < lenth_cutoff:
                    continue

                if float(abs(qend-qstart))/float(abs(target_end-target_start)) < lenth_cutoff:
                    continue

                qfields=qname.split("_")
                id=qfields[0] #site id, start from 0
                dir=qfields[1]#indicates left or right flank

                #need to set cutoff here to filter out aligned in a long range of reads
                #need to save whether reverse complementary of the alignment
                align_id=read_id+strand

                if m_algn.has_key(id)==False:
                    m_algn[id]={}
                if m_algn[id].has_key(dir)==False:
                    m_algn[id][dir]={}
                m_algn[id][dir][align_id]=(target_start, target_end)

        m_hit={}
        for id in m_algn: ##for each candidate site
            if len(m_algn[id])<2:#only have left or right or none reads hit
                continue
            for rid in m_algn[id]["left"]:
                if m_algn[id]["right"].has_key(rid):

                    lstart=m_algn[id]["left"][rid][0]
                    lend=m_algn[id]["left"][rid][1]
                    rstart=m_algn[id]["right"][rid][0]
                    rend=m_algn[id]["right"][rid][1]

                    gntp=0
                    pstart=0
                    pend=0
                    if lstart<rstart:
                        if lend<rstart:
                            pstart=lend
                            pend=rstart
                            if (pend-pstart)>=dist:
                                gntp=1
                        else:#overlap
                            if (lend-rstart)>dist:
                                continue

                    elif rstart<lstart:
                        if rend<lstart:
                            pstart=rend
                            pend=lstart
                            if (pend-pstart)>=dist:
                                gntp=1
                        else:#overlap
                            if (rend-lstart)>dist:
                                continue

                    read_id=rid[:-1] ##reduce the strand information
                    if m_hit.has_key(id)==False:
                        m_hit[id]={}
                    m_hit[id][read_id]=(pstart,pend,gntp) ## gntp 1 indicate the insertion exist
        return m_hit


    def gnrt_reverse_complementary(self, s):
        lth=len(s)
        s_rc=""
        for i in range(lth-1,-1,-1):
            if s[i]=="A" or s[i]=="a":
                s_rc=s_rc+"T"
            elif s[i]=="T" or s[i]=="t":
                s_rc=s_rc+"A"
            elif s[i]=="C" or s[i]=="c":
                s_rc=s_rc+"G"
            elif s[i]=="G" or s[i]=="g":
                s_rc=s_rc+"C"
            else:
                s_rc=s_rc+s[i]
        return s_rc


    def align_to_long_reads(self, sf_flank, sf_reads_index, sf_out):
        #sf_idx="/n/data1/hms/dbmi/park/simon_chu/projects/data/NA12878_pacbio/NA12878_pacbio_40x.mmi"
        cmd="minimap2 -x sr --cs {0} {1} > {2}".format(sf_reads_index, sf_flank, sf_out)
        Popen(cmd, shell = True, stdout = PIPE).communicate()


#define some global variables
sf_out_folder=""
sf_ref=""
allowed_dist=0
flank_lenth=300
align_cutoff=0.75
sf_reads_prefix="/n/data1/hms/dbmi/park/simon_chu/projects/data/NA12878_pacbio/NA12878_pacbio_40x"
def run_validate_for_chrm(chrm):
    global sf_out_folder
    global sf_ref
    global allowed_dist
    global flank_lenth
    global sf_reads_prefix
    lrv=LongReadValidation(chrm, sf_out_folder)

    #generate flank regions
    sf_pos=sf_out_folder+"{0}_sites.txt".format(chrm)
    if os.path.exists(sf_pos)==False:
        return (None,None)
    sf_out_flank=sf_out_folder+"{0}_flank_regions.fa".format(chrm)
    lrv.gnrt_flank(sf_pos, sf_ref, chrm, flank_lenth, sf_out_flank)

    #align to long reads
    #sf_reads_prefix="/data2/chongchu/pacbio/human.polished"
    #sf_reads_prefix="/data2/chongchu/1000G_pacbio_high_cov/NA19239_pacbio/NA19239_unique_id"
    sf_reads_index=sf_reads_prefix+".mmi"
    sf_out_algnmt=sf_out_folder+"{0}_align.mini".format(chrm)
##########uncomment for temporary usage
    #lrv.align_to_long_reads(sf_out_flank, sf_reads_index, sf_out_algnmt)

    #check alignment
    m_hit=lrv.check_alignment(sf_out_algnmt, allowed_dist, align_cutoff)
    #print chrm, len(m_hit)

    return (chrm, m_hit)

##generate the genotype
#also call out the repeat copies
def get_gntp_copy_seqs(m_rslts, sf_reads):
    n_chrm=len(m_rslts)
    m_read_id=[]
    m_gntp=[]
    m_copies=[]
    for i in range(n_chrm):
        m_read_id.append({})
        m_gntp.append({})
        m_copies.append({})

    for i in range(n_chrm):
        if m_rslts[i][1]==None:
            continue
        for id in m_rslts[i][1]:
            g0=0
            g1=0
            for rid in m_rslts[i][1][id]:
                if m_rslts[i][1][id][rid][2]==1: # indicates the copy exists
                    g1=g1+1
                    if m_read_id[i].has_key(rid)==False:
                        m_read_id[i][rid]={}
                    m_read_id[i][rid][id]=(m_rslts[i][1][id][rid][0], m_rslts[i][1][id][rid][1])
                else:
                    g0=g0+1
            m_gntp[i][id]=(g0, g1)

    for record in SeqIO.parse(sf_reads, "fasta"):
        rid=str(record.id)

        for i in range(n_chrm):
            if m_read_id[i].has_key(rid):
                for id in m_read_id[i][rid]:
                    start=m_read_id[i][rid][id][0]
                    end=m_read_id[i][rid][id][1]
                    seq=str(record.seq[start:end])

                    if m_copies[i].has_key(id)==False:
                        m_copies[i][id]={}
                    m_copies[i][id][rid]=seq

    return m_gntp, m_copies


def parallel_validate(sf_ref_fai, n_threads):
    global sf_reads_prefix
    l_chrm=[]
    with open(sf_ref_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            chrm=fields[0]
            l_chrm.append(chrm)

    bunch_size=1 #every time just run one command
    pool = Pool(n_threads)
    rslts=pool.map(run_validate_for_chrm, l_chrm, bunch_size) #[{}]
    pool.close()
    pool.join()

    sf_reads=sf_reads_prefix+".fa"
    m_all_gntp, m_all_copies=get_gntp_copy_seqs(rslts, sf_reads)

    n_chrm=len(m_all_copies)
    for i in range(n_chrm):
        #generate copies
        chrm=rslts[i][0]
        if chrm==None:
            continue
        m_gntp=m_all_gntp[i]
        m_copies=m_all_copies[i]

        sf_copies_folder=sf_out_folder+"{0}_copies".format(chrm)
        if os.path.exists(sf_copies_folder)==False:
            cmd="mkdir {0}".format(sf_copies_folder)
            Popen(cmd, shell = True, stdout = PIPE).communicate()

        for id in m_copies:
            g0=m_gntp[id][0]
            g1=m_gntp[id][1]
            with open(sf_copies_folder+"/"+str(id)+"_"+str(g0)+"_"+str(g1)+".fa","w") as fout_copy:
                for rid in m_copies[id]:
                    fout_copy.write(">"+rid+"\n")
                    fout_copy.write(m_copies[id][rid]+"\n")

        sf_site_gntf=sf_out_folder+"{0}_sites_gntp.txt".format(chrm)
        with open(sf_site_gntf,"w") as fout_gntp:
            for id in m_gntp:
                fout_gntp.write(str(id)+" "+str(m_gntp[id][0])+" "+str(m_gntp[id][1])+"\n")


def call_final_validated_sites(sf_fai, sf_out_folder):
    sf_final="{0}final_validated.vcf".format(sf_out_folder)
    sf_final_copy="{0}/final_validate_copies/".format(sf_out_folder)
    if os.path.exists(sf_final_copy)==False:
        cmd="mkdir {0}".format(sf_final_copy)
        Popen(cmd, shell = True, stdout = PIPE).communicate()

    with open(sf_fai) as fin_fai, open(sf_final,"w") as fout_final:
        for line in fin_fai:
            fields=line.split()
            chrm=fields[0]

            sf_pos="{0}{1}_sites.txt".format(sf_out_folder,chrm)
            sf_gntp="{0}{1}_sites_gntp.txt".format(sf_out_folder,chrm)
            l_pos=[]
            if os.path.isfile(sf_pos)==False:
                continue
            with open(sf_pos) as fin_pos:
                for pos_line in fin_pos:
                    l_pos.append(pos_line.rstrip())
            with open(sf_gntp) as fin_gntp:
                for gntp_line in fin_gntp:
                    gntp_fields=gntp_line.split()
                    gindex=int(gntp_fields[0])
                    g0=int(gntp_fields[1])
                    g1=int(gntp_fields[2])

                    pos_fileds=l_pos[gindex].split()
                    if g1>1 and g1>g0:#here for temporary cutoff################################################################################
                        sf_copy_seq="{0}{1}_copies/{2}_{3}_{4}.fa".format(sf_out_folder, chrm, gindex, g0, g1)
                        sf_new_seq="{0}final_validate_copies/{1}_{2}_{3}_{4}.fa".format(sf_out_folder, chrm, pos_fileds[0], g0, g1)
                        cmd="cp {0} {1}".format(sf_copy_seq, sf_new_seq)
                        Popen(cmd, shell = True, stdout = PIPE).communicate()

                        fout_final.write(chrm+" " + l_pos[gindex]+" "+str(g0)+" "+str(g1)+"\n")


def load_vcf(sf_vcf, itype, m_chrm_pos):
    with open(sf_vcf) as fin_alu:
        for line in fin_alu:
            if line[0]=="#":
                continue
            fields=line.split()
            chrm=fields[0]
            pos=int(fields[1])
            if m_chrm_pos.has_key(chrm)==False:
                m_chrm_pos[chrm]={}
            m_chrm_pos[chrm][pos]=itype

def gnrt_pos_from_vcf(sf_vcf):
    # sf_alu=sf_out_folder+"ALU.final_comp.vcf"
    # sf_L1=sf_out_folder+"LINE1.final_comp.vcf"
    # sf_sva=sf_out_folder+"SVA.final_comp.vcf"
    #sf_pos=sf_out_folder+"ref_filtered_alu.ref_filtered_L1.filtered_SVA.vcf"
    m_chrm_pos={}
    load_vcf(sf_vcf, 0, m_chrm_pos)
    # load_vcf(sf_alu, 0, m_chrm_pos)
    # load_vcf(sf_L1, 1, m_chrm_pos)
    # load_vcf(sf_sva, 2, m_chrm_pos)

    for chrm in m_chrm_pos:
        sf_sites=sf_out_folder+"{0}_sites.txt".format(chrm)
        with open(sf_sites,"w") as fout_sites:
            for pos in m_chrm_pos[chrm]:
                fout_sites.write(str(pos)+" "+"-1 -1 -1 -1 -1 -1 \n")


if __name__ == "__main__":
    flank_lenth=400
    sf_ref_fai=sys.argv[1]
    sf_ref=sys.argv[2]
    sf_sites=sys.argv[3]
    allowed_dist=int(sys.argv[4]) #Max allowed distance between the two aligned flanks
    n_threads=int(sys.argv[5])
    sf_out_folder=sys.argv[6]

    #align_cutoff=0.75
    #sf_reads_prefix="/data2/chongchu/1000G_pacbio_high_cov/NA19239_pacbio/NA19239_unique_id"

    if sf_out_folder[-1]!="/":
        sf_out_folder=sf_out_folder+"/"

    gnrt_pos_from_vcf(sf_sites)
    parallel_validate(sf_ref_fai, n_threads)
    call_final_validated_sites(sf_ref_fai, sf_out_folder)

    # sf_p1=sys.argv[1]
    # sf_p2=sys.argv[2]
    # sf_c=sys.argv[3]
    # check_trio_consistent(sf_p1, sf_p2, sf_c)
