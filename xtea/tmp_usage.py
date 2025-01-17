import pysam

def get_centromere_gap_flanks():
    m_pos={}
    l_centromere=[]
    with open("centromeres.txt") as fin_centromere:
        for line in fin_centromere:
            fields=line.split()
            chrm=fields[1]
            if chrm not in m_pos:
                m_pos[chrm]=[]
            i_start=int(fields[2])
            i_end=int(fields[3])
            m_pos[chrm].append(i_start)
            m_pos[chrm].append(i_end)

        for chrm in m_pos:
            l_tmp=m_pos[chrm]
            l_tmp.sort()
            i_start=l_tmp[0]
            i_end=l_tmp[-1]
            l_centromere.append((chrm, i_start, i_end))

    print(l_centromere)
    ####
    sf_ref="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    sf_out="ref_centromere_anchor.fa"
    with open(sf_out,"w") as fout:
        f_fa = pysam.FastaFile(sf_ref)
        for rcd in l_centromere:
            chrm=rcd[0]
            i_l_end=rcd[1]
            i_l_start=i_l_end-20000
            i_r_start=rcd[2]
            i_r_end=i_r_start+20000
            s_l_seq = f_fa.fetch(chrm, i_l_start, i_l_end)
            fout.write(">"+chrm+"_left\n")
            fout.write(s_l_seq+"\n")

            s_r_seq = f_fa.fetch(chrm, i_r_start, i_r_end)
            fout.write(">" + chrm + "_right\n")
            fout.write(s_r_seq + "\n")
        f_fa.close()
####
# #def get_hit_pos():
# samfile = pysam.AlignmentFile(sf_algnmt, s_open_fmt)  # read in the sam file
# for algnmt in samfile.fetch():  # check each alignment, and find "left" and "right" flank #