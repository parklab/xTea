#This is to do snapshot for given sites and bams using BamSnap
class XBamSnap():
    def prepare_bamsnap_scripts_multi_bams(self, sf_sites, sf_bams, s_wfolder, i_extnd, sf_ref, sf_out):#
        if len(s_wfolder)<=0:
            s_wfolder="./"
        elif s_wfolder[-1]!="/":
            s_wfolder+="/"

        m_sites = {}
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields = line.split()
                chrm = fields[0]
                start = fields[1]
                end = str(int(start)+1)
                s_site = chrm + ":" + start + "-" + end
                s_title = chrm + "_" + start + "_" + end
                m_sites[s_site] = s_title

        ####
        m_bams = {}
        with open(sf_bams) as fin_bams:
            for line in fin_bams:
                fields = line.split()
                if len(fields) < 2:
                    print(("Wrong bam file list ", line))
                    continue
                l_sf_bams = []
                for sf_tmp_bam in fields[1:]:
                    l_sf_bams.append(sf_tmp_bam)
                m_bams[fields[0]] = l_sf_bams

        ####
        with open(sf_out, "w") as fout_cmd:
            for s_site in m_sites:
                s_title = m_sites[s_site]#
                for s_sample in m_bams:
                    sf_bam = " ".join(m_bams[s_sample])
                    sf_name=s_wfolder + s_sample + "_" + s_title + ".png"
                    cmd = "bamsnap -bam %s -title %s -pos %s -bamplot coverage read -margin %d -no_target_line " \
                          "-show_soft_clipped -read_color_by interchrom -save_image_only -ref %s -out %s\n" % \
                          (sf_bam, s_sample + "_" + s_title, s_site, i_extnd, sf_ref, sf_name)
                    sf_name2 = s_wfolder + s_sample + "_" + s_title + "_for_transfer_learning.png"
                    cmd2 = "bamsnap -bam %s -pos %s -draw bamplot -bamplot read -margin 400 -plot_margin_top 0 " \
                           "-plot_margin_bottom 0 -separator_height 5 -no_title " \
                          "-show_soft_clipped -read_color_by interchrom -save_image_only -ref %s -out %s\n" % \
                          (sf_bam, s_site, sf_ref, sf_name2)
                    fout_cmd.write(cmd)
                    fout_cmd.write(cmd2)

                    l_trio_id=["proband","father","mother"]##
                    i_idx=0
                    if len(m_bams[s_sample])<=3:#for de novo
                        for sf_tmp_bam in m_bams[s_sample]:
                            s_tmp_id=l_trio_id[i_idx]#this is for each individual
                            sf_name_tmp = s_wfolder + s_sample + "_" + s_title + "_%s_for_transfer_learning.png" % (s_tmp_id)
                            cmd_tmp = "bamsnap -bam %s -pos %s -draw bamplot -bamplot read -margin 400 -plot_margin_top 0 " \
                                   "-plot_margin_bottom 0 -separator_height 5 -no_title " \
                                   "-show_soft_clipped -read_color_by interchrom -save_image_only -ref %s -out %s\n" % \
                                   (sf_tmp_bam, s_site, sf_ref, sf_name_tmp)
                            fout_cmd.write(cmd_tmp)
                            i_idx+=1
####
    def prepare_bamsnap_scripts_multi_processed_bams(self, sf_sites, sf_bams, s_wfolder, i_extnd, sf_ref, sf_out):#
        if len(s_wfolder)<=0:
            s_wfolder="./"
        elif s_wfolder[-1]!="/":
            s_wfolder+="/"

        m_sites = {}
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields = line.split()
                chrm = fields[0]
                start = fields[1]
                end = str(int(start)+1)
                s_site = chrm + ":" + start + "-" + end
                s_title = chrm + "_" + start + "_" + end
                m_sites[s_site] = s_title

        ####
        m_bams = {}
        with open(sf_bams) as fin_bams:
            for line in fin_bams:
                fields = line.split()
                if len(fields) < 2:
                    print(("Wrong bam file list ", line))
                    continue
                l_sf_bams = []
                for sf_tmp_bam in fields[1:]:
                    l_sf_bams.append(sf_tmp_bam)
                m_bams[fields[0]] = l_sf_bams

        ####
        with open(sf_out, "w") as fout_cmd:
            for s_site in m_sites:
                s_title = m_sites[s_site]#
                for s_sample in m_bams:
                    sf_bam = " ".join(m_bams[s_sample])
                    sf_name=s_wfolder + s_sample + "_" + s_title + ".processed.png"
                    cmd = "bamsnap -bam %s -title %s -pos %s -bamplot coverage read -margin %d -height 500 -no_target_line " \
                          "-show_soft_clipped -read_color_by interchrom -save_image_only -ref %s -out %s\n" % \
                          (sf_bam, s_sample + "_" + s_title, s_site, i_extnd, sf_ref, sf_name)
                    sf_name2 = s_wfolder + s_sample + "_" + s_title + "_for_transfer_learning.processed.png"
                    cmd2 = "bamsnap -bam %s -pos %s -draw bamplot -bamplot read -margin 400 -height 500 -plot_margin_top 0 " \
                           "-plot_margin_bottom 0 -separator_height 150 -no_title " \
                          "-show_soft_clipped -read_color_by interchrom -save_image_only -ref %s -out %s\n" % \
                          (sf_bam, s_site, sf_ref, sf_name2)
                    fout_cmd.write(cmd)
                    fout_cmd.write(cmd2)

                    l_trio_id=["proband","father","mother"]##
                    i_idx=0
                    if len(m_bams[s_sample])<=3:#for de novo
                        for sf_tmp_bam in m_bams[s_sample]:
                            s_tmp_id=l_trio_id[i_idx]#this is for each individual
                            sf_name_tmp = s_wfolder + s_sample + "_" + s_title + "_%s_for_transfer_learning.processed.png" % (s_tmp_id)
                            cmd_tmp = "bamsnap -bam %s -pos %s -draw bamplot -bamplot read -margin 400 -plot_margin_top 0 " \
                                   "-plot_margin_bottom 0 -separator_height 5 -no_title " \
                                   "-show_soft_clipped -read_color_by interchrom -save_image_only -ref %s -out %s\n" % \
                                   (sf_tmp_bam, s_site, sf_ref, sf_name_tmp)
                            fout_cmd.write(cmd_tmp)
                            i_idx+=1