import sys
import os
from subprocess import *
from multiprocessing import Pool
from Bio import SeqIO

def run_cmd(cmd):
    Popen(cmd, shell = True, stdout = PIPE).communicate()

def tabix_sv(sf_fai, s_individual):
    #grep the header
    # cmd="tabix -H ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20150224_integrated_sv_map/v0/ALL.wgs.integrated_sv_map.20130502.svs.genotypes.vcf.gz " \
    #     ">> {0}_SV.vcf".format(s_individual)
    # run_cmd(cmd)
    bhead=True
    sf_sv="{0}_SV.vcf".format(s_individual)
    with open(sf_fai) as fin_fai:
        for line in fin_fai:
            fields=line.split()
            chrm=fields[0]
            lenth=int(fields[1])

            sf_chrm_sv="{0}_{1}_SV.vcf".format(s_individual, chrm)
            s_reg="{0}:0-{1}".format(chrm, lenth+1000)
            cmd="tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20150224_integrated_sv_map/v0/ALL.wgs.integrated_sv_map.20130502.svs.genotypes.vcf.gz " \
                "{0} | perl /data2/chongchu/tools/vcftools-vcftools-1d27c24/src/perl/vcf-subset -c {1} " \
                "> {2}".format(s_reg, s_individual, sf_chrm_sv)
            print cmd
            run_cmd(cmd)
            if bhead==True:
                bhead=False
                cmd="mv {0} {1}".format(sf_chrm_sv, sf_sv)
                run_cmd(cmd)
            else:
                with open(sf_chrm_sv)as fin_chrm_sv, open(sf_sv,"a") as fout_sv:#
                    for line in fin_chrm_sv:
                        if line[0]!="#":
                            fout_sv.write(line)

if __name__ == "__main__":
    sf_fai=sys.argv[1]
    s_individual=sys.argv[2]
    tabix_sv(sf_fai, s_individual)