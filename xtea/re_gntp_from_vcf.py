import os
import sys
from x_genotype_classify import *
#this is a separate module for re-genotyping the called insertions from the final VCF
#The goal is to correct some wrongly genotyped TE insertions that predicted from an old version of ML model

#Input is a list of VCFs, and each row with two columns: ID full-path-VCF

def re_gntp_vcf(sf_model, sf_vcf, sf_new_vcf):
    gntp_classifer=GntpClassifier_DF21()
    gntp_classifer.predict_for_site_vcf(sf_model, sf_vcf, sf_new_vcf)

def re_gntp_vcfs(sf_model, sf_list, sf_out_prefix):
    if len(sf_out_prefix)<=0:
        sf_out_prefix="./"
    elif sf_out_prefix[-1]!="/":
        sf_out_prefix+="/"
    with open(sf_list) as fin_list:
        for line in fin_list:
            l_fields=line.rstrip().split()
            sf_vcf=l_fields[1]
            s_id=l_fields[0]####
            sf_new_vcf=sf_out_prefix+s_id+"_re_gnpted.vcf"
            re_gntp_vcf(sf_model, sf_vcf, sf_new_vcf)

if __name__ == '__main__':
    sf_model=sys.argv[1]
    sf_vcf_list=sys.argv[2]
    sf_out_prefix=sys.argv[3]
    re_gntp_vcfs(sf_model, sf_vcf_list, sf_out_prefix)
####