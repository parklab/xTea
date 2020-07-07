import os
import sys
from subprocess import *
from multiprocessing import Pool

def run_cvt(cmd):
    Popen(cmd, shell=True, stdout=PIPE).communicate()

def parse_bam(sf_list, n_jobs):
    with open(sf_list) as fin_list:
        l_cmd=[]
        for line in fin_list:
            fields=line.split(".")
            sf_sra=line.rstrip()
            sf_out=fields[0]+".bam"
            cmd="sam-dump {0} | samtools view -bS - > {1}".format(sf_sra, sf_out)
            l_cmd.append(cmd)

        pool = Pool(n_jobs)
        pool.map(run_cvt, l_cmd, 1)
        pool.close()
        pool.join()

if __name__ == '__main__':
    sf_list=sys.argv[1]
    n_jobs=int(sys.argv[2])
    parse_bam(sf_list, n_jobs)