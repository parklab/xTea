import os
import sys
from subprocess import *
from multiprocessing import Pool

HG19="/n/data1/hms/dbmi/park/SOFTWARE/LongRanger/refdata-b37-2.1.0/fasta/genome.fa"
GRCh38="/n/data1/hms/dbmi/park/simon_chu/projects/data/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
def run_align(record):
    sf_fa=record[0]
    s_id=record[1]
    sf_ref=record[2]
    #
    sf_out=s_id+"_2_ref.bam"
    #align, here use 8 cores for alignment
    cmd="minimap2 -t 8 -ax map-pb {0} {1} | samtools view -h -S -b - > {2}".format(sf_ref, sf_fa, sf_out)
    Popen(cmd, shell=True, stdout=PIPE).communicate()
    #sort
    sf_tmp="./tmp"
    if os.path.exists(sf_tmp)==False:
        cmd="mkdir {0}".format(sf_tmp)
        Popen(cmd, shell=True, stdout=PIPE).communicate()
    cmd="sambamba sort -m 45G -t 8 --tmpdir={0} {1}".format(sf_tmp, sf_out)
    Popen(cmd, shell=True, stdout=PIPE).communicate()

def align_long_read(sf_list, sf_ref, n_jobs):
    l_list=[]
    with open(sf_list) as fin_list:
        for line in fin_list:
            fields=line.split()
            sf_fa=fields[0]
            s_id=fields[1]
            l_list.append((sf_fa, s_id, sf_ref))
    ####
    pool = Pool(n_jobs)
    pool.map(run_align, l_list, 1)
    pool.close()
    pool.join()

if __name__ == '__main__':
    sf_list=sys.argv[1]
    sf_ref=sys.argv[2]
    n_jobs=int(sys.argv[3])
    align_long_read(sf_list, sf_ref, n_jobs)