import os
import sys
from subprocess import *
from multiprocessing import Pool

#
def run_extract(record):
    sf_file = record
    sf_out = record+".fastq"
    #cmd = "cat {0} | xargs dextract -faq -o./{1}".format(sf_file, sf_out) #output in one file
    cmd = "cat {0} | xargs dextract -faq ".format(sf_file)
    Popen(cmd, shell=True, stdout=PIPE).communicate()

####merge the files
def merge_files(file_suffix, sf_prefix):
    cmd = "cat *.{0} > {1}.{2}".format(file_suffix, sf_prefix, file_suffix)
    Popen(cmd, shell=True, stdout=PIPE).communicate()

# split the list to several partial lists
def split_file_list(sf_list, n_jobs, sf_out_prefix):
    m_list = {}
    with open(sf_list) as fin_list:
        for line in fin_list:
            path_fields = line.split("/")
            fields = path_fields[-1].split(".")
            m_list[fields[0]] = 1
    # split the list to nparts
    l_list = m_list.keys()
    n_all = len(m_list)
    n_each_part = n_all / n_jobs
    l_file_list=[]
    for i in range(n_jobs):
        istart = i * n_each_part
        iend = istart + n_each_part
        l_part_list = l_list[istart:iend]
        if i == n_jobs - 1:
            l_part_list = l_list[istart:]
        sf_partial = sf_out_prefix + "_p{0}".format(i)
        with open(sf_partial, "w") as fout_partial:
            for hdf5_prefix in l_part_list:
                s1=hdf5_prefix + ".1.bax.h5"
                if os.path.isfile(s1):
                    fout_partial.write(s1 + "\n")
                else:
                    print s1, "doesn't exist!!!!"
                s2=hdf5_prefix + ".2.bax.h5"
                if os.path.isfile(s2):
                    fout_partial.write(s2 + "\n")
                else:
                    print s2, "doesn't exist!!!!"
                s3=hdf5_prefix + ".3.bax.h5"
                if os.path.isfile(s3):
                    fout_partial.write(s3 + "\n")
                else:
                    print s3, "doesn't exist!!!!"
        l_file_list.append(sf_partial)

    pool = Pool(n_jobs)
    pool.map(run_extract, l_file_list, 1)
    pool.close()
    pool.join()

    merge_files("fasta", sf_out_prefix)
    merge_files("quiva", sf_out_prefix)
    merge_files("arrow", sf_out_prefix)

    # #merge the fasta files
    # s_files="{0}_p*.fasta".format(sf_out_prefix)
    # cmd="cat {0} > {1}.fasta".format(s_files, sf_out_prefix)
    # Popen(cmd, shell=True, stdout=PIPE).communicate()
    # #merge the quiva files
    # s_files = "{0}_p*.quiva".format(sf_out_prefix)
    # cmd = "cat {0} > {1}.quiva".format(s_files, sf_out_prefix)
    # Popen(cmd, shell=True, stdout=PIPE).communicate()
    # #merge the arrow files
    # s_files = "{0}_p*.arrow".format(sf_out_prefix)
    # cmd = "cat {0} > {1}.arrow".format(s_files, sf_out_prefix)
    # Popen(cmd, shell=True, stdout=PIPE).communicate()


if __name__ == '__main__':
    sf_list=sys.argv[1]
    n_jobs=int(sys.argv[2])
    sf_out_prefix=sys.argv[3]
    split_file_list(sf_list, n_jobs, sf_out_prefix)

