# Simple script for generate sh scripts to run longranger wgs
# see: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/using/wgs
# By Carl Vitzthum

import argparse
import os

# global variables
ref_path = '/n/data1/hms/dbmi/park/SOFTWARE/LongRanger/refdata-b37-2.1.0'
software_path = '/n/data1/hms/dbmi/park/SOFTWARE/LongRanger/longranger-2.1.4/longranger'
# 0 is female, 1 is male
sex_by_sample = {
    '1465QA25_A12_S6_L006': 1,
    '1465QA25_B08_S4_L004': 1,
    '1465QA25_D09_S5_L005_R1_001': 1,
    '1465_S1_L001': 1,
    '4638_S2_L002': 0,
    '4643_S3_L003': 0
}

def walk_fastq_dir(fastq_dir):
    samples = []
    for root, dirs, files in os.walk(fastq_dir, topdown=True):
        for fname in files:
            # IMPORTANT: these slices may need to changed based on filenames
            if fname[-9:] != '.fastq.gz':
                continue
            sample_id = fname[:-24]  # eliminate _S#_L00#_R#_001.fastq.gz for sample name
            if sample_id not in samples:
                samples.append(sample_id)
    print('>>SAMPLES>>\n', samples)
    return samples


def one_script(fastq_dir, sample, cores, mem, sex, somatic):
    """
    Generate a single sh script for the longranger wgs command for one sample
    """
    write_list = ['#!/bin/bash']
    lr_str = (software_path + ' wgs --id=' + sample + ' --reference=' + ref_path +
        ' --sample=' + sample + ' --fastqs=' + fastq_dir + ' --localcores=' +
        str(cores) + ' --localmem=' + str(mem))
    if sex and sample in sex_by_sample:
        if sex_by_sample[sample] == 1:
            lr_str += ' --sex=male'
        else:
            lr_str += ' --sex=female'
    if somatic:
        lr_str += ' --somatic'
    write_list.append(lr_str)
    write_out('longranger_wgs_'+sample+'.sh', write_list)
    o_str = ('bsub -q mcore -n ' + str(cores) + ' -R "rusage[mem='
             + str(mem*1028) + ']"'
             ' -o longranger_wgs_'+sample+'.lsf -N -W 240:00 "bash ./longranger_wgs_'+sample+'.sh"')
    return o_str


def write_out(filename, w_list):
    with open(filename, 'w') as f:
        for line in w_list:
            f.write(line + '\n')


def build_scripts(fastq_dir, cores, mem, sex, somatic):
    samples = walk_fastq_dir(fastq_dir)
    master_wl = ['#!/bin/bash']
    for sample in samples:
        master_wl.append(one_script(fastq_dir, sample, cores, mem, sex, somatic))
    write_out('longranger_wgs_master.txt', master_wl)


def main():
    """
    Execute the program from the command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_dir", help="Path to directory of fastq files")
    parser.add_argument("-c", "--cores",help="Number of cores to use. Defaults to 1", type=int, default=1)
    parser.add_argument("-m", "--mem",help="Amount of memory to use PER CORE. Defaults to 4G", type=int, default=4)
    parser.add_argument("-s", "--sex", help="If used, use internal sex_by_sample variable to determine --sex parameter. Else, let wgs detect automatically", action="store_true")
    parser.add_argument("-o", "--somatic", help="If used, sets the --somatic parameter for longranger wgs", action="store_true")
    args = parser.parse_args()
    build_scripts(args.fastq_dir, args.cores, args.mem, args.sex, args.somatic)


if __name__ == '__main__':
    main()

