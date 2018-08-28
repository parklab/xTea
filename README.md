# xTEA

## Dependency
1. bwa (version 0.7 or later), which can be downloaded from https://github.com/lh3/bwa.
2. samtools (v1.0 or later), which can be downloaded from https://github.com/samtools.
3. Python 2.7 or later version
	+ pysam (https://github.com/pysam-developers/pysam, v0.12 or later) is required to be installed.
	
		+ In detail, first install Anaconda:
		```
		wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh
		sh Anaconda2-5.2.0-Linux-x86_64.sh
		```
		
		+ Install pysam:
		```
		conda config --add channels r
		conda config --add channels bioconda
		conda install pysam
		```
4. Note: bwa and samtools need to be added to the $PATH.

## Run xTEA
1. **Input**
	+ A sample id file, e.g. a file named `sample_id.txt` with content:
	
		```
		NA12878
		NA12877
		```
	
	+ An Illumina bam file (sorted and indexed) list, e.g. a file named `illumina_bam_list.txt` with content:

		```
		NA12878 /path/na12878_illumina_1_sorted.bam
		NA12878 /path/na12878_illumina_2_sorted.bam
		NA12877 /path/na12877_illumina_1_sorted.bam
		NA12877 /path/na12877_illumina_2_sorted.bam
		```
	
	+  An 10X bam file (sorted and indexed, find [TXtools](https://github.com/parklab/TXtools) for barcode based index) list, e.g. a file named `10X_bam_list.txt` with content:
	
		```
		NA12878 /path/na12878_10X_1_sorted.bam /path/na12878_10X_1_barcode_indexed.bam
		NA12878 /path/na12878_10X_2_sorted.bam /path/na12878_10X_2_barcode_indexed.bam
		NA12877 /path/na12877_10X_1_sorted.bam /path/na12877_10X_1_barcode_indexed.bam
		NA12877 /path/na12877_10X_2_sorted.bam /path/na12877_10X_2_barcode_indexed.bam
		```
		
2. **Run the pipeline**
		
	2.1 Generate the running script.	

	
	+ Run on Amazon AWS
		+ Run for one given bam/cram file
		```
		python ./xTEA/gnrt_pipeline_cloud.py -b input.bam -p /home/ec2-user/results2/ -o run_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -x /home/ec2-user/xTEA/ --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 5 --flklen 3000 -f 19 -y 7
		```
		
		+ Parameters for Amazon AWS
		
		```
			-D: if set, then decompress the files(.tar.gz) -l and -r specify			-b: the input bam/cram file (sorted and indexed);
			-p: working folder, where the results and temporary files will be saved;
			-o: temporary running scripts under the working folder;
			-n: number of cores
			-l: repeat library folder (if -D is set, then this is directly the `rep_lib_annotation.tar.gz` file) (decompressed from `s3://leelab-datafiles/rep_lib_annotation.tar.gz`);
			-r: reference genome file (-f -D is set, then this is the `hg19_decoy.tar.gz`) (indexed by `bwa index`, and decompressed from `s3://leelab-datafiles/hg19_decoy.tar.gz`)
			-x: path of xTEA folder (xTEA can be downloaded with command `git clone https://github.com/parklab/xTEA`)
			-y: type of repeats will work on (1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrion)
			Other parameters can be keep unchanged.
		```
		
	+ Run on O2 (slurm) cluster
		+ Only with Illumina data
			```
			python ./xTEA/gnrt_pipeline_local.py -i sample_id.txt -b illumina_bam_list.txt -x null -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -x /home/ec2-user/xTEA/ --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 5 --flklen 3000 -f 19 -y 1 			```

		+ Only with 10X data
			```
				python ./xTEA/gnrt_pipeline_local.py -i sample_id.txt -b null -x 10X_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -x /home/ec2-user/xTEA/ --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 5 --flklen 3000 -f 19 -y 1 			```
		
		+ Working with hybrid data of 10X and Illumina
			```
			python ./xTEA/gnrt_pipeline_local.py -i sample_id.txt -b illumina_bam_list.txt -x 10X_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -x /home/ec2-user/xTEA/ --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 5 --flklen 3000 -f 19 -y 1
			```
		+ Parameters:
			
			```
				-i: samples id list file (each sample id per line);
				-b: Illumna bam/cram file list file (sorted and indexed, each file per line);
				-x: 10X bam file list (sorted and indexed, each file per line);
				-p: working folder, where the results and temporary files will be saved;
				-o: temporary running scripts under the working folder;
				-n: number of cores;
				-l: repeat library folder;	
				-r: reference genome file;
				-y: type of repeats will work on (1-L1, 2-Alu, 4-SVA, 8-HERV, 16-Mitochondrion)
				Other parameters can be keep unchanged or adjust according.
			```
		
	2.2 The previous step will generate a shell script called `run_xTEA_pipeline.sh` under `WFOLDER/sample_id/L1(or other type of repeats)`.
		
	+ To run on the cluster, like slurm system: `sbatch < run_xTEA_pipeline.sh`
	+ To run on a single node: `sh run_xTEA_pipeline.sh`
	
3. **Output**

	3.1 For Illumina, `candidate_disc_filtered_cns.txt` is the final output.
	