## xTea

xTea (comprehensive Transposable element analyzer) is designed to identify TE insertions from paired-end Illumina reads, barcode Linked-Reads, long reads (PacBio or Nanopore), or hybrid data (WGS/WES). 

Copyright (c) 2019 - President and Fellows of Harvard College. All rights reserved.

##Download
1. short reads (Illumina and Linked-Reads)

	+ 1.1 Latest version

	```
	git clone https://github.com/parklab/xTea.git
	```

	+ 1.2 cloud binary version

	```
	git clone --single-branch --branch release_xTea_cloud_1.0.0-beta  https://github.com/parklab/xTea.git
	```

2. long reads (PacBio or Nanopore)

	```
	git clone --single-branch --branch xTea_long_release_v0.1.0 https://github.com/parklab/xTea.git
	```
	
## Dependency
1. bwa (version **0.7.17** or later, require the **-o** option), which can be downloaded from https://github.com/lh3/bwa.
2. samtools (v1.0 or later), which can be downloaded from https://github.com/samtools.
3. minimap2 (for long reads only), which can be downloaded from https://github.com/lh3/minimap2.
4. wtdbg2 (for long reads only), which can be downloaded from https://github.com/ruanjue/wtdbg2.
5. Python 2.7 or later version
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
			conda install pysam -y
			```
	+ sortedcontainers
		+ Install sortedcontainers
		`conda install sortedcontainers -y`

4. Note: bwa and samtools need to be added to the $PATH.

## Run xTea
1. **Input**
	+ A sample id file, e.g. a file named `sample_id.txt` with content:
	
		```
		NA12878
		NA12877
		```
	
	+ An Illumina bam file (sorted and indexed) list, e.g. a file named `illumina_bam_list.txt` with content:

		```
		NA12878 /path/na12878_illumina_1_sorted.bam
		NA12877 /path/na12877_illumina_1_sorted.bam
		```
	
	+  A 10X bam file (sorted and indexed, find [TXtools](https://github.com/parklab/TXtools) for barcode based index) list, e.g. a file named `10X_bam_list.txt` with content:
	
		```
		NA12878 /path/na12878_10X_1_sorted.bam /path/na12878_10X_1_barcode_indexed.bam
		NA12877 /path/na12877_10X_1_sorted.bam /path/na12877_10X_1_barcode_indexed.bam
		```
		
2. **Run the pipeline**
	
	
	2.1 Generate the running script.	
			
	+ Run on a slurm cluster or a single node (Please replace `gnrt_pipeline_local.pyc` with `gnrt_pipeline_local_v38.pyc` for running on **GRCh38**)
		+ Only with Illumina data
			```
			python ./xTEA/gnrt_pipeline_local.py -i sample_id.txt -b illumina_bam_list.txt -x null -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa --xtea /home/ec2-user/xTEA/ --flklen 3000 -y 7 			```

		+ Only with 10X data
			```
			python ./xTEA/gnrt_pipeline_local.py -i sample_id.txt -b null -x 10X_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -xtea /home/ec2-user/xTEA/ --flklen 3000 -y 7 			```
		
		+ Working with hybrid data of 10X and Illumina 
			```
			python ./xTEA/gnrt_pipeline_local.pyc -i sample_id.txt -b illumina_bam_list.txt -x 10X_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -xtea /home/ec2-user/xTEA/ --flklen 3000 -y 7
			```
			
		+ Parameters:
			
			```
			Required:
				-i: samples id list file (each sample id per line);
				-b: Illumna bam/cram file list file (sorted and indexed, each file per line);
				-x: 10X bam file list (sorted and indexed, each file per line);
				-p: working folder, where the results and temporary files will be saved;
				-l: repeat library folder;
				-r: reference genome file;
				-y: type of repeats will work on (1-L1, 2-Alu, 4-SVA, 8-HERV; sum all selected as one value)
				--xtea: xTEA full path
				-o: generated running scripts under the working folder;
			Optional:
				-n: number of cores (by default 8, should be an integer);
				-m: maximum memory in GB (by default 25, should be an integer);
				
			Other parameters can be keep unchanged or adjust accordingly.
			```
		
	2.2 The previous step will generate a shell script called `run_xTEA_pipeline.sh` under `WFOLDER/sample_id/L1(or other type of repeats)`, where `WFOLDER` is specified by `-p` option.
		
	+ To run on the cluster, like slurm system: `sbatch < run_xTEA_pipeline.sh`
	+ To run on a single node: `sh run_xTEA_pipeline.sh`
	
	
	2.3 Run from the Cloud
	
	+ A docker file and cwl file are provided for running on AWS/GCP/Fire Cloud.
	
		
	
3. **Output**

	A gVCF file will be generated for each sample.


## Contact
Requests for use of the Software for or on behalf of for-profit entities or for any commercial purposes, please contact:

Office of Technology Development  
Harvard University  
Smith Campus Center, Suite 727E  
1350 Massachusetts Avenue  
Cambridge, MA 02138 USA  
Telephone: (617) 495-3067  
E-mail: otd@harvard.edu