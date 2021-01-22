## This is a pre-release version of xTea for paper review only. This software is provided ``as is`` without warranty of any kind. We will formally relese xTea very soon. Before that, please contact us if you have any questions or want to collaborate.


## xTea

xTea (comprehensive Transposable element analyzer) is designed to identify TE insertions from paired-end Illumina reads, barcode Linked-Reads, long reads (PacBio or Nanopore), or hybrid data (WGS/WES). 

![alt text](./xTea_workflow.png)


## Download

1. short reads (Illumina and Linked-Reads)

	+ 1.1 Latest version (master branch)

	```
	git clone https://github.com/parklab/xTea.git
	```

	+ 1.2 cloud binary version (branch: release_xTea_cloud_1.0.0-beta)

	```
	git clone --single-branch --branch release_xTea_cloud_1.0.0-beta  https://github.com/parklab/xTea.git
	```

2. long reads (PacBio or Nanopore, branch to be merged: xTea_long_release_v0.1.0)

	```
	git clone --single-branch --branch xTea_long_release_v0.1.0 https://github.com/parklab/xTea.git
	```
3. pre-processed repeat library used by xTea (The same library is used for both short and long reads)  
	The file size is large, and please use git lfs (https://git-lfs.github.com/)  
	```
	git lfs get 
	```
	
	Or directly download through `wget https://github.com/parklab/xTea/raw/master/rep_lib_annotation.tar.gz`. 
4. Gene annotation file are downloaded from GENCODE (https://www.gencodegenes.org/human/release_33.html)

## Dependency

1. bwa (version **0.7.17** or later, require the **-o** option), which can be downloaded from https://github.com/lh3/bwa.
2. samtools (v1.0 or later), which can be downloaded from https://github.com/samtools.
3. minimap2 (for long reads only), which can be downloaded from https://github.com/lh3/minimap2.
4. wtdbg2 (for long reads only), which can be downloaded from https://github.com/ruanjue/wtdbg2.
5. Python 2.7+/3.6+
	+ pysam (https://github.com/pysam-developers/pysam, v0.12 or later) is required to be installed.
		+ Install pysam:

			```
			conda config --add channels r
			conda config --add channels bioconda
			conda install pysam -y
			```
	+ sortedcontainers
		+ Install sortedcontainers
		`conda install sortedcontainers -y`

	+ numpy, scikit-learn, and pandas
			+ Install numpy, scikit-learn and pandas
			`conda install numpy scikit-learn pandas -y`

6. Note: bwa and samtools need to be added to the $PATH.


## Install

1. **Use Conda**
	```
	conda install xtea
	```
2. **Install free**
	```
	git clone https://github.com/parklab/xTea
	```


## Run xTea
1. **Input**
	+ A sample id file, e.g. a file named `sample_id.txt` with content (each line is one unique sample id):
	
		```
		NA12878
		NA12877
		```
	
	+ A file of listed alignments:

		+ An Illumina bam/cram file (sorted and indexed) list, e.g. a file named `illumina_bam_list.txt` with content (two columns separated by space or tab: sample-id bam-path):

			```
			NA12878 /path/na12878_illumina_1_sorted.bam
			NA12877 /path/na12877_illumina_1_sorted.bam
			```
		
		+  A 10X bam/cram file (sorted and indexed, find [BarcodeMate](https://github.com/simoncchu/BarcodeMate) for barcode based index) list, e.g. a file named `10X_bam_list.txt` with content (three columns separated by space or tab: sample-id bam-path barcode-index-bam-path):
		
			```
			NA12878 /path/na12878_10X_1_sorted.bam /path/na12878_10X_1_barcode_indexed.bam
			NA12877 /path/na12877_10X_1_sorted.bam /path/na12877_10X_1_barcode_indexed.bam
			```
		
		+  A list of cast-ctrl bam/cram files. Three columns per line separated by space or tab: sample-id case-bam-path ctrl-bam-path
			```
			DO0001 /path/DO001_case_sorted.bam /path/DO001_ctrl_sorted.bam
			DO0002 /path/DO002_case_sorted.bam /path/DO002_ctrl_sorted.bam
			```


2. **Run the pipeline**
	
	
	2.1 Generate the running script (if it is install-free, then use `bin/xtea` instead.)
			
	+ Run on a slurm cluster or a single node (by default `xtea` assume the reference genome is **GRCh38** or **hg38**, for `hg19` or `GRCh37`, please use `xtea_hg19`)
		+ Only with Illumina data
			```
			xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/rep_lib_annotation/ -r /home/reference/genome.fa -g /home/gene_annotation_file.gff3 --xtea /home/ec2-user/xtea/ -f 5907 -y 7 			```

		+ Only with 10X data
			```
			xtea -i sample_id.txt -b null -x 10X_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -g /home/gene_annotation_file.gff3  --xtea /home/ec2-user/xtea/ -y 7 -f 5907 			```
		
		+ Working with hybrid data of 10X and Illumina 
			```
			xtea -i sample_id.txt -b illumina_bam_list.txt -x 10X_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -g /home/gene_annotation_file.gff3 --xtea /home/ec2-user/xtea/ -y 7 -f 5907
			```
		+ Working with case-ctrl mode
			```
			xtea --case_ctrl --tumor -i sample_id.txt -b case_ctrl_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -g /home/gene_annotation_file.gff3 --xtea /home/ec2-user/xtea/ -y 7 -f 5907
			```
		+ Working with long reads (non case-ctrl; more detailed steps please check the "xTea_long_release_v0.1.0" branch)
			```
			xtea_long -i sample_id.txt -b long_read_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -n 8 -m 32 --rmsk ./rep_lib_annotation/LINE/hg38/hg38_L1_larger2K_with_all_L1HS.out -r /home/ec2-user/reference/genome.fa --cns ./rep_lib_annotation/consensus/LINE1.fa --rep /home/ec2-user/rep_lib_annotation/ -f 31 -y 7 --clean
			```

		+ Parameters:
			
			```
			Required:
				-i: samples id list file (each sample id per line);
				-b: Illumna bam/cram file list file (sorted and indexed, each file per line);
				-x: 10X bam file list (sorted and indexed, each file per line);
				-p: working folder, where the results and temporary files will be saved;
				-l: repeat library folder (folder contain files decompressed from the downloaded "rep_lib_annotation.tar.gz");
				-r: reference genome fasta/fa file;
				-y: type of repeats will work on (1-L1, 2-Alu, 4-SVA, 8-HERV; sum all selected as one value. 
				    For example, if want to check L1 and SVA only, then set `-y 5`. 
				    Each repeat type will be separately processed, however some of the early steps are shared. 
				    Thus, if the user has a large cohort, to improve the efficiency (and save money on cloud), 
				    we highly recommend to run on one repeat type first, and then on the rest. 
				    For example, first set '-y 1', and for the second run set '-y 6');
				-f: steps to run. (5907 means run all the steps);
				--xtea: this is the full path of the xTea/xtea folder (or the xTea_long_release_v0.1.0 folder for long reads module), 
				        where the python scripts reside in.
				-g: gene annotation file in gff3 format;
				-o: generated running scripts under the working folder;
			Optional:
				-n: number of cores (by default 8, should be an integer);
				-m: maximum memory in GB (by default 25, should be an integer);
				-q: partition name for a cluster;
				-t: job running time;
				--flklen: flanking region length;
				--lsf: add this option if this is for LSF system (by default slurm system);
				--tumor: indicates tumor case-ctrl samples;
				--purity: tumor purity (by default 0.45);
				--blacklist: black list file in bed format, and candidates fall in the regions will be filtered out;
				--slurm: add script header for slurm system;
				--lsf: add script header for LSF system
			
			Cutoffs will be automatically set based on the read depth (and also the purity if it is a tumor sample); 
			parameters have been thoroughly tuned based on the test on benchmark data and also on large cohort analysis. 
			For advanced users (optional major cutoffs):
				--user: by default, this is turned off. If this option is set, then user specific cutoff will be used;
				--nclip: minimum number of clipped reads;
				--cr: minimum number of clipped reads whose mates map in repetitive regions;
				--nd: minimum number of discordant pair;

			Specific parameters for long reads module:
			    --rmsk: this is reference full length L1 annotation file from RepeatMasker only for the "ghost" L1 detection module. 
			            One file named "hg38_L1_larger2K_with_all_L1HS.out" within the downloaded library could be directly used;
			    --cns: this is the L1 concensus sequence needed only by the "ghost" L1 detection module. 
			           One file named "LINE1.fa" within the downloaded library could be directly used;
			    --rep: repeat library folder (folder contain files decompressed from the downloaded "rep_lib_annotation.tar.gz");
			    --clean: clean the intermediate files

			```
		
	2.2 The previous step will generate a shell script called `run_xTea_pipeline.sh` under `WFOLDER/sample_id/L1(or other type of repeats)`, where `WFOLDER` is specified by `-p` option.
		
	+ To run on the script: `sh run_xTea_pipeline.sh` or users can submit the jobs (each line one job) to a cluster.
	
	
	2.3 Run from the Cloud
	
	+ A docker file and cwl file are provided for running on AWS/GCP/Fire Cloud.
	
			
3. **Output**

	A gVCF file will be generated for each sample.
	+ For germline TE insertion calling on short reads, the `orphan transduction` module usually has higher false positive rate, thus users could filter out those events with command like `grep -v "orphan" out.vcf > new_out.vcf` to get those higher confident events.

