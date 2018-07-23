# xTEA

## Dependency
1. bwa (version 0.7 or later), which can be downloaded from https://github.com/lh3/bwa.
2. samtools (v1.0 or later), which can be downloaded from https://github.com/samtools.
3. Python 2.7 or later version
	+ pysam (https://github.com/pysam-developers/pysam, v0.13 or later) is required to be installed.
	
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

	+ Only with Illumina data
		```
		python ./../xTEA/gnrt_pipeline_cloud.py -b input.bam -p /home/ec2-user/results2/ -o run_jobs.sh -n 8 -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -x /home/ec2-user/xTEA/ --nclip 4 --cr 2 --nd 5 --nfclip 3 --nfdisc 5 --flklen 3000 -f 19
		```

			+ For the parameters:
			```
				-D: if set, then decompress the files(.tar.gz) -l and -r specify
				-b: the input bam/cram file (sorted and indexed);
				-p: working folder, where the results and temporary files will be saved;
				-o: temporary running scripts under the working folder;
				-n: number of cores
				-l: repeat library folder (if -D is set, then this is directly the `rep_lib_annotation.tar.gz` file) (decompressed from `s3://leelab-datafiles/rep_lib_annotation.tar.gz`);
				-r: reference genome file (-f -D is set, then this is the `hg19_decoy.tar.gz`) (indexed by `bwa index`, and decompressed from `s3://leelab-datafiles/hg19_decoy.tar.gz`)
				-x: path of xTEA folder (xTEA can be downloaded with command `git clone https://github.com/parklab/xTEA`)
				Other parameters can be keep unchanged.
			```

	+ Only with 10X data
		```
		sh run_gnrt_pipeline.sh sample_id.txt null 10X_bam_list.txt ./path_work_folder ./rep_lib_annotation/hg19_LINE1_lib_config.txt
		```
		
	+ Working with hybrid data of 10X and Illumina
		```
		sh run_gnrt_pipeline.sh sample_id.txt illumina_bam_list.txt 10X_bam_list.txt ./path_work_folder ./rep_lib_annotation/hg19_LINE1_lib_config.txt
		```
		
	2.2 The previous step will generate a shell script called `run_xTEA_pipeline.sh` under `WFOLDER/sample_id` with content like:
		
		```
		#!/bin/bash

		#SBATCH -n 16
		#SBATCH -t 1-5:00
		#SBATCH --mem=70G
		#SBATCH -p park
		#SBATCH -o hostname_%j.out
		#SBATCH --mail-type=END
		#SBATCH --mail-user=chong.simonchu@gmail.com
		#SBATCH --account=park_contrib
		####
		PREFIX=/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/1000G/hybrid/L1/NA19239/
		############
		############
		SF_FLANK=/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/rep_lib_annotation/LINE/hg38/hg38_FL_L1_flanks_3k.fa
		L1_CNS=/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/rep_lib_annotation/consensus/LINE1.fa
		BAM1=${PREFIX}"10X_phased_possorted_bam.bam"
		L1_COPY_WITH_FLANK=/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/rep_lib_annotation/LINE/hg38/hg38_L1HS_copies_larger_5K_with_flank.fa
		TMP_CNS=${PREFIX}"tmp/cns/"
		BARCODE_BAM=${PREFIX}"10X_barcode_indexed.sorted.bam"
		ANNOTATION=/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/rep_lib_annotation/LINE/hg38/hg38_L1_larger2K_with_all_L1HS.out
		BAM_LIST=${PREFIX}"bam_list.txt"
		TMP=${PREFIX}"tmp/"
		TMP_CLIP=${PREFIX}"tmp/clip/"
		REF=/n/data1/hms/dbmi/park/simon_chu/projects/data/reference/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
		XTEA_PATH=/n/data1/hms/dbmi/park/simon_chu/projects/XTEA/
		############
		############	
		python ${XTEA_PATH}"x_TEA_main.py" -C -i ${BAM_LIST} --lc 5 --rc 5 --cr 3  -r ${L1_COPY_WITH_FLANK}  -a ${ANNOTATION} --ref ${REF} -p ${TMP} -o ${PREFIX}"candidate_list_from_clip.txt"  -n 8
		python ${XTEA_PATH}"x_TEA_main.py"  -D -i ${PREFIX}"candidate_list_from_clip.txt" --nd 6 --ref ${REF} -a ${ANNOTATION} -b ${BAM_LIST} -p ${TMP} -o ${PREFIX}"candidate_list_from_disc.txt" -n 8
		python ${XTEA_PATH}"x_TEA_main.py" -N --cr 3 --nd 5 -b ${BAM_LIST} -p ${TMP_CNS} --fflank ${SF_FLANK} --flklen 3000 -n 8 -i ${PREFIX}"candidate_list_from_disc.txt" -r ${L1_CNS} --ref ${REF} -a ${ANNOTATION} -o ${PREFIX}"candidate_disc_filtered_cns.txt"
		```
	+ To run on the cluster, like slurm system: `sbatch < run_xTEA_pipeline.sh`
	+ To run on a single node: `sh run_xTEA_pipeline.sh`
	
3. **Output**

	3.1 For Illumina, `candidate_disc_filtered_cns.txt` is the final output.