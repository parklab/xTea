
## xTea-mosaic

xTea-mosaic is derived from xTea to identify mosaic TE insertions from bulk high-depth PE short WGS data. The current version is for testing only, not guranteed for stable output.

![alt text](./xTea_workflow.png)


## Download

1. short reads (Illumina; High-depth)

	```
	git clone --single-branch --branch xtea_mosaic https://github.com/parklab/xTea.git
	```
2. pre-processed repeat library used by xTea (this library is used for both short and long reads)  
	
	```
	wget https://github.com/parklab/xTea/raw/master/rep_lib_annotation.tar.gz
	```
	
3. gene annotation files are downloaded from GENCODE. Decompressed gff3 files are required.
	+ For GRCh38 (or hg38), gff3 files are downloaded and decompressed from https://www.gencodegenes.org/human/release_33.html ;
	+ For GRCh37 (or hg19), gff3 files are downloaded and decompressed from https://www.gencodegenes.org/human/release_33lift37.html ;
	+ For CHM13v2, gff3 files are downloaded from https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3;
	+ Or use the latest version

## Dependencies

1. bwa (version **0.7.17** or later, which requires the **-o** option), which can be downloaded from https://github.com/lh3/bwa.
2. samtools (version 1.0 or later), which can be downloaded from https://github.com/samtools.
3. minimap2 (for long reads only), which can be downloaded from https://github.com/lh3/minimap2.
4. wtdbg2 (for long reads only), which can be downloaded from https://github.com/ruanjue/wtdbg2.
5. Python 2.7+/3.6+
	+ For the following packages, only a conda-based installation in shown. You may also install these in other ways, such as pip. 
	+ pysam (https://github.com/pysam-developers/pysam, version 0.12 or later) is required.
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
		`conda install numpy scikit-learn=0.18.1 pandas -y`

	+ DF21 (this is used to replease scikit-learn, which is complained by several users for version incompatible)
		+ Install DF21
		`pip install deep-forest`

6. Note: bwa and samtools need to be added to the $PATH.



## Run xTea
1. **Input**
	+ A sample id file list, e.g. a file named `sample_id.txt` with content as follows (each line represents one unique sample id):
	
		```
		NA12878
		NA12877
		```
	
	+ A file of listed alignments:

		+ An Illumina bam/cram file (sorted and indexed) list, e.g. a file named `illumina_bam_list.txt` with content as follows (two columns separated by a space or tab: sample-id bam-path):

			```
			NA12878 /path/na12878_illumina_1_sorted.bam
			NA12877 /path/na12877_illumina_1_sorted.bam
			```


2. **Run the pipeline from local cluster or machine**
	

	2.1 Generate the running script
			
	+ Run on a cluster or a single node (by default xTea_mosaic assumes the reference genome is **GRCh38** or **hg38**. For `hg19` or `GRCh37`, please use `gnrt_pipeline_local.py`; for `CHM13`, please use `gnrt_pipeline_local_chm13.py`)
		+ Here, the slurm system is used as an example. If using LSF, replace `--slurm` with `--lsf`. For those using clusters other than slurm or LSF, users must adjust the generated shell script header accordingly. Users also must adjust the number of cores (`-n`) and memory (`-m`) accordingly. For very high depth bam files, runtime (denoted by `-t`) may take longer.
		+ **Note that `--xtea` is a required option that points to the *exact folder* containing python scripts.**

		+ Using high-depth Illumina data
			```
			python /home/xTea_mosaic/gnrt_pipeline_local_v38.py -M -U -i sample_id.txt -b illumina_bam_list.txt -x null -p ./path_work_folder/ -o submit_jobs.sh -l /home/rep_lib_annotation/ -r /home/reference/genome.fa -g /home/gene_annotation_file.gff3 --xtea /home/ec2-user/xTea/xtea/ -f 5907 -y 1 --slurm -t 10-12:00 -q short -n 8 -m 96 --nclip 2 --cr 0 --nd 1 --nfclip 1 --nfdisc 1 --blacklist /home/germline_insertions.bed
			```

		+ Parameters:
			
			```
			Required:
			    -M: indicates this is for mosaic mode;
			    -U: indicates this is to use user-specific cutoff, rather than automatic parameter setting.
				-i: samples id list (one sample id per line);
				-b: Illumina bam/cram file list (sorted and indexed — one file per line);
				-x: 10X bam file list (sorted and indexed — one file per line);
				-p: working directory, where the results and temporary files will be saved;
				-l: repeat library directory (directory which contains decompressed files from "rep_lib_annotation.tar.gz");
				-r: reference genome fasta/fa file;
				-y: type of repeats to process (1-L1, 2-Alu, 4-SVA, 8-HERV; sum the number corresponding to the repeat type to process multiple repeats. 
				    For example, to run L1 and SVA only, use `-y 5`. 
				    Each repeat type will be processed separately, however some of the early processing steps are common to multiple repeat types.
				    Thus, when analyzing a large cohort, to improve the efficiency (and save money on the cloud), 
				    it is highly recommended to run the tool on one repeat type first, and subsequently on the rest. 
				    For example, first use '-y 1', and for then use '-y 6' in a second run);
				-f: steps to run. (5907 means run all the steps);
				--xtea: this is the full path of the xTea/xtea folder, where the python scripts reside in;
				-g: gene annotation file in gff3 format;
				-o: generated running scripts under the working folder;
				--nclip: minimum number of clipped reads;
				--cr: minimum number of clipped reads whose mates map to repetitive regions;
				--nd: minimum number of discordant pairs;
				--blacklist: for mosaic calling, this one is the germline insertions from population data (e.g. gnomAD). File should be in bed format. Regions do not want to consider can also be added to this file.

			Optional:
				-n: number of cores (default: 8, should be an integer);
				-m: maximum memory in GB (default: 25, should be an integer; for mosaic mode, large memory is needed);
				-q: cluster partition name;
				-t: job runtime (for mosaic mode, long running time is needed);
				--flklen: flanking region length;
				--lsf: add this option if using an LSF cluster (by default, use of the slurm scheduler is assumed);
				--slurm: runs using the slurm scheduler. Generates a script header fit for this scheduler;


			```
		
	2.2 The previous step will generate a shell script called `run_xTea_pipeline.sh` under `WFOLDER/sample_id/L1(or other types of repeats)`, where `WFOLDER` is specified by `-p` option.
		
	+ To run on the script: `sh run_xTea_pipeline.sh` or users can submit the jobs (where each line corresponds to one job) to a cluster.
	
			
3. **Output**

	A gVCF file will be generated for each sample. Sensitivity is more important in detecting mosaic TE insertions, thus "intermediate" output usually are higher recommended:

	```
	First step output only using clipped reads: candidate_list_from_clip.txt
	Combined with discordant checking: candidate_list_from_disc.txt
	Extra filtering in checking alignment patterns: candidate_disc_filtered_cns.txt
	Added transduction calling: candidate_disc_filtered_cns2.txt
	Added several post-filters: candidate_disc_filtered_cns_post_filtering.txt

	```

4. **Update log**
	+ 06/10/25 Fixed bugs in post-filtering step.  

	+ 10/13/23 Release the version for mosaic insertion calling from high-depth WGS.

	+ 06/11/23 Add `gnrt_pipeline_local_chm13.py` for CHM13_v2.0 reference genome .

	+ 06/09/22 Update the Dockerfile and cwl for germline module (hg38).

	+ 04/20/22 A fatal error was noticed at the genotyping step. The machine learing model was trained with features extracted with a old version of xTea, and this will introduce bias to predict the features extracted with the latest version of xTea. A new model is uploaded for non-conda version.
	
	+ 04/20/22 The scikit-learn version issue is complained by several users. To solve this issue, the new genotype classification model is trained with DF21 (https://github.com/LAMDA-NJU/Deep-Forest). Users need to install with command `pip install deep-forest`. For now, this is only for the non-conda version. I'll update the conda version soon.
