## xTea repeat library preparation (Version for testing only)

xTea is designed to identify TE insertions from paired-end Illumina reads, barcode linked-reads, long reads (PacBio or Nanopore), or hybrid data from different sequencing platforms and takes whole-exome sequencing (WES) or whole-genome sequencing (WGS) data as input. 

![alt text](./xTea_workflow.png)


## Dependencies

1. bwa (version **0.7.17** or later, which requires the **-o** option), which can be downloaded from https://github.com/lh3/bwa.
2. samtools (version 1.0 or later), which can be downloaded from https://github.com/samtools.
3. Python 2.7+/3.6+
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

4. Note: bwa and samtools need to be added to the $PATH.



## Run xTea library preparation (Note: You need to prepare for each repeat type separately)
1. **Input**
	+ A TE consensus file in `fasta` format.
	
		```
		path-to-rep-lib-folder/TE-type.consensus.fa
		```
	
		Make sure it is indexed with command `bwa index path-to-rep-lib-folder/TE-type.consensus.fa`

	+ A RepeatMasker output file for the working on TE type:
		```
		path-to-rep-lib-folder/TE-type_rmsk.out
		```
	
	+ A fasta file with full length copies within the reference genome (or copies longer than the specific cutoff):
		To generate this file:

		```
		xtea -P -K -p ./ -r path-of-reference-genome.fa -a path-to-rep-lib-folder/full-length-TE-type_rmsk.out -o path-output-folder/TE_copies_with_flank.fa -e 100 
		```
		Note, if there are transduction events reported for this type to TE before, you may consider set a larger `-e`, like 3000. Also, you only need to keep those full length copies in the RepeatMasker annotation file `path-to-rep-lib-folder/full-length-TE-type_rmsk.out`.


2. **Run the pipeline from local cluster or machine**

	2.1 Adjust the wrapper code (Temorary step for now)
		Adjust the following paths in `gnrt_pipeline_local_v38.py` to the ones you prepared:

		```
		sf_anno = "ANNOTATION " + sf_folder_rep + "MSTA/hg38/hg38_MSTA.out\n"
        sf_anno1 = "ANNOTATION1 " + sf_folder_rep + "MSTA/hg38/hg38_MSTA.out\n"
        sf_copy_with_flank = "L1_COPY_WITH_FLANK " + sf_folder_rep + "MSTA/hg38/hg38_MSTA_copies_with_flank.fa\n"
        sf_cns = "L1_CNS " + sf_folder_rep + "consensus/MSTA.fa\n"

		```

	2.2 Generate the running script (if it is install-free, then use the full path of the downloaded `bin/xtea` instead.)
			
	+ Run on a cluster or a single node (Note `-y` should be 32)
		+ Here, the slurm system is used as an example. If using LSF, replace `--slurm` with `--lsf`. For those using clusters other than slurm or LSF, users must adjust the generated shell script header accordingly. Users also must adjust the number of cores (`-n`) and memory (`-m`) accordingly. In general, each core will require 2-3G memory to run. For very high depth bam files, runtime (denoted by `-t`) may take longer.
		+ **Note that `--xtea` is a required option that points to the *exact folder* containing python scripts.**

		+ Using only Illumina data
			```
			xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p ./path_work_folder/ -o submit_jobs.sh -l /home/rep_lib_annotation/ -r /home/reference/genome.fa -g /home/gene_annotation_file.gff3 --xtea /home/ec2-user/xTea/xtea/ -f 5907 -y 32 --slurm -t 0-12:00 -q short -n 8 -m 25
			```

		+ Using case-ctrl mode
			```
			xtea --case_ctrl --tumor -i sample_id.txt -b case_ctrl_bam_list.txt -p ./path_work_folder/ -o submit_jobs.sh -l /home/ec2-user/rep_lib_annotation/ -r /home/ec2-user/reference/genome.fa -g /home/gene_annotation_file.gff3 --xtea /home/ec2-user/xTea/xtea/ -y 32 -f 5907 --slurm -t 0-12:00 -q short -n 8 -m 25
			```

		

