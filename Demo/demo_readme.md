## Note: Here we use `NA12878` downloaded from Illumina (Platinum Genomes) as the demo input bam. Please replace it with your own bam if you already have a processed one.



## 1. Download and preprocess the bam file
### 1.1 Download
From `https://github.com/Illumina/PlatinumGenomes` (`https://www.ebi.ac.uk/ena/data/view/PRJEB3381`), specifically with wget:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194147/NA12878_S1.bam
```
### 1.2 Preprocess (make sure the bam is sorted and indexed)
```
samtools index NA12878_S1.bam
```

## 3. Prepare the running command
Change the abosolute path in file `run_gnrt_pipeline.sh`, and then run  
`sh run_gnrt_pipeline.sh`  
This will generate a shell script named `submit_jobs.sh`


## 4. Submit the job to cluster or run locally
`sh submit_jobs.sh `

## 5. Demo output
The demo ouput are 3 gVCF files (for L1, Alu and SVA respecitively) under the `output` folder.