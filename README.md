### cli branch

$${\color{red}This branch of xTea is currently in development. Results are not guaranteed to match the xTea-master branch. Use with caution.}$$

This branch of xTea refactors installation of xTea as well as the command line interface. It currently does not support case-control running or 10X barcoded reads.

It does support Illumina short-reads and can be run in germline or mosaic mode. Additionally only hg38 is supported at this time.


## xTea

xTea (comprehensive transposable element analyzer) is designed to identify TE insertions from paired-end Illumina reads, barcode linked-reads, long reads (PacBio or Nanopore), or hybrid data from different sequencing platforms and takes whole-exome sequencing (WES) or whole-genome sequencing (WGS) data as input. 

![alt text](./xTea_workflow.png)


## Download

```
git clone --single-branch --branch cli  https://github.com/parklab/xTea.git
```
	
Gene annotation files are downloaded from GENCODE. GFF3 files are required.
	+ For GRCh38 (or hg38), gff3 files are downloaded from https://www.gencodegenes.org/human/release_33.html;

## Dependencies

1. bwa (version **0.7.17** or later, which requires the **-o** option), which can be downloaded from https://github.com/lh3/bwa.
2. samtools (version 1.0 or later), which can be downloaded from https://github.com/samtools.
3. Python 3.6+

***Note: bwa and samtools need to be added to the $PATH.***


## Install

To build and install xTea, run:
```
# only if you need to install building packages
# pip install --upgrade build
# pip install --upgrade setuptools

# within the xTea/ directory
python -m build
pip install .
``` 

## Run xTea

Once xTea is installed you can view all the available options by running:

```
xtea -h

usage: xtea [-h] -c CONFIG --input_bams INPUT_BAMS [INPUT_BAMS ...] --sample_name SAMPLE_NAME
            [--repeat_type REPEAT_TYPE [REPEAT_TYPE ...]] [-m MODE] [-g GENOME] [--cr CR] [--clip_cutoff CLIP_CUTOFF] [--nd ND]
            [--tumor TUMOR] [--purity PURITY] [--single SINGLE] [--extend EXTEND] --rep_lib_annot_dir REP_LIB_ANNOT_DIR
            --onnx_model_file ONNX_MODEL_FILE --genome_reference GENOME_REFERENCE --genome_gff3 GENOME_GFF3 -o OUTPUT_DIR
            [--tmp_dir TMP_DIR] [-n CORES] [--resume RESUME] [--save-intermediate-files] [-v]

options:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        TOML config file path
  --input_bams INPUT_BAMS [INPUT_BAMS ...], -i INPUT_BAMS [INPUT_BAMS ...]
                        input bam(s)
  --sample_name SAMPLE_NAME, -s SAMPLE_NAME
                        Sample identifier
  --repeat_type REPEAT_TYPE [REPEAT_TYPE ...]
                        Type of repeats to detect. Options include: L1, ALU, SVA, HERV, Mitochondrial (Default = ALU L1 SVA)
  -m MODE, --mode MODE  Which mode to run. Options: germline, mosaic, case-control, denovo (Default: germline)
  -g GENOME, --genome GENOME
                        Genome version used. Options: hg38 (default), hg19, chm13
  --cr CR               When specified, override default automatic calculation: cutoff of minimum # of clipped parts fall in repeats
  --clip_cutoff CLIP_CUTOFF
                        When specified, override default automatic calculation: cutoff of minimum # of clipped reads (used for both
                        left & right sides)
  --nd ND               When specified, override default automatic calculation: cutoff of minimum # of discordant pair
  --tumor TUMOR         Working on tumor samples
  --purity PURITY       Tumor purity
  --single SINGLE       Call clip positions from single-end reads
  --extend EXTEND       Number of bp to extend gff annotation for annotations
  --rep_lib_annot_dir REP_LIB_ANNOT_DIR
                        Path to rep_lib_annotation/ directory
  --onnx_model_file ONNX_MODEL_FILE
                        Path to .onnx model file
  --genome_reference GENOME_REFERENCE
                        Path to genome fasta file
  --genome_gff3 GENOME_GFF3
                        Path to genome annotation file (.gff3 format)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory
  --tmp_dir TMP_DIR     Temporary directory
  -n CORES              number of cores
  --resume RESUME       resume previous run if available
  --save-intermediate-files
                        Save all intermediate files
  -v                    version

Args that start with '--' can also be set in a config file (params.toml or specified via -c). Uses configparser module to parse an
INI file which allows multi-line values. See https://docs.python.org/3/library/configparser.html for details. This parser includes
support for quoting strings literal as well as python list syntax evaluation.  In general, command-line values override config file
values which override defaults.
```

xTea utilizes a TOML parameter file for specifying run parameters. Every parameter specified in the TOML file can also be specified on the command line. Command line specifications take precedence over the TOML file. An example TOML file can be found at `xTea/params/default_params.toml`. 

			
4. **Output**

	A gVCF file will be generated with insertion calls underneath the specified `output_dir`.
	+ For germline TE insertion calling on short reads, the `orphan transduction` module usually has a higher false positive rate. Users can filter out false positive events with a command such as `grep -v "orphan" out.vcf > new_out.vcf` to retrieve higher confidence events.


5. **Citation and accompanying scripts**
	If you are using xTea for your project, please cite:
	
	```
	Chu, C., Borges-Monroy, R., Viswanadham, V.V. et al. Comprehensive identification of transposable element insertions using multiple sequencing technologies. Nat Commun 12, 3836 (2021). https://doi.org/10.1038/s41467-021-24041-8
	```

	The accompanying scripts for reproduce the results in the paper can be found here: `https://github.com/parklab/xTea_paper`
