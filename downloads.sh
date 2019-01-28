#!/bin/sh


## SOFTWARE: bwa
## VERSION: 0.7.17
## TYPE: aligner
## SOURCE_URL: https://github.com/lh3/bwa
wget https://github.com/lh3/bwa/archive/v0.7.17.tar.gz
tar -xzf v0.7.17.tar.gz
cd bwa-0.7.17
make
cd ..
ln -s bwa-0.7.17 bwa


## SOFTWARE: samtools
## VERSION: 1.2
## TYPE: file format converter
## SOURCE_URL: https://github.com/samtools/samtools
wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
tar -xjf samtools-1.2.tar.bz2
cd samtools-1.2
make
cd ..
ln -s samtools-1.2 samtools
