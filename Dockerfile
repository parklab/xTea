FROM ubuntu:16.04
MAINTAINER Soo Lee (duplexa@gmail.com)

# 1. general updates & installing necessary Linux components
RUN apt-get update -y && apt-get install -y \
    bzip2 \
    gcc \
    git \
    less \
    libncurses-dev \
    make \
    time \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    liblz4-tool

# conda and pysam
RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && bash Miniconda2-latest-Linux-x86_64.sh -p /miniconda2 -b
ENV PATH=/miniconda2/bin:$PATH
RUN conda update -y conda \
    && rm Miniconda2-latest-Linux-x86_64.sh
RUN conda config --add channels r \
    && conda config --add channels bioconda \
    && conda install pysam==0.14.1 -y \
    && conda install sortedcontainers -y

# download tools
WORKDIR /usr/local/bin
COPY downloads.sh .
RUN . downloads.sh

# set path
ENV PATH=/usr/local/bin/bwa/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# wrapper
COPY *.pyc *.sh ./
RUN chmod +x *.pyc

# default command
CMD ["ls /usr/local/bin"]

