FROM ubuntu:16.04
MAINTAINER Chong Simon Chu (chong.simon.chu@gmail.com) (Initially by Soo Lee (duplexa@gmail.com)) 

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
    liblz4-tool \
    python3-pip

# conda and pysam
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -p /miniconda -b
ENV PATH=/miniconda/bin:$PATH
RUN conda update -y conda \
    && rm Miniconda3-py38_4.10.3-Linux-x86_64.sh
RUN conda config --add channels r \
    && conda config --add channels bioconda \
    && conda install -c conda-forge libgcc-ng \
    && conda install -c bioconda samtools \
    && conda install -c bioconda bwa \
    && conda install pysam sortedcontainers numpy pandas scikit-learn -y

#install bamsnap
RUN pip install --no-cache-dir bamsnap

#install deep-forest
RUN pip install deep-forest

# download tools
WORKDIR /usr/local/bin
#COPY downloads.sh .
#RUN . downloads.sh

# set path
#ENV PATH=/usr/local/bin/bwa/:$PATH
#ENV PATH=/usr/local/bin/samtools/:$PATH

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# wrapper
COPY *.py *.sh ./
RUN chmod +x *.py

# copy the trained model for genotyping
COPY genotyping ./genotyping

# default command
CMD ["ls /usr/local/bin"]

#
