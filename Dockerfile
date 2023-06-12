#######################################################################
#     Basic image
#######################################################################
FROM ubuntu:16.04
MAINTAINER Chong Simon Chu (chong.simon.chu@gmail.com), Michele Berselli (berselli.michele@gmail.com)

#######################################################################
#     General updates & installing necessary Linux components
#######################################################################
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

#######################################################################
#     Setting working env
#######################################################################
WORKDIR /usr/local/bin

#######################################################################
#     Software
#######################################################################
## conda and pysam
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -p /miniconda -b
ENV PATH=/miniconda/bin:$PATH
RUN conda config --add channels r \
    && conda config --add channels bioconda \
    && conda install -c conda-forge libgcc-ng \
    && conda install -c bioconda samtools==1.6 \
    && conda install -c bioconda bwa==0.7.17 \
    && conda install pysam==0.17.0 sortedcontainers==2.4.0 numpy==1.22.3 pandas==1.4.2 scikit-learn==1.0.2 -y

## deep-forest
RUN pip install deep-forest==0.1.5

#######################################################################
#     Scripts
#######################################################################
## xTea
RUN git clone https://github.com/parklab/xTea.git && \
    cd xTea && \
    git checkout v0.1.9 && \
    cd .. && \
    cp -r xTea/xtea/* .
RUN rm -rf xTea
RUN chmod +x *.py

#######################################################################
#     Setting env variables
#######################################################################
## Supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

CMD ["ls /usr/local/bin"]
