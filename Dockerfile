FROM debian:bullseye

LABEL org.opencontainers.image.authors="Chong Simon Chu chong.simon.chu@gmail.com"
LABEL org.opencontainers.image.contributors="Soo Lee duplexa@gmail.com\
  Alexander Solovyov alexander.solovyov@gmail.com"

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    bzip2 gcc git less libncurses-dev make time unzip vim wget zlib1g-dev \
    liblz4-tool python3 python3-pysam samtools bwa python3-pip \
    python3-sortedcontainers python3-dev python3-setuptools \
    python-is-python3

RUN pip install --no-cache-dir bamsnap numpy scipy pandas scikit-learn deep-forest

# clone the code
RUN mkdir -p /opt/xtea/annotation
WORKDIR /opt/xtea
ADD . /opt/xtea
RUN rm rep_lib_annotation.tar.gz && \
  wget https://github.com/parklab/xTea/raw/master/rep_lib_annotation.tar.gz && \
  tar -C /opt/xtea/annotation -xf rep_lib_annotation.tar.gz && \
  rm rep_lib_annotation.tar.gz

# supporting UTF-8
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV TERM=xterm
ENV PATH="/opt/xtea/bin:${PATH}"

# default command
CMD ["/bin/bash"]
