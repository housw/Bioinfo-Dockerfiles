FROM ubuntu:12.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make cmake build-essential libtool zlib1g-dev emboss wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# blast
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz && \
    tar zxf blast-2.2.26-x64-linux.tar.gz && rm blast-2.2.26-x64-linux.tar.gz && chmod -R 777 blast-2.2.26
ENV PATH=${PATH}:/blast-2.2.26/bin

# emboss
RUN apt-get -qq update && \
    apt-get install -y emboss-lib && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    wget ftp://emboss.open-bio.org/pub/EMBOSS/old/6.3.1/EMBOSS-6.3.1.tar.gz && \
    tar zxf EMBOSS-6.3.1.tar.gz && rm EMBOSS-6.3.1.tar.gz && cd EMBOSS-6.3.1 && \ 
    sh ./configure --without-x --prefix=/bin && make && make install && \
    cd / && chmod -R 777  EMBOSS-6.3.1 && \
    cp -r /usr/lib/emboss/ /bin/lib
ENV PATH=${PATH}:/bin/bin:/EMBOSS-6.3.1/emboss

# primer3
RUN apt-get -qq update && \
    apt-get install -y primer3 hmmer ncbi-blast+ littler && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# clustalw
RUN wget ftp://ftp.ebi.ac.uk/pub/software/clustalw2/2.0.12/clustalw-2.0.12.tar.gz && \
    tar zxf clustalw-2.0.12.tar.gz && rm clustalw-2.0.12.tar.gz && cd clustalw-2.0.12 && \
    sh ./configure && make && make install && cd / && \
    chmod -R 777 /clustalw-2.0.12 && ln -s /clustalw-2.0.12/src/clustalw2 /usr/bin/clustalw

# growthpred v1.08, custom version based on growthpred-v1.07
# 1) change fasta name length via "#define FASTA_NAMLEN  256" in cov_fasta.h and cub_fasta.h
# 2) add an option to parse in input dir "--dirinput" in growthpred-v1.08.py 
# 3) remove dirR, use system R
ADD growthpred-v1.08 /growthpred
ENV PATH=${PATH}:/growthpred

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]


