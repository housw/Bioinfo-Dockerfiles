FROM ubuntu:18.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gcc make cmake build-essential libtool zlib1g-dev libbz2-dev ca-certificates git libseqan2-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install bam2tcc
ENV bam2tcc_DIR="/bam2tcc"
WORKDIR /
RUN git clone https://github.com/pachterlab/bam2tcc.git && \
    cd bam2tcc &&  mkdir build && cd build && \
    cmake .. -DCMAKE_PREFIX_PATH=/usr/share/cmake -DCMAKE_INCLUDE_PATH=/usr/include && \
    make && \
    chmod -R 777 ${bam2tcc_DIR}
ENV PATH=${PATH}:${bam2tcc_DIR}/build/src

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

