# transdecoder 
#
# VERSION       latest

FROM shengwei/ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update && \
    apt-get install -y gcc make build-essential zlib1g-dev wget liburi-escape-xs-perl littler cpanminus ncbi-blast+ hmmer && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# perl lib
RUN cpanm URI::Escape

# variables
ENV TRANSDECODER_DIR /transdecoder
ENV VERSION v5.3.0
RUN mkdir $TRANSDECODER_DIR
WORKDIR $TRANSDECODER_DIR

# install from github
RUN wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-${VERSION}.tar.gz
RUN tar zxfv TransDecoder-${VERSION}.tar.gz
RUN rm TransDecoder-${VERSION}.tar.gz
RUN ln -s TransDecoder-TransDecoder-${VERSION} current

# add to PATH 
ENV PATH="$PATH":$TRANSDECODER_DIR/current:$TRANSDECODER_DIR/current/util

# Entry
VOLUME /mnt
WORKDIR /mnt
CMD [ "/bin/bash" ]

