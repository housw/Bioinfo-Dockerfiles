# transdecoder 
#
# VERSION       latest

FROM shengwei/ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gcc
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y make 
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y build-essential
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y zlib1g-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y wget
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y liburi-escape-xs-perl
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y cpanminus
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y ncbi-blast+ 
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y hmmer
RUN DEBIAN_FRONTEND=noninteractive apt-get clean

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

