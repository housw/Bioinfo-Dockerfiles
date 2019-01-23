# sickle 
#
# VERSION       1.33

FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential zlib1g-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*    

# Variables
ENV SICKLE_DIR /sickle
ENV VERSION 1.33

# Download and uncompress sickle
RUN mkdir $SICKLE_DIR
WORKDIR $SICKLE_DIR
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends wget && \
    wget https://github.com/najoshi/sickle/archive/v${VERSION}.tar.gz --no-check-certificate && \
    tar -pxvzf v${VERSION}.tar.gz && rm v${VERSION}.tar.gz && \
    apt-get remove -y --auto-remove wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install sickle
WORKDIR $SICKLE_DIR/sickle-${VERSION}
RUN make && chmod +x sickle && cp sickle /bin 

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

