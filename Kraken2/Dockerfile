FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# duskmasker from blast+ is required 
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential libtool libssl-dev perl zlib1g-dev cpanminus ncbi-blast+ ca-certificates git rsync wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# perl lib
RUN cpanm LWP::Simple

# kraken
WORKDIR /
RUN git clone https://github.com/DerrickWood/kraken2.git && \
    cd kraken2 && ./install_kraken2.sh /usr/bin && chmod -R 777 /kraken2 

ENV PATH=${PATH}:/kraken2/scripts

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]


