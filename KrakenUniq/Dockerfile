FROM ubuntu:14.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# duskmasker from blast+ is required 
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential libtool libssl-dev perl zlib1g-dev cpanminus ncbi-blast+ && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# perl lib
RUN cpanm LWP::Simple

# krakenuniq
WORKDIR /tmp
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends git wget ca-certificates && \
    git clone https://github.com/fbreitwieser/krakenuniq && \
    cd krakenuniq && ./install_krakenuniq.sh -j /usr/bin/ && \
    rm -rf /tmp/krakenuniq && \
    apt-get remove -y --auto-remove git wget ca-certificates && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]


