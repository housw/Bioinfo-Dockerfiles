FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gcc make build-essential libtool zlib1g-dev wget ca-certificates && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# install rscape
ENV Rscape_DIR /rscape
WORKDIR /
RUN wget http://eddylab.org/software/rscape/rscape.tar.gz && \
    tar zxfv rscape.tar.gz && mv rscape_v* ${Rscape_DIR} && \
    cd ${Rscape_DIR} && ./configure && make && make install && \
    cd / && rm rscape.tar.gz 

ENV PATH=${PATH}:${Rscape_DIR}/bin

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
