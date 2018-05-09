# centrifuge 
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
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y git

# variables
ENV CENTRIFUGE_DIR /centrifuge
RUN mkdir $CENTRIFUGE_DIR
WORKDIR $CENTRIFUGE_DIR

# install from github
RUN git clone https://github.com/infphilo/centrifuge
WORKDIR $CENTRIFUGE_DIR/centrifuge
# fix gcc flags
# https://github.com/infphilo/centrifuge/issues/108 
RUN sed -i 's/^RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3$/RELEASE_FLAGS  = -O3 $(BITS_FLAG) $(SSE_FLAG) -funroll-loops -g3 -std=c++14 -fpermissive/' Makefile
RUN make
RUN make install prefix=/usr/local

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

