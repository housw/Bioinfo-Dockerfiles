# sickle 
#
# VERSION       1.33

FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gcc
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y make 
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y build-essential
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y zlib1g-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y wget

# Variables
ENV SICKLE_DIR /sickle
ENV VERSION 1.33

# Download and uncompress sickle
RUN mkdir $SICKLE_DIR
WORKDIR $SICKLE_DIR
RUN wget https://github.com/najoshi/sickle/archive/v${VERSION}.tar.gz --no-check-certificate
RUN tar -pxvzf v${VERSION}.tar.gz
RUN rm v${VERSION}.tar.gz

# Install sickle
WORKDIR $SICKLE_DIR/sickle-${VERSION}
RUN make
RUN chmod +x sickle
RUN cp sickle /bin 

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

