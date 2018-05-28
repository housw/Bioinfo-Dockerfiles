# InterProScan 
#
# VERSION       5.28-67.0

FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y default-jre
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y python2.7 
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y wget

# Variables
ENV IPR_DIR /interproscan
ENV VERSION 5.28-67.0
ENV BIT 64
ENV PANTHER_VER 12.0

# Download and uncompress InterProScan
RUN mkdir $IPR_DIR
WORKDIR $IPR_DIR
RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.28-67.0/interproscan-${VERSION}-${BIT}-bit.tar.gz
RUN tar -pxvzf interproscan-${VERSION}-${BIT}-bit.tar.gz
RUN rm interproscan-${VERSION}-${BIT}-bit.tar.gz

# Install Panther Models
WORKDIR $IPR_DIR/interproscan-${VERSION}/data/
RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-${PANTHER_VER}.tar.gz
RUN tar -pxvzf panther-data-${PANTHER_VER}.tar.gz
RUN rm panther-data-${PANTHER_VER}.tar.gz

# Reinstall rpsblast (https://github.com/ebi-pf-team/interproscan/issues/44)
# https://github.com/ebi-pf-team/interproscan/issues/35
# Or replace blast/ncbi-blast-2.6.0+/rpsblast by /blast/2.2.24/bin/rpsblast
WORKDIR $IPR_DIR/interproscan-${VERSION}/bin/blast
RUN cp ./2.2.24/bin/rpsblast ./ncbi-blast-2.6.0+/

# update PATH
WORKDIR /$IPR_DIR/interproscan-${VERSION}/
RUN ln -s interproscan.sh interproscan
RUN chmod +x interproscan
ENV PATH $IPR_DIR/interproscan-${VERSION}:$PATH

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

