# das_tool 
#
# VERSION       latest

FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
RUN DEBIAN_FRONTEND=noninteractive apt-get -qq update
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gcc
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y make 
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y build-essential
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y libtool
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y automake
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y autoconf
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y autopoint
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y bison
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gettext
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y gperf
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y texinfo
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y rsync
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y coreutils
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y libpcre3
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y libpcre3-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y zlib1g-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y wget
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y liburi-escape-xs-perl
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y littler
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y ruby
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y cpanminus
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y ncbi-blast+ 
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y hmmer
RUN DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y git
RUN DEBIAN_FRONTEND=noninteractive apt-get clean

# variables
ENV DAS_TOOL_DIR /DAS_Tool
ENV VERSION 1.1.0
RUN mkdir $DAS_TOOL_DIR
WORKDIR $DAS_TOOL_DIR

# install R packages
ADD install_dependences.R $DAS_TOOL_DIR
RUN Rscript install_dependences.R
RUN rm install_dependences.R

# install from github
RUN wget https://github.com/cmks/DAS_Tool/archive/${VERSION}.tar.gz
RUN tar zxfv ${VERSION}.tar.gz
RUN rm ${VERSION}.tar.gz
RUN ln -s DAS_Tool-${VERSION} current
RUN chmod +x current/DAS_Tool
WORKDIR $DAS_TOOL_DIR/current
RUN unzip db.zip
RUN rm db.zip
RUN R CMD INSTALL ./package/DASTool_1.1.0.tar.gz

# install pullseq
RUN wget https://github.com/bcthomas/pullseq/archive/master.zip
RUN unzip master.zip
WORKDIR $DAS_TOOL_DIR/current/pullseq-master
RUN autoreconf -i
RUN ./configure
RUN make
RUN make install
WORKDIR $DAS_TOOL_DIR/current
RUN rm -rf pullseq-master && rm master.zip

# install prodigal
RUN wget https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz
RUN tar zxfv v2.6.3.tar.gz 
WORKDIR $DAS_TOOL_DIR/current/Prodigal-2.6.3 
RUN make
RUN make install
WORKDIR $DAS_TOOL_DIR/current
RUN rm -rf Prodigal-2.6.3 && rm v2.6.3.tar.gz

# install diamond
RUN wget http://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz
RUN tar zxfv diamond-linux64.tar.gz 
RUN rm diamond-linux64.tar.gz
RUN mv diamond /usr/local/bin/

# add to PATH 
ENV PATH="$PATH":$DAS_TOOL_DIR/current
ADD usearch /usr/local/bin/

# Entry
VOLUME /mnt
WORKDIR /mnt
CMD [ "/bin/bash" ]

