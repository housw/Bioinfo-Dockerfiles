FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gcc cpp make build-essential libtool libncurses5-dev libbz2-dev liblzma-dev wget zlib1g-dev zip unzip bzip2 \
                    hmmer libtbb-dev python-pip python-setuptools && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install python packages
RUN python -m pip install --upgrade pip && \
    pip install wheel && \
    pip install nose numpy pandas pysam biopython

# install vsearch
WORKDIR /
RUN wget https://github.com/torognes/vsearch/releases/download/v2.11.0/vsearch-2.11.0-linux-x86_64.tar.gz && \
    tar xzf vsearch-2.11.0-linux-x86_64.tar.gz && rm vsearch-2.11.0-linux-x86_64.tar.gz && \
    mv vsearch-2.11.0-linux-x86_64 vsearch && chmod -R 777 vsearch

# install MIDAS
# don't run setup.py from /, distutils will do a listdir('.') and hang there
# https://github.com/moby/moby/issues/10025
RUN wget https://github.com/snayfach/MIDAS/archive/master.zip && \
    unzip master.zip && rm master.zip && \
    mv MIDAS-master MIDAS && chmod -R 777 MIDAS && \
    cd MIDAS && python setup.py install && cd / 

# moutnt midas_db_v1.2 to /midas_db_v1.2 while run this image
ENV MIDAS_DB=/midas_db_v1.2
RUN mkdir ${MIDAS_DB} && chmod -R 777 ${MIDAS_DB}

# install samtools v1.4
RUN wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 && \
    bunzip2 samtools-1.4.tar.bz2 && tar xvf samtools-1.4.tar && rm samtools-1.4.tar && \
    cd samtools-1.4 && ./configure && make && make install && \
    cd /MIDAS/bin/Linux && rm samtools && ln -s /usr/local/bin/samtools . 

# update PATH 
ENV PATH=$PATH:/vsearch/bin:/MIDAS/bin/Linux:/MIDAS/scripts
ENV PYTHONPATH=$PYTHONPATH:/MIDAS 

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

