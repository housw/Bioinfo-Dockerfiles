FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gcc make build-essential zlib1g-dev wget git openjdk-8-jdk && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install FastViromeExplorer
ENV FastViromeExplorer_DIR /FastViromeExplorer
RUN mkdir $FastViromeExplorer_DIR 
WORKDIR $FastViromeExplorer_DIR
RUN git clone https://code.vt.edu/saima5/FastViromeExplorer.git
WORKDIR ${FastViromeExplorer_DIR}/FastViromeExplorer
RUN cp tools-linux/kallisto /usr/local/bin/ && \
    cp tools-linux/samtools /usr/local/bin/ && \
    javac -d bin src/*.java && \
    chmod -R 777 $FastViromeExplorer_DIR
 
# add to PATH 
ENV PATH="$PATH":${FastViromeExplorer_DIR}/FastViromeExplorer/bin:${FastViromeExplorer_DIR}/FastViromeExplorer/utility-scripts

# Entry
VOLUME /mnt
WORKDIR /mnt
CMD [ "/bin/bash" ]

