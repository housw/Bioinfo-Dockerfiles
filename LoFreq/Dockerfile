FROM shengwei/ubuntu1604py27:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y libncurses5-dev libncursesw5-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /

# install samtools v1.1 
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.1/samtools-1.1.tar.bz2 && \
    bunzip2 samtools-1.1.tar.bz2 && tar xf samtools-1.1.tar && rm samtools-1.1.tar && \
    cd samtools-1.1 && make && make install && \
    cd / && chmod -R 777 samtools-1.1 

# install lofreq 
RUN wget https://sourceforge.net/projects/lofreq/files/lofreq_star-2.1.2_linux-x86-64.tgz && \
    tar zxf lofreq_star-2.1.2_linux-x86-64.tgz && rm lofreq_star-2.1.2_linux-x86-64.tgz && \
    chmod -R 777 /lofreq_star-2.1.2
ENV PATH=${PATH}:/lofreq_star-2.1.2/bin

# install GATK
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y zip unzip default-jre && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.1.0/gatk-4.1.1.0.zip && \
    unzip gatk-4.1.1.0.zip && rm gatk-4.1.1.0.zip && chmod -R 777 gatk-4.1.1.0 
ENV PATH=${PATH}:/gatk-4.1.1.0

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

