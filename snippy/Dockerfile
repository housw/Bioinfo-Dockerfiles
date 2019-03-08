FROM shengwei/ubuntu1604brew:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get update && \
    apt-get install -yq --no-install-recommends zlib1g-dev libbz2-dev bioperl git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install snippy dependences
RUN brew install brewsci/bio/snippy && \
    brew cleanup 

# get the latest snippy version
RUN git clone https://github.com/tseemann/snippy.git /snippy && \
    chmod -R 777 /snippy && rm /home/linuxbrew/bin/snippy* && \
    rm /home/linuxbrew/bin/bcftools /home/linuxbrew/bin/bwa \
       /home/linuxbrew/bin/samtools /home/linuxbrew/bin/vcfallelicprimitives \
       /home/linuxbrew/bin/vcfstreamsort /home/linuxbrew/bin/vcfuniq 
ENV PATH=$PATH:/snippy/bin:/snippy/binaries/linux

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
