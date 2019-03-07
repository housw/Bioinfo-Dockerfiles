FROM shengwei/ubuntu1604brew:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get update && \
    apt-get install -yq --no-install-recommends zlib1g-dev libbz2-dev bioperl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install snippy
RUN brew install brewsci/bio/snippy && \
    brew cleanup 

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
