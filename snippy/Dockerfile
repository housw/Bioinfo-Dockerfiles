FROM shengwei/ubuntu1604py37
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get update && \
    apt-get install -yq --no-install-recommends zlib1g-dev libbz2-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install snippy
RUN conda install -c conda-forge -c bioconda snippy && \
    conda clean -ya

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
