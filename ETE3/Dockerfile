FROM shengwei/ubuntu1604py27
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y xvfb libopenblas-base && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install ete3
RUN conda install --yes -c etetoolkit ete3 ete_toolchain && \
    conda clean -ya 

# check installation
RUN xvfb-run ete3 build check

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

