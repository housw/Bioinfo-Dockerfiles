FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gcc make build-essential libtool zlib1g-dev wget libopenblas-dev ca-certificates xvfb && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install miniconda, ete3, r351
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && rm ~/miniconda.sh && \
    /opt/conda/bin/conda install --yes -c etetoolkit ete3 ete_toolchain && \
    /opt/conda/bin/conda install --yes -c r r-base==3.5.1 && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
ENV PATH=/opt/conda/bin:$PATH

# check ete3 dependencies
RUN xvfb-run ete3 build check

# install r packages
RUN Rscript -e "install.packages(c('BiocManager', 'argparse'), repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install('ggtree', version = '3.8')" 
RUN Rscript -e "BiocManager::install(c('Biostrings', 'seqinr'))"

# add mview
RUN mkdir /mview && chmod -R 777 /mview
ADD mview-1.64 /mview/
ENV PATH=${PATH}:/mview/bin

# add scripts
RUN mkdir /scripts && chmod -R 777 /scripts
ADD scripts /scripts/
ENV PATH=${PATH}:/scripts

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
