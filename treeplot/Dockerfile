FROM shengwei/ubuntu1604py37r351
MAINTAINER Shengwei Hou, housw2010@gmail.com

# install argparse and ggtree
RUN conda install --yes -c bioconda r-argparse && \
    conda clean -ya
RUN Rscript -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "BiocManager::install('ggtree', version = '3.8')" 
RUN Rscript -e "BiocManager::install('Biostrings')"

RUN mkdir /scripts && chmod -R 777 /scripts
ADD scripts /scripts/
ENV PATH=${PATH}:/scripts

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
