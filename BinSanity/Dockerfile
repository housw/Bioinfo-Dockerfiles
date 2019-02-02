FROM shengwei/checkm:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com


# install dependences 
RUN conda install --yes -c bioconda subread bedtools bowtie2 samtools && \
    conda clean -ya

# install python packages
RUN pip install scikit-learn pandas biopython BinSanity

# update PATH
ENV PATH $PATH:/opt/conda/bin

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

