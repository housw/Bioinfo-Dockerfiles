FROM shengwei/checkm:1.0.13
MAINTAINER Shengwei Hou, housw2010@gmail.com

# install dependences 
RUN conda install --yes -c bioconda diamond krona && \
    conda clean -ya

# install python packages
RUN pip install refinem

# update PATH
ENV PATH $PATH:/opt/conda/bin

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

