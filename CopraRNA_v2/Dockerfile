# CopraRNA 
#
# VERSION       2.1.2

FROM shengwei/ubuntu:16.04

MAINTAINER Shengwei Hou, housw2010@gmail.com

ENV DEBIAN_FRONTEND noninteractive

RUN chown root:root /tmp && chmod 1777 /tmp

# install miniconda3
RUN apt-get -qq update && apt-get -qq -y install vim curl bzip2 \ 
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \ 
    && bash /tmp/miniconda.sh -b -u -p /anaconda \
    && rm -rf /tmp/miniconda.sh
    
# update PATH and conda
ENV PATH /anaconda/bin:$PATH

# install CopraRNA2
#WORKDIR /home/
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
#RUN conda install -c bioconda intarna --yes
RUN conda install -c bioconda coprarna --yes
RUN conda clean --all --yes

CMD [ "/bin/bash" ]

