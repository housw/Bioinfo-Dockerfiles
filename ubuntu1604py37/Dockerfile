FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gcc make build-essential libtool zlib1g-dev wget ca-certificates && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install miniconda 
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
/bin/bash ~/miniconda.sh -b -p /opt/conda && \
rm ~/miniconda.sh && \
/opt/conda/bin/conda clean -tipsy && \
ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh
ENV PATH=/opt/conda/bin:$PATH

# create py37 env
#ENV PATH=/opt/conda/bin:$PATH
#RUN conda install conda-build python=3.7 \
#    && conda create -y --name py37 python=3.7 \
#    && conda clean -ya
#ENV CONDA_DEFAULT_ENV=py37
#ENV CONDA_PREFIX=/opt/conda/envs/$CONDA_DEFAULT_ENV
#ENV PATH=$CONDA_PREFIX/bin:$PATH

# install required tools and libraries
#RUN conda install -c bioconda -y numpy && \
#    conda clean -ya

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
