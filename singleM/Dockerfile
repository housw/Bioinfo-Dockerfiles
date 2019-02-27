FROM shengwei/ubuntu1604py27:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y hmmer mafft fasttree git && \
    ln -s /usr/bin/fasttree /usr/bin/FastTree && \
    ln -s /usr/bin/fasttreeMP /usr/bin/FastTreeMP && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install graftm
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install -c bioconda diamond pplacer krona orfm graftm vsearch && \
    conda clean -ya

# install fxtract
WORKDIR /
RUN git clone --recursive https://github.com/ctSkennerton/fxtract.git && \
    cd fxtract && \
    make NO_PCRE=1
ENV PATH=$PATH:/fxtract/

# install smafa
WORKDIR /
RUN wget https://github.com/wwood/smafa/releases/download/v0.1.2/smafa-static-x86_64-unknown-linux-musl-0.1.2.tar.gz && \
    tar zxf smafa-static-x86_64-unknown-linux-musl-0.1.2.tar.gz && \
    rm smafa-static-x86_64-unknown-linux-musl-0.1.2.tar.gz && \
    mv smafa-static-x86_64-unknown-linux-musl-0.1.2 smafa
ENV PATH=$PATH:/smafa/

# install python packages
RUN pip install matplotlib==2.0.2

# install singlem
#RUN /opt/conda/bin/pip install singlem
RUN git clone --recursive https://github.com/wwood/singlem.git && \
    chmod -R 777 singlem && cd singlem && \
    python setup.py install && \
    python setup.py clean --all

# run test 
#RUN cd /singlem/test && python test_pipe.py 

WORKDIR /mnt
CMD [ "/bin/bash" ]

