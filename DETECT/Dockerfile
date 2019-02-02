FROM shengwei/ubuntu1604py27
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# add conda channels and install required tools and libraries
RUN conda config --add channels bioconda && \
    conda install -c bioconda -y emboss blast && \
    conda clean -ya 

# install python packages
RUN pip install biopython

# install DETECT v2
ENV DETECT_DIR /DETECT
RUN mkdir -p $DETECT_DIR
WORKDIR $DETECT_DIR
RUN wget http://compsysbio.org/projects/DETECTv2/DETECTv2.tar.gz && \
    tar zxfv DETECTv2.tar.gz && \
    rm DETECTv2.tar.gz && chmod -R 777 $DETECT_DIR
ENV PATH $PATH:${DETECT_DIR}/DETECTv2  

#test detect v2
#WORKDIR ${DETECT_DIR}/DETECTv2
#RUN detect.py sample.fasta --output_file sample.out --top_predictions_file sample.top --num_threads 8 --fbeta_file sample.f1.out --beta 1

# Entry
WORKDIR $DETECT_DIR/DETECTv2
CMD [ "/bin/bash" ]

