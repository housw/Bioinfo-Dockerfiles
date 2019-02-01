FROM shengwei/ubuntu1604py27
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y default-jre && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install blast
RUN conda install -c bioconda blast-legacy && \
    conda clean -ya 

# install priam database 
ENV PRIAM_DIR /PRIAM
RUN mkdir -p $PRIAM_DIR && chmod -R 777 $PRIAM_DIR
WORKDIR $PRIAM_DIR
COPY PRIAM_search.jar Distribution/ $PRIAM_DIR/

# test priam search
#ENV threads 5
#ENV session_name priam_job
#ADD sample.fasta .
#RUN java -jar -Xmx4096m PRIAM_search.jar -np $threads -n $job_name -i sample.fasta -p PRIAM_JAN18 -od priam_output_dir --min_proba 0.5 --min_proportion 70 --check_catalytic

# Entry
WORKDIR $PRIAM_DIR
CMD [ "/bin/bash" ]

