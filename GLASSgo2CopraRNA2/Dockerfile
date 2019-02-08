FROM shengwei/ubuntu1604py37:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y ncbi-blast+ clustalo bioperl software-properties-common && \
    add-apt-repository ppa:j-4/vienna-rna && apt-get -qq update && apt-get install -y vienna-rna && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install python package
RUN pip install numpy biopython ete3

# add scripts
ADD scripts /GLASSgo/
RUN chmod -R 777 /GLASSgo
ENV PATH $PATH:/GLASSgo/

# update taxa.sqlite
RUN /GLASSgo/createJSON-1_0.py -u true


WORKDIR /mnt
CMD [ "/bin/bash" ]

