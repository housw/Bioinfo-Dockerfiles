FROM shengwei/ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive

# install required tools
RUN /root/miniconda3/bin/conda install --yes python=3.5 && \
    /root/miniconda3/bin/conda install --yes -c bioconda seqtk spades amos cd-hit && \
    ln -s /root/miniconda3/bin/show-coords /usr/local/bin/show-coords && \
    ln -s /root/miniconda3/bin/blat /usr/local/bin/blat
COPY opt /opt/

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

