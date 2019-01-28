FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com


# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gfortran llvm libncurses5-dev build-essential cmake python python-dev python-pip samtools bowtie2 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# install SPARSE
RUN pip install setuptools && \
    pip install wheel && \
    pip install  meta-sparse


#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

