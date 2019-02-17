FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com


# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gcc build-essential cmake python3 python3-dev python3-pip r-base && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# install ppanggolin
RUN pip3 install setuptools && \
    pip3 install wheel && \
    pip3 install cython && \
    pip3 install ppanggolin

# add scripts
ADD install_dependences.R /usr/bin
RUN Rscript /usr/bin/install_dependences.R


#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

