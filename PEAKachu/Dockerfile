# Peakachu
#
# VERSION       0.1.0

FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

ENV DEBIAN_FRONTEND noninteractive


RUN chown root:root /tmp && chmod 1777 /tmp

# install miniconda3
RUN apt-get -qq update && apt-get -qq -y install vim curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh \
    && rm -rf /tmp/miniconda.sh \
    && conda update conda

# install peakachu
WORKDIR /home/
COPY ./peakachu2.yaml ./
RUN conda env create --yes -n peakachu -f peakachu2.yaml
RUN conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH

CMD [ "/bin/bash" ]

