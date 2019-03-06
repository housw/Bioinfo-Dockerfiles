FROM shengwei/ubuntu1604py27:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y git cmake libhdf5-serial-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install karp
ENV Karp_DIR="/Karp"
WORKDIR /
RUN git clone https://github.com/mreppell/Karp.git && \
    cd Karp &&  mkdir build && cd build && cmake .. && make && \
    chmod -R 777 ${Karp_DIR}
ENV PATH=${PATH}:${Karp_DIR}/build/src

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

