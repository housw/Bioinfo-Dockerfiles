FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 LC_CTYPE=en_US.UTF-8 LANGUAGE=en_US.UTF-8
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y build-essential curl file git locales ca-certificates && \
    locale-gen en_US en_US.UTF-8 && \
    dpkg-reconfigure locales && \
    update-locale LANG=en_US.UTF-8 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install linuxbrew
RUN mkdir /home/linuxbrew/
RUN git clone https://github.com/Homebrew/brew /home/linuxbrew/ && \
    chmod -R 777 /home/  
ENV PATH=${PATH}:/home/linuxbrew/bin
RUN brew upgrade

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

