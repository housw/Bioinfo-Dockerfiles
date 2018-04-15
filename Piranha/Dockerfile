# Piranha 1.2.1
#
# VERSION       0.1

FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

RUN chown root:root /tmp && chmod 1777 /tmp

# install basic packages
RUN apt-get -qq update \ 
    && apt-get -qq -y install vim \
    && apt-get -qq -y install curl \
    && apt-get -qq -y install bzip2 \
    && apt-get -qq -y install zlib1g zlib1g-dev \
    && apt-get -qq -y install wget \ 
    && apt-get -qq -y install build-essential \
    && apt-get -qq -y install pkg-config \
    && apt-get -qq -y install cmake \
    && apt-get -qq -y install cmake-data \
    && apt-get install bedtools

# Install R and packages
#RUN apt-get -qq -y install r-base-core
#RUN Rscript -e 'install.packages("evir", dependencies=TRUE, repos="http://cran.us.r-project.org")'
#RUN Rscript -e 'install.packages("seqinr", dependencies=TRUE, repos="http://cran.us.r-project.org")'

# Install Python and packages
RUN apt-get -qq -y install python2.7

# Install bioperl and cpanminus
#RUN apt-get -qq -y install bioperl cpanminus
#RUN wget http://search.cpan.org/CPAN/authors/id/F/FA/FANGLY/Statistics-R-0.26.tar.gz
#RUN cpanm Statistics-R-0.26.tar.gz
#RUN rm Statistics-R-0.26.tar.gz

# append /usr/local/lib to LD_IBRARY_PATH 
ENV LD_LIBRARY_PATH /usr/local/lib:/usr/lib:${LD_LIBRARY_PATH}
ENV PATH /usr/local/bin:/usr/bin:${PATH}

# Setup
RUN mkdir piranha
WORKDIR piranha

# Install gsl
RUN wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.14.tar.gz
RUN tar zxf gsl-1.14.tar.gz
RUN cd gsl-1.14/ \
    && ./configure \
    && make \
    && make install \
    && cd ../

# Install bamtools
RUN wget https://github.com/pezmaster31/bamtools/archive/v2.4.2.tar.gz
RUN tar zxf v2.4.2.tar.gz
RUN cd bamtools-2.4.2 \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=/usr/local/bin .. \
    && make install \
    && cd ../../

# Install piranha
RUN wget http://smithlabresearch.org/downloads/piranha-1.2.1.tar.gz \
    && tar zxf piranha-1.2.1.tar.gz \
    && cd piranha-1.2.1 \
    && ./configure \
    && make all \
    && make install \
    && cp ./src/progs/Piranha ./src/progs/Simulate /usr/local/bin/

ENTRYPOINT ["/bin/bash"]

# Cleanup
RUN rm -rf gsl-1.14.tar.gz v2.4.2.tar.gz \
           piranha-1.2.1.tar.gz bamtools-2.4.2 gsl-1.14
    
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

