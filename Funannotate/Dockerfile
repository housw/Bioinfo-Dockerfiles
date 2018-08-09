# funannotate 
#
# VERSION       latest

FROM ubuntu:16.04
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN chown root:root /tmp && chmod 1777 /tmp
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update \
    && apt-get -qq -y install gcc \
                           make \
                           cmake \
                           build-essential \
                           zlib1g-dev \
                           wget \
                           liburi-escape-xs-perl \
                           curl \
                           bzip2 \
                           mysql-server \
                           libdb-dev \
    && apt-get clean
RUN /etc/init.d/mysql start

# install miniconda2
RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /anaconda \
    && rm -rf /tmp/miniconda.sh
ENV PATH /anaconda/bin:$PATH
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda config --add channels etetoolkit

# python lib
RUN  conda update conda \
     && conda install -y python=2 \
                         numpy \
                         matplotlib \
                         seaborn \
                         pandas \
                         scipy \
                         scikit-learn \
                         biopython \
                         goatools \
                         natsort \
                         psutil \
                         requests \
    && conda install -y -c etetoolkit ete3 ete_toolchain \
    && conda clean --all --yes

# external dependences
RUN conda install -y -c bioconda blast \
                                 rmblast \
                                 goatools \
                                 fisher \
                                 bedtools \
                                 blat \
                                 hmmer \
                                 exonerate \
                                 diamond \
                                 tbl2asn \
                                 hisat2 \
                                 ucsc-pslcdnafilter \
                                 samtools \
                                 raxml \
                                 trimal \
                                 mafft \
                                 kallisto \
                                 bowtie2 \
                                 infernal \
                                 mummer \
                                 minimap2 \
                                 perl-app-cpanminus \
                                 repeatmasker \
                                 repeatmodeler \
    && conda clean --all --yes

# perl lib
RUN /anaconda/bin/cpanm URI::Escape \
                        Bio::Perl \
                        Carp \
                        Clone \
                        Data::Dumper \
                        File::Basename \
                        File::Which \
                        Getopt::Long \
                        Hash::Merge \
                        JSON \
                        LWP::UserAgent \
                        Logger::Simple \
                        POSIX \
                        Parallel::ForkManager \
                        Pod::Usage \
                        Scalar::Util::Numeric \
                        Storable \
                        Text::Soundex \
                        Thread::Queue \
                        Tie::File \
                        URI::Escape \
                        YAML \
                        threads \
                        threads::shared \
                        Text::Soundex \
                        DBD::SQLite \
                        DBD::mysql \
                        DBI \
                        DB_File --force

# GeneMask
ENV GENEMASK_DIR /GeneMask
RUN mkdir $GENEMASK_DIR
WORKDIR $GENEMASK_DIR
COPY gm_et_linux_64.tar.gz gm_key_64.gz ${GENEMASK_DIR}/
RUN zcat gm_key_64.gz > /root/.gm_key
RUN tar zxfv gm_et_linux_64.tar.gz && rm gm_et_linux_64.tar.gz
ENV PATH ${GENEMASK_DIR}/gm_et_linux_64/gmes_petap:$PATH

# SignalP 
ENV SignalP_DIR /SignalP
RUN mkdir $SignalP_DIR
WORKDIR $SignalP_DIR
COPY signalp-4.1f.Linux.tar.gz ${SignalP_DIR}/
RUN tar -zxvf signalp-4.1f.Linux.tar.gz \
    && sed -i 's,/usr/cbs/bio/src/signalp-4.1,${SignalP_DIR}/signalp-4.1,g' signalp-4.1/signalp \
    && sed -i 's,#!/usr/bin/perl,#!/usr/bin/env perl,g' signalp-4.1/signalp
ENV PATH ${SignalP_DIR}/signalp-4.1:$PATH

# bamtools
#RUN apt-get -qq -y install cmake
ENV BamTools_DIR /BamTools
RUN mkdir $BamTools_DIR
WORKDIR $BamTools_DIR
RUN wget https://github.com/pezmaster31/bamtools/archive/v2.5.0.tar.gz \
    && tar -zxvf v2.5.0.tar.gz \
    && rm v2.5.0.tar.gz \
    && mv bamtools-2.5.0 bamtools \
    && cd bamtools \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make \
    && make install \
    && cd /usr/include \
    && ln -f -s ../local/include/bamtools/ \
    && cd /usr/lib/ \
    && ln -f -s /usr/local/lib/bamtools/libbamtools.* .

# augustus
ENV AUGUSTUS_DIR /augustus
RUN mkdir $AUGUSTUS_DIR
WORKDIR $AUGUSTUS_DIR
RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.3.tar.gz \
    && tar -zxvf augustus-3.2.3.tar.gz \
    && rm augustus-3.2.3.tar.gz \
    && mv augustus-3.2.3 augustus \
    && cd augustus \
    && make clean \
    && make
ENV PATH ${AUGUSTUS_DIR}/augustus/bin:$PATH

# EVidence modeler
ENV EVidenceModeler_DIR /EVidenceModeler
RUN mkdir $EVidenceModeler_DIR
WORKDIR $EVidenceModeler_DIR
RUN wget https://github.com/nextgenusfs/EVidenceModeler/archive/0.1.3.tar.gz \
    && tar -zxvf 0.1.3.tar.gz \
    && rm 0.1.3.tar.gz \
    && mv EVidenceModeler-0.1.3 evidencemodeler
ENV PATH ${EVidenceModeler_DIR}/evidencemodeler:$PATH

# Trinity
ENV TRINITY_DIR /Trinity
RUN mkdir $TRINITY_DIR
WORKDIR $TRINITY_DIR
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.5.1.tar.gz \
    && tar -zxvf Trinity-v2.5.1.tar.gz \
    && rm Trinity-v2.5.1.tar.gz \
    && mv trinityrnaseq-Trinity-v2.5.1 Trinity \
    && cd Trinity \
    && make \
    && make plugins
ENV PATH ${TRINITY_DIR}/Trinity:$PATH

# config repeatmasker and repeatmodeler
WORKDIR /anaconda/share/RepeatMasker
COPY RepBaseRepeatMaskerEdition-20170127.tar.gz /anaconda/share/RepeatMasker/
COPY repeatmasker.txt /anaconda/share/RepeatMasker/
RUN tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz \
    && rm RepBaseRepeatMaskerEdition-20170127.tar.gz \
    && cp /anaconda/share/RepeatMasker/util/* /anaconda/bin/ \
    && cd /anaconda/share/RepeatMasker \
    && perl ./configure < repeatmasker.txt
WORKDIR /anaconda/share/RepeatModeler
COPY nseg.tar.gz /anaconda/share/RepeatModeler/ 
COPY repeatmodeler.txt /anaconda/share/RepeatModeler/
RUN tar zxvf nseg.tar.gz \
    && rm nseg.tar.gz \
    && cd nseg \
    && make \
    && mv nmerge nseg /usr/bin \
    && cd ../ && rm -rf nseg \
    && perl ./configure < repeatmodeler.txt

# PASA
ENV PASA_DIR /PASA
RUN mkdir $PASA_DIR
WORKDIR $PASA_DIR
RUN wget https://github.com/PASApipeline/PASApipeline/archive/pasa-v2.2.0.tar.gz \
    && tar -zxvf pasa-v2.2.0.tar.gz \
    && rm pasa-v2.2.0.tar.gz \
    && mv PASApipeline-pasa-v2.2.0 PASApipeline \
    && cd PASApipeline \
    && make clean \
    && make
ENV PATH ${PASA_DIR}/PASApipeline/bin:$PATH

# funannotate
ENV FUNANNOTATE_DIR /Funannotate
ENV FUNANNOTATE_DB /db/Funannotate_db
RUN mkdir -p $FUNANNOTATE_DIR $FUNANNOTATE_DB
WORKDIR $FUNANNOTATE_DIR
COPY funannotate-1.0.0.tar.gz ${FUNANNOTATE_DIR}/
RUN tar -zxvf funannotate-1.0.0.tar.gz \
    && rm funannotate-1.0.0.tar.gz \
    && mv funannotate-1.0.0 funannotate
ENV PATH ${FUNANNOTATE_DIR}/funannotate:$PATH
RUN funannotate setup -d $FUNANNOTATE_DB
ENV FUNANNOTATE_DB $FUNANNOTATE_DB

# variables
ENV PASAHOME ${PASA_DIR}/PASApipeline
ENV TRINITYHOME ${TRINITY_DIR}/Trinity
ENV EVM_HOME ${EVidenceModeler_DIR}/evidencemodeler
ENV AUGUSTUS_CONFIG_PATH ${AUGUSTUS_DIR}/augustus/config
ENV GENEMARK_PATH ${GENEMASK_DIR}/gm_et_linux_64/gmes_petap
ENV BAMTOOLS_PATH /usr/local/bin

# Entry
VOLUME /mnt
WORKDIR /mnt
CMD [ "/bin/bash" ]

