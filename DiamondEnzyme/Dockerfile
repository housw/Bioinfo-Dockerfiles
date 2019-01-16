FROM shengwei/ubuntu1604py27
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y gzip zip unzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# add conda channels
RUN conda config --add channels r && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda

# install required tools and libraries
RUN conda install -c bioconda -y diamond && \
    conda clean -ya 

# install diamond database
ENV EC_Annotation_DIR /EC_annotation
ENV Statistics_Analysis_DIR /Statistics_analysis
COPY Metabolic-pathways-in-microbiomes /
WORKDIR $EC_Annotation_DIR
RUN python Download_EC_Annotated_SwissProt.py && \
    diamond makedb -p 6 --in uniprot_sprot_annotated.fasta -d uniprot_sprot_annotated && \
    rm uniprot_sprot.fasta && chmod -R 777 $EC_Annotation_DIR

# test diamond
#RUN diamond blastp --query /DETECT/DETECTv2/sample.fasta --db uniprot_sprot_annotated --outfmt 6 qseqid sseqid qstart qend sstart send evalue bitscore qcovhsp slen pident --out diamond_output_file.txt --evalue 1e-10 --max-target-seqs 1 && \
#    python extract_diamond.py diamond_output_file.txt diamond_inferred_EC_file.tsv SwissProt_EC_Mapping.tsv

# Entry
WORKDIR $EC_Annotation_DIR
CMD [ "/bin/bash" ]

