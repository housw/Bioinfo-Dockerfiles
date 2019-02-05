FROM shengwei/ubuntu1604py37
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get update && \
    apt-get install -yq --no-install-recommends bc curl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install GToTree
WORKDIR /
RUN curl -L https://github.com/AstrobioMike/GToTree/archive/v1.1.3.tar.gz -o GToTree-1.1.3.tar.gz && \
    tar zxf GToTree-1.1.3.tar.gz && rm GToTree-1.1.3.tar.gz && \
    chmod -R 777 GToTree-1.1.3 && \
    cd /GToTree-1.1.3 && /bin/bash conda-setup.sh && \
    conda clean -ya 

# env
ENV CONDA_DEFAULT_ENV=gtotree
ENV PATH=/opt/conda/envs/$CONDA_DEFAULT_ENV/bin:/GToTree-1.1.3/bin:$PATH
ENV GToTree_HMM_dir=/GToTree-1.1.3/hmm_sets

# test 
#WORKDIR /GToTree-1.1.3/test_data
#RUN GToTree -a ncbi_accessions.txt -g genbank_files.txt -f fasta_files.txt -H Bacteria -m genome_to_id_map.tsv -t -j 4 -o GToTree_test

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]
