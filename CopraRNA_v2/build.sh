#!/bin/bash

# update supported organisms
git clone https://github.com/PatrickRWright/CopraRNA.git && \
cd CopraRNA/update_kegg2refseq && mkdir -p run && cd run && ../build_kegg2refseq.pl 

# build image
docker build -t 'shengwei/coprarna:2.1.2' .

# clean up 
rm -rf CopraRNA
