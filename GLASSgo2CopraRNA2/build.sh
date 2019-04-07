#!/bin/bash

BASE_DIR="$PWD"

# update supported organisms
#git clone https://github.com/PatrickRWright/CopraRNA.git && \
#cd CopraRNA/update_kegg2refseq && mkdir -p run && cd run && ../build_kegg2refseq.pl && \
#chmod 777 CopraRNA_available_organisms.txt && cp CopraRNA_available_organisms.txt ${BASE_DIR}/scripts && \
#cd ${BASE_DIR} && rm -rf CopraRNA

# update taxid_to_refseq
#cd scripts && ./update_taxid_to_refseq.r

# build image
cd ${BASE_DIR}
chmod -R 777 scripts
docker build -t 'shengwei/glassgo2coprarna2:1.5.1b' .
