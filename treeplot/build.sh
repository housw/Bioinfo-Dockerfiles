#!/bin/bash 

# download mview
#:<<'COMMENT'
wget https://sourceforge.net/projects/bio-mview/files/bio-mview/mview-1.64/mview-1.64.tar.gz && \
tar zxfv mview-1.64.tar.gz && rm mview-1.64.tar.gz && chmod -R 777 mview-1.64 && \
sed -i -e 's/\/path\/to\/mview\/lib/\/mview\/lib/g' mview-1.64/bin/mview
#COMMENT

# build docker image
docker build -t 'shengwei/treeplot:latest' .
