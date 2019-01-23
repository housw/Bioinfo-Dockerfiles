#!/bin/bash

# clone esom and fix the interpreter
git clone https://github.com/CK7/esom.git
sed --in-place -e 's/\/usr\/bin\/perl/\/usr\/bin\/env\ perl/g' esom/prepare_esom_files.pl

# build docker image
docker build -t 'shengwei/abawaca:1.0.7' .

# clean up
rm -rf esom
