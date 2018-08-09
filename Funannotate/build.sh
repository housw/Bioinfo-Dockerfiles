#!/bin/bash

# copy source files
cp /home/shengwei/software/src/GeneMark/gm_key_64.gz .
cp /home/shengwei/software/src/GeneMark/gm_et_linux_64.tar.gz .
cp /home/shengwei/software/src/SignalP/signalp-4.1f.Linux.tar.gz .
cp /home/shengwei/software/src/RepeatMasker/RepBaseRepeatMaskerEdition-20170127.tar.gz .
cp /home/shengwei/software/src/RepeatMasker/nseg.tar.gz .
cp /home/shengwei/software/src/Funannotate/funannotate-1.0.0.tar.gz .

# build 
docker build -t shengwei/funannotate -f Dockerfile .

