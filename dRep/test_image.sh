#!/bin/bash

CONTAINER="dRep"
INPUT_BINs="genomes"
OUTPUT_DIR=${INPUT_BINs}_dRep


# initialize a container in background
docker run -t -d --rm \
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
           -u $(id -u $USER):$(id -g $USER) \
           -v $PWD:/mnt:rw \
           --name ${CONTAINER} shengwei/drep:latest

# 1) run compare 
cmd1="dRep compare /mnt/$OUTPUT_DIR/compare -g /mnt/$INPUT_BINs/*.fasta"
echo $cmd1
docker exec $CONTAINER /bin/bash -c "$cmd1"

# 2) dereplicate
cmd2="dRep dereplicate $OUTPUT_DIR/derep -g $INPUT_BINs/*.fasta"
echo $cmd2
docker exec $CONTAINER /bin/bash -c "$cmd2"

# kill the container
docker kill $CONTAINER

