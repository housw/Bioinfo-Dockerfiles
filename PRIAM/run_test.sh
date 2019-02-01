#!/bin/bash

CONTAINER="priam"
INPUT_FA="sample.fasta"
FILESTEM=${INPUT_FA%.fasta}


# initialize a container in background
docker run -t -d --rm \
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
           -u $(id -u $USER):$(id -g $USER) \
           -v $PWD:/mnt:rw \
           --name ${CONTAINER} shengwei/priam:latest

# 1) split scaffolds into two sets 
cmd1="java -jar -Xmx4096m PRIAM_search.jar -np 4  --in /mnt/sample.fasta --priam /PRIAM/PRIAM_JAN18 --out /mnt/priam_test_output --job_name $FILESTEM --min_proba 0.5 --min_proportion 70 --check_catalytic"
echo $cmd1
docker exec $CONTAINER /bin/bash -c "$cmd1"

docker kill $CONTAINER
