#!/bin/bash

CONTAINER="test_KofamKOALA"
profiles_dir="/kraken_home/shengwei/GitHub/docker/KofamKOALA/profiles"
ko_list_file="/kraken_home/shengwei/GitHub/docker/KofamKOALA/ko_list"

# initialize a container in background
docker run -t -d --rm \
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
           -u $(id -u $USER):$(id -g $USER) \
           -v ${profiles_dir}:/db/profiles \
	   -v ${ko_list_file}:/db/ko_list \
           -v $PWD:/mnt:rw \
           --name ${CONTAINER} shengwei/kofamkoala:latest

cmd="/kofamscan/bin/exec_annotation -o test_ko.txt --profile /db/profiles --ko-list /db/ko_list --cpu 10 --format mapper test.faa"
echo $cmd
docker exec $CONTAINER /bin/bash -c "$cmd"

# kill the container
docker kill $CONTAINER


