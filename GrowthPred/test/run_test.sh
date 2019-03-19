#!/bin/bash


GrowthPred_CONTAINER="my_growthpred"
INPUT_DIR="input_dir"
OUTPUT_DIR="output_dir"
mkdir -p ${OUTPUT_DIR}

# initialize a container in background
docker run -t -d --rm \
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
           -u $(id -u $USER):$(id -g $USER) \
           -v $PWD:/mnt \
           --name ${GrowthPred_CONTAINER} shengwei/growthpred:latest

# with heg
cmd_test1="/growthpred/growthpred-v1.08.py -d ${INPUT_DIR} -f ecoli_ribosomal_genes.txt -g ecoli_complete_genome.txt \
	   -s -S -c 0 -t -m -o ${OUTPUT_DIR}/test_with_heg"
echo ${cmd_test1}
docker exec ${GrowthPred_CONTAINER} /bin/bash -c "$cmd_test1"

# without heg
cmd_test2="/growthpred/growthpred-v1.08.py -d ${INPUT_DIR} -g ecoli_complete_genome.txt \
           -s -S -c 0 -t -m -b -o ${OUTPUT_DIR}/test_without_heg"
echo ${cmd_test2}
docker exec ${GrowthPred_CONTAINER} /bin/bash -c "$cmd_test2"

# kill container
docker container kill ${GrowthPred_CONTAINER}



