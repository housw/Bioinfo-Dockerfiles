#!/bin/bash

CONTAINER="GLASSgo2"
sRNA_DIR="00_sRNA_folder"
GLASSgo_OUT_FOLDER='01_GLASSgo_Results'
mkdir -p $GLASSgo_OUT_FOLDER

# parameters
nt_DIR="/mnt/data/db/blastdb/nt/"
acc_list="acc_list/Firmicutes.acc"
CONTAINER="glassgo2coprarna2"

# initialize a container in background
docker run -t -d --rm \
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
           -u $(id -u $USER):$(id -g $USER) \
           -v $nt_DIR:/db/nt/ \
           -v $PWD:/mnt:rw \
           --name ${CONTAINER} shengwei/glassgo2coprarna2:1.5.1


for sRNA in `ls $sRNA_DIR/*.fa`; do
    echo "Running GLASSgo for sRNA: " $sRNA
    BASENAME=$(basename $sRNA)
    output_sRNAs=GLASSgo_output_${BASENAME}
    output_json=${output_sRNAs%.*}.json
    # run glassgo
    cmd1="/GLASSgo/GLASSgo.py --in_seek_seq /mnt/$sRNA --db /db/nt/nt --eValue 1 --acc_pos /mnt/$acc_list --num_threads 40 --output_name $output_sRNAs"
    echo $cmd1
    docker exec $CONTAINER /bin/bash -c "$cmd1"
    # generate json
    cmd2="/GLASSgo/createJSON-1_0.py -g $output_sRNAs -o $output_json"
    echo $cmd2
    docker exec $CONTAINER /bin/bash -c "$cmd2"
    # genereate html
    cmd3="/GLASSgo/generate_GLASSgo_html.py $output_json"
    echo $cmd3
    docker exec $CONTAINER /bin/bash -c "$cmd3"
done
mv GLASSgo_output_* $GLASSgo_OUT_FOLDER

# kill the container
docker kill $CONTAINER

