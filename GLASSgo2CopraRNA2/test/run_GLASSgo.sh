#!/bin/bash

# ----------
# parameters
# ----------

sRNA_DIR="00_sRNA_folder"
nt_DIR="/mnt/data/db/blastdb/nt/"
acc_list="acc_list/Firmicutes.acc"
Wildcard_RefSeq="NC_007795,NC_004461,NC_013893,NC_020164,NZ_CP007601,NC_014925,NC_007350"
Interested_RefSeq="NZ_CP018205"
CONTAINER="glassgo2coprarna2"


# initialize a container in background
docker run -t -d --rm \
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro \
           -u $(id -u $USER):$(id -g $USER) \
           -v $nt_DIR:/db/nt/ \
           -v $PWD:/mnt:rw \
           --name ${CONTAINER} shengwei/glassgo2coprarna2:1.5.1a # shengwei/glassgo2coprarna2:1.5.1


# ------------------------------
# GLASSgo sRNA homolog searching
# ------------------------------

GLASSgo_OUT_FOLDER='01_GLASSgo_Results'
mkdir -p $GLASSgo_OUT_FOLDER

:<<'COMMENT'
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
COMMENT


# --------------------
# run GLASSgo2CopraRNA 
# --------------------

GLASSgo2CopraRNA="02_GLASSgo2CopraRNA"
mkdir -p $GLASSgo2CopraRNA

:<<"COMMENT"
for f in `ls $GLASSgo_OUT_FOLDER`; do
    if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
        echo "Running GLASSgo2CopraRNA for sRNA: " $f
        fullpath=$GLASSgo_OUT_FOLDER/$f
        filestem=${f%.*}
        ID=${filestem#GLASSgo_output_}
        # fasta generation
        cmd1="R --slave -f  /GLASSgo/GLASSgo2CopraRNA_fasta_generation_by_accession.r --args input_file=$fullpath refpath=/GLASSgo/taxid_to_refseq cop_path=/GLASSgo/CopraRNA_available_organisms.txt output_file=${ID}_coprarna_candidates.txt"
        echo $cmd1
        docker exec $CONTAINER /bin/bash -c "$cmd1"
        # exclusion 
        cmd2="R --slave -f /GLASSgo/GLASSgo2CopraRNA_exclusion_script.r --args datapath=full_GLASSgo_table.Rdata wildcard=${Wildcard_RefSeq} ooi=${Interested_RefSeq}"
        echo $cmd2
        docker exec $CONTAINER /bin/bash -c "$cmd2"
        # balanced selection
        cmd3="R --slave -f /GLASSgo/GLASSgo2CopraRNA_balanced_ooi_selection.r --args wildcard=${Wildcard_RefSeq} ooi=${Interested_RefSeq} max_number=20 outfile_prefix=${ID}_sRNA"
        echo $cmd3
        docker exec $CONTAINER /bin/bash -c "$cmd3"
        # move results to its folder
        mkdir -p $GLASSgo2CopraRNA/$ID
        mv ${ID}_* $GLASSgo2CopraRNA/$ID
        # clean up
        rm full_GLASSgo_table.Rdata refined_GLASSgo_table.Rdata temp_fasta_wo_wild_w_ooi.fasta temp_fasta_wo_wild.fasta un_fasta 
    fi
done
COMMENT


# kill the container
docker kill $CONTAINER

