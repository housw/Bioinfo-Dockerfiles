#!/bin/bash

# download priam and database
wget http://priam.prabi.fr/utilities/PRIAM_search.jar  
wget http://priam.prabi.fr/REL_JAN18/Distribution.zip
unzip Distribution.zip -d Distribution
rm Distribution.zip

# test run to build database 
source activate blast-legacy
cd Distribution
threads=5
session_name=priam_job
input_faa=../sample.faa
java -jar -Xmx50g ../PRIAM_search.jar -np $threads -n $session_name -i $input_faa -p PRIAM_JAN18 -od priam_output_dir --min_proba 0.5 --min_proportion 70 --check_catalytic
cd ../
source deactivate

# change permission
chmod -R 777 PRIAM_search.jar Distribution/

# build the docker image 
docker build -t 'shengwei/priam:latest' .

# clean up
rm -rf Distribution/ PRIAM_search.jar 
