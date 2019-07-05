
## initialize a container named dastool in background
```
docker run -t -d --rm 
           -v /etc/group:/etc/group:ro -v /etc/passwd:/etc/passwd:ro 
           -u $(id -u $USER):$(id -g $USER) 
           -v $PWD:/mnt 
           --name dastool shengwei/das_tool:latest
```

## run das_tool 

#### copy sample_input and sample_output to /mnt

The default working directory is set to /mnt, which is the location your local working directory mount to. 

```
docker exec dastool /bin/bash -c "cp -r /DAS_Tool/current/sample_* ." 
```

#### run sample commands
```
# sample 1
cmd1="DAS_Tool -i sample_data/sample.human.gut_concoct_scaffolds2bin.tsv,sample_data/sample.human.gut_maxbin2_scaffolds2bin.tsv,sample_data/sample.human.gut_metabat_scaffolds2bin.tsv,sample_data/sample.human.gut_tetraESOM_scaffolds2bin.tsv -l concoct,maxbin,metabat,tetraESOM -c sample_data/sample.human.gut_contigs.fa -o sample_output/DASToolRun1
         --write_bin_evals 1 
         --create_plots 1 
         --write_bins 1 
         --threads 20"
echo $cmd1
docker exec dastool /bin/bash -c "$cmd1"

# sample2 
CMD2="DAS_Tool -i sample_data/sample.human.gut_concoct_scaffolds2bin.tsv,sample_data/sample.human.gut_maxbin2_scaffolds2bin.tsv,sample_data/sample.human.gut_metabat_scaffolds2bin.tsv,sample_data/sample.human.gut_tetraESOM_scaffolds2bin.tsv -l concoct,maxbin,metabat,tetraESOM -c sample_data/sample.human.gut_contigs.fa -o sample_output/DASToolRun2 --threads 2 --score_threshold 0.6 --proteins sample_output/DASToolRun1_proteins.faa
         --write_bin_evals 1 
         --create_plots 1 
         --write_bins 1 
         --threads 20"
docker exec dastool /bin/bash -c "$cmd2"
```

## kill container
docker container kill dastool

## Disclaimer: 

The `usearch` might not work properly with your account, please request a linux version `usearch` from [here](https://www.drive5.com/usearch/download.html), rename it to `usearch` and put it into `/usr/local/bin`. Or you can build a new image from scratch with your own `usearch`, the Dockerfile is [here](https://github.com/housw/docker/tree/master/DAS_Tool). 
