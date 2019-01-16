#!/bin/bash

#:<<'COMMENT'
git clone https://github.com/ParkinsonLab/Metabolic-pathways-in-microbiomes.git
mv Metabolic-pathways-in-microbiomes/EC\ annotation Metabolic-pathways-in-microbiomes/EC_annotation
mv Metabolic-pathways-in-microbiomes/Statistics\ analysis Metabolic-pathways-in-microbiomes/Statistics_analysis

# convert mac linebreaks to unix style
for file in Metabolic-pathways-in-microbiomes/EC_annotation/*.py ; do
    vi +':w ++ff=unix' +':q' ${file}
done
for file in Metabolic-pathways-in-microbiomes/Statistics_analysis/*.R ; do
    vi +':w ++ff=unix' +':q' ${file}
done

# fix the input bug in extract_diamond.py
sed --in-place -e 's/sys.argv\[4\]/sys.argv\[3\]/g' Metabolic-pathways-in-microbiomes/EC_annotation/extract_diamond.py
#COMMENT

docker build -t 'shengwei/diamondenzyme:latest' .

rm -rf Metabolic-pathways-in-microbiomes
