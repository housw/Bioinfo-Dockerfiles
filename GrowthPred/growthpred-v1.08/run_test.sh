#!/bin/bash

mkdir -p test_output
./growthpred-v1.08.py -d ./ -f ecoli_ribosomal_genes.txt -g ecoli_complete_genome.txt -s -S -c 0 -o test_output/test_ecoli -t -m 
#./growthpred-v1.08.py -d ./ -g ecoli_complete_genome.txt -s -S -c 0 -o test_ecoli -t -m -b
#./growthpred-v1.08.py -f PTR0Bin1Ribosomal.fasta -g PTR0Bin1nucgenes.fna -s -S -c 0 -o test_PTR0Bin1 -t -m
#./growthpred-v1.08.py -f PTR1BBin29ribo.fasta -g PTR1BBin29nucgenes.fna -s -S -c 0 -o test_PTR1Bin29 -t -m


