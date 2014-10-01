#! /bin/bash

CS_FILE=/home/skimmy/data/ecoli.solid/ecoli_F3_200k_sample.csfasta
QUAL_FILE=/home/skimmy/data/ecoli.solid/ecoli_F3_200k_sample_QV.qual
GENOME_FILE=/home/skimmy/data/ecoli.solid/DH10B_NCBI.fna

./seq $CS_FILE $QUAL_FILE $GENOME_FILE
