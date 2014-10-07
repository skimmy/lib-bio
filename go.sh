#!/bin/bash

BASE_IN_DIR=/home/skimmy/data/ecoli.solid
BASE_OUT_DIR=/home/skimmy/data/gpu

INPUT="-g /home/skimmy/data/ecoli.solid/genome.dat -f 0 -r /home/skimmy/data/ecoli.solid/reads.dat -F 0"
OUTPUT="-G /home/skimmy/data/ecoli.solid/genome.dat -R /dev/null -o /home/skimmy/data/gpu/align.out"
PREPROC="-p 8 -c 8"
OPERATION="-A 0 -k 3"

./seq.o ${INPUT} ${PREPROC} ${OPERATION} ${OUTPUT}
