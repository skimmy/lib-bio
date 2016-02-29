#!/bin/bash

N=$1
M=$2
m=$3

rm -f results100.txt
for i in `seq 1 100`
do
    ./simulator.out -N ${N} -M ${M} -m ${m} -e 0.01 -o -p >> results100.txt
done
