#!/bin/bash

rm -r results.txt
#for M in `seq 500 25 1000`
for M in `seq 1 200`
do
    #./simulator.out -N 10000 -m 100 -M $M -e 0 | grep -i fail >> results.txt
    ./simulator.out -N 10000 -m 100 -M 1000 -e 0 | grep -i fail >> results.txt
done
