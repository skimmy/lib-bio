#!/bin/bash

k=$1
P=$2
c=$3

out_file="/tmp/ed_sample.out"
rm -f "${out_file}"

for N in `seq 11 100`
do
    printf "%s\t" ${N} >> "${out_file}"
    ../simulator.out -N ${N} -O 6 -f 32 -k $k -P $P -c $c >> ${out_file}
done


