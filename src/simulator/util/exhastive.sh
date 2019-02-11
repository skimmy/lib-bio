#!/bin/bash

declare -a sigmas=("01" "0123" "0123" "0123456789abcdef")

declare -a Ns=(4 5 6 7 8 9)

printf "n,e(n),alpha(n)\n"
for sigma in "${sigmas[@]}"
do
    s=${#sigma}
    for N in "${Ns[@]}"
    do
	time ./simulator.out -O 6 -f 1 -N ${N} -t 2 -a ${sigma} 2>/dev/null
    done
done

