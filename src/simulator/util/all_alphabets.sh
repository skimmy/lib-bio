#!/bin/bash

#declare -a sigmas=("01" "012")

declare -a sigmas=("01" "012" "01234" "012345" "01234567"
		   "0123456789ABCDEF" "0123456789ABCDEFabcdefghijklmnop")

declare -a Ns=(1024 2048 4096 8192 16384)

k=5000
printf "sigma,N,s_mean,s_var\n"
for sigma in "${sigmas[@]}"
do
    s=${#sigma}
    for N in "${Ns[@]}"
    do
	printf "%s,%d," ${s} ${N}
	printf "sigma=%s\tN=%d" ${s} ${N} >&2
	time ./simulator.out -O 6 -N ${N} -k ${k} -a ${sigma} 2>/dev/null
	printf "\n" >&2
    done
done

