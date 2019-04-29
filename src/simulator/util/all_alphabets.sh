#!/bin/bash

declare -a sigmas=("01" "0123")

#declare -a sigmas=("01" "012" "0123" "01234" "012345" "01234567"
#		   "0123456789ABCDEF" "0123456789ABCDEFabcdefghijklmnop")

declare -a Ns=(1024 2048)

#declare -a Ns=(16 32 64 128 256 512 1024 2048 4096 8192 16384)

k=5000
printf "sigma,N,s_mean,s_var\n"
for sigma in "${sigmas[@]}"
do
    s=${#sigma}
    for N in "${Ns[@]}"
    do
	printf "%s,%d," ${s} ${N}
	printf "sigma=%s\tN=%d\n" ${s} ${N} >&2
	time ./simulator.out -O 6 -N ${N} -k ${k} -a ${sigma} -v 1 -V "outputs/verbose_sigma${s}_n${N}_k${k}.txt" 2>/dev/null
	printf "\n" >&2
    done
done

