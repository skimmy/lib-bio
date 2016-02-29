#!/bin/bash



N=$1
m=$2
e=0.01
#coverage=`seq 1 20`

for c in `seq 5 5 20`
do
    M=$(((N*c)/m))
    ./run100.sh $N $M $m $e
    mean=`python mean.py results100.txt`
    echo -e "${m}\t${c}\t${mean}" >> results.txt
    cp results100.txt results100_${m}_${c}.txt
done


#for M in `seq 500 25 1000`
#for M in `seq 1 200`
#do
    #./simulator.out -N 10000 -m 100 -M $M -e 0 | grep -i fail >> results.txt
#    ./simulator.out -N 10000 -m 100 -M 1000 -e 0 | grep -i fail >> results.txt
#done
