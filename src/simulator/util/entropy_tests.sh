#!/bin/bash

for i in `seq 1 40`
do
    python entropy.py $i 0.25
done
