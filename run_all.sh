#!/bin/bash

for kl in {1,2,4,6};do
    for a in {0.05,0.075,0.1};do
        echo kl = ${kl}, a = ${a}
        python wavenumber.py $kl $a 110 &
    done
done
