#!/bin/bash

for kl in {0.5,1,2,4,6};do
    for a in {0.001,0.01,0.1};do
        echo kl = ${kl}, a = ${a}
        python wavenumber.py $kl $a 60
    done
done
