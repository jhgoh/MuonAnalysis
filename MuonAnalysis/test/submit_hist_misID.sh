#!/bin/bash

SUBMIT='bsub -q 1nh run_hist_misID.sh'
#SUBMIT='./run_hist_misID.sh'

for i in `seq 1 161`; do 
    $SUBMIT RD ks $i
    $SUBMIT RD phi $i
    $SUBMIT RD lamb $i
done

for i in `seq 1 628`; do 
    $SUBMIT MC ks $i
    $SUBMIT MC phi $i
    $SUBMIT MC lamb $i
done

