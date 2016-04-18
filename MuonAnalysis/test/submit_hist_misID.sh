#!/bin/bash

if [ ! -d RD ]; then
  mkdir -p RD/ks
  mkdir -p RD/phi
  mkdir -p RD/lamb
fi
if [ ! -d MC ]; then
  mkdir -p MC/ks
  mkdir -p MC/phi
  mkdir -p MC/lamb
fi

SUBMIT='bsub -q 1nh'

for i in `seq 1 161`; do 
    $SUBMIT -oo RD/ks/log_$i.log run_hist_misID.sh RD ks $i
    $SUBMIT -oo RD/phi/log_$i.log run_hist_misID.sh RD phi $i
    $SUBMIT -oo RD/lamb/log_$i.log run_hist_misID.sh RD lamb $i
done

for i in `seq 1 628`; do 
    $SUBMIT -oo MC/ks/log_$i.log run_hist_misID.sh MC ks $i
    $SUBMIT -oo MC/phi/log_$i.log run_hist_misID.sh MC phi $i
    $SUBMIT -oo MC/lamb/log_$i.log run_hist_misID.sh MC lamb $i
done

