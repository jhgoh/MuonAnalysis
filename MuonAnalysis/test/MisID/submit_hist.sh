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

for i in `seq 1 161`; do
    bsub -q 1nh -oo RD/ks/log_$i.log run_hist.sh RD ks $i
    bsub -q 8nh -oo RD/phi/log_$i.log run_hist.sh RD phi $i
    bsub -q 1nh -oo RD/lamb/log_$i.log run_hist.sh RD lamb $i
done

for i in `seq 1 628`; do
    bsub -q 1nh -oo MC/ks/log_$i.log run_hist.sh MC ks $i
    bsub -q 8nh -oo MC/phi/log_$i.log run_hist.sh MC phi $i
    bsub -q 1nh -oo MC/lamb/log_$i.log run_hist.sh MC lamb $i
done

