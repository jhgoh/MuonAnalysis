#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: ./runFit.sh {MC,RD}"
    exit 1
fi

for i in Loose Medium Tight2012 LooseWithRPC LooseDefault RPC{"",LSt,SSt,TSt}{Loose,Tight}; do
    cmsRun tp_fit_${1}_cfg.py $i
done
