#!/bin/bash

for i in RPC{"",LSt,SSt,TSt}{Loose,Tight}; do
    echo cmsRun tp_fit_RD_cfg.py $i
done
