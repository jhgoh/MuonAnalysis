#!/bin/bash

for i in Loose Medium Tight2012 LooseWithRPC LooseDefault RPC{"",LSt,SSt,TSt}{Loose,Tight}; do
    cmsRun tp_fit_RD_cfg.py $i
done
