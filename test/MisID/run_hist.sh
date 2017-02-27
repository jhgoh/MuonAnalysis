#!/bin/bash

cd $CMSSW_BASE/src/MuonAnalysis/MuonIdentification/test/MisID
eval `scram runtime -sh`

python hist.py $1 $2 $3
