#!/bin/bash

cd /afs/cern.ch/user/j/jhgoh/work/MuonPOG/MisID/CMSSW_7_6_3_patch2/src/SKKU/MuonAnalysis/test/MisID
eval `scram runtime -sh`

python hist.py $1 $2 $3
