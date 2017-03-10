#!/usr/bin/env python

from ROOT import *
from glob import glob
import os

gROOT.SetBatch()
TProof.Open("")
  
basedir = "20170303_1"

for sample in ["JetHT_2015D", "TT_powheg", "QCD"]:
    outdir = "%s/%s" % ("hist", sample)
    if not os.path.exists(outdir): os.makedirs(outdir)
    files = glob("%s/%s*/*.root" % (basedir, sample))

    for mode in ["ks", "phi", "lamb", "Bp"]:
        chain = TChain("%s/tree" % mode)
        for f in files: chain.Add(f)
        chain.SetProof();
        chain.Process("Selector.C+", mode, 1000)
    os.system("mv *.root %s" % outdir)

