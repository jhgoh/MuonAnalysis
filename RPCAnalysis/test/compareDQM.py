#!/usr/bin/env python

import sys, os
if len(sys.argv) < 3:
    print "Usage: compareRoot.py File1.root File2.root"
    sys.exit(1)

def fetchDirStructure(f, dirName="/"):
    histNames = []
    tdir = f.GetDirectory(dirName)
    if tdir == None: return histNames

    for key in [x.GetName() for x in tdir.GetListOfKeys()]:
        pathName = dirName+"/"+key

        obj = f.Get(pathName)
        if obj == None: continue

        if obj.InheritsFrom("TDirectory"):
            histNames.extend(fetchDirStructure(f, pathName))
        elif obj.InheritsFrom("TH1"):
            histNames.append(pathName)

    return histNames

from ROOT import *
f1 = TFile(sys.argv[1])
f2 = TFile(sys.argv[2])

hNames1 = fetchDirStructure(f1, "DQMData/Run 1/RPC/Run summary")
hNames2 = fetchDirStructure(f2, "DQMData/Run 1/RPC/Run summary")

commonHists= set(hNames1).intersection(set(hNames2))
print "Comparing %d histograms" % len(commonHists)
diffHists = []
for hPath in commonHists:
    h1 = f1.Get(hPath)
    h2 = f2.Get(hPath)

    ## Check difference between two histograms.
    isDiff = False
    for b in range(h1.GetNbinsX()):
        y1 = h1.GetBinContent(b)
        y2 = h2.GetBinContent(b)
        if y1 != y2:
            isDiff = True
            break
    if isDiff: diffHists.insert(hPath)

if len(diffHists) == 0:
    print "Identical :-)"
else:
    print "Got %d histograms with different bin contents" % len(diffHists)
    objs = {}
    for i, hPath in enumerate(diffHists):
        if h1.InheritsFrom("TH2"):
            c = TCanvas("c%d" % i, hPath, 800, 400)
            c.Divide(2,1)
            c.cd(1)
            h1.Draw()
            c.cd(2)
            h2.Draw()
        else:
            c = TCanvas("c%d" % i, hPath, 400, 400)
            h1.SetMaximum(max(h1.GetMaximum(), h2.GetMaximum())*1.1)
            h1.SetLineColor(kRed)
            h2.SetLineColor(kBlue)
            h1.Draw()
            h2.Draw("same")

        objs[hPath] = [c]

