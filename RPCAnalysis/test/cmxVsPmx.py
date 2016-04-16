#!/usr/bin/env python
from ROOT import *
from math import *
from SKKU.CommonTools.tdrStyle import *
#gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)

#wget 'https://github.com/hep-skku/SKKU/blob/49e3f08cd002f7df016e8e25256002d35250e323/RPCAnalysis/test/CMX.root?raw=true -O CMX.root'
#wget 'https://github.com/hep-skku/SKKU/blob/49e3f08cd002f7df016e8e25256002d35250e323/RPCAnalysis/test/PMX_DET.root?raw=true -O PMX_DET.root'
#wget 'https://github.com/hep-skku/SKKU/blob/49e3f08cd002f7df016e8e25256002d35250e323/RPCAnalysis/test/PMX_NONDET.root?raw=true -O PMX_NONDET.root'

refInput = ["classic", kBlack, 1, TFile("CMX.root")]
cmpInputs = [
  ["premix-det", kBlue, 1, TFile("PMX_DET.root")],
  ["premix-nondet", kRed, 1, TFile("PMX_NONDET.root")]
]

def getHistNames(tdir):
    hNames = []
    for key in tdir.GetListOfKeys():
        keyName = key.GetName()
        obj = tdir.Get(keyName)
        if obj.IsA().InheritsFrom("TH1"):
            path = '/'.join(tdir.GetPath().split('/')[1:])
            hNames.append(path+'/'+obj.GetName())
        elif obj.IsA().InheritsFrom("TDirectory"):
            hNames.extend(getHistNames(obj))
        else:
            print obj.IsA().GetName()
    return hNames

hNames = getHistNames(refInput[3])

allObjs = []
for i, hName in enumerate(hNames):
    gROOT.cd()

    c = TCanvas("c%d" % i, hName, 500, 600)
    pad1 = TPad("pad%d_1" % i, hName, 0, 0.3, 1, 1)
    pad2 = TPad("pad%d_2" % i, hName, 0, 0, 1, 0.3)
    pad1.SetMargin(0.15, 0.08, 0.05, 0.15)
    pad2.SetMargin(0.15, 0.08, 0.2, 0)
    pad1.SetGridx()
    pad1.SetGridy()
    pad2.SetGridx()
    pad2.SetGridy()

    hRef = refInput[3].Get(hName).Clone()
    sumRef = hRef.Integral()

    pad1.cd()
    hRef.SetMarkerSize(0.7)
    hRef.SetMarkerStyle(20)
    hRef.SetLineColor(refInput[1])
    hRef.SetMarkerColor(refInput[1])
    hRef.GetXaxis().SetLabelSize(0.03)
    hRef.Draw("e")

    fRatio = hRef.Clone()
    fRatio.Reset()
    fRatio.SetTitle("")
    fRatio.SetMinimum(0.8)
    fRatio.SetMaximum(1.2)
    fRatio.GetXaxis().SetLabelSize(0.10)
    fRatio.GetYaxis().SetLabelSize(0.10)
    fRatio.GetYaxis().SetNdivisions(505)
    pad2.cd()
    fRatio.Draw()

    leg = TLegend(0.65, 0.6, 0.9, 0.8)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(hRef, refInput[0], "lp")

    objs = [c, pad1, pad2, hRef, fRatio, leg]

    for cmpInput in cmpInputs:
        fCmp = cmpInput[3]
        hCmp = fCmp.Get(hName)
        hCmp.Scale(sumRef/hCmp.Integral())

        grpRatio = TGraphErrors()
        for b in range(1, hRef.GetNbinsX()):
            x = hRef.GetBinCenter(b)
            ex = hRef.GetBinWidth(b)/2
            y1, ey1 = hRef.GetBinContent(b), hRef.GetBinError(b)
            y2, ey2 = hCmp.GetBinContent(b), hCmp.GetBinError(b)
            r, er = 0, 0
            if y1 > 0: r = y2/y1
            if y1 > 0 and y2 > 0: er = sqrt( (ey1/y1)**2 + (ey2/y2)**2 )*r

            grpRatio.SetPoint(b-1, x, r)
            grpRatio.SetPointError(b-1, ex, er)

        pad1.cd()
        hCmp.SetMarkerSize(0.7)
        hCmp.SetMarkerStyle(20)
        hCmp.SetMarkerColor(cmpInput[1])
        hCmp.SetLineColor(cmpInput[1])
        hCmp.Draw("samee")
        pad2.cd()
        grpRatio.SetMarkerSize(0.7)
        grpRatio.SetMarkerStyle(20)
        grpRatio.SetMarkerColor(cmpInput[1])
        grpRatio.SetLineColor(cmpInput[1])
        grpRatio.Draw("ep")

        leg.AddEntry(hCmp, cmpInput[0], "lp")

        objs.extend([grpRatio, hCmp])

    pad1.cd()
    leg.Draw()

    pad1.Modified()
    pad2.Modified()
    c.Modified()

    c.cd()
    pad1.Draw()
    c.cd()
    pad2.Draw()

    c.Update()
    c.Draw()

    allObjs.extend(objs)

