#!/usr/bin/env python
from ROOT import *
gROOT.ProcessLine(".x CMS/rootlogon.C")
#gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)

inputs = [
  ("classic", kRed, 1, TFile("CMX.root")),
  ("premix-nondet", kBlue, 1, TFile("PMX_NONDET.root"))
]

valDir = ["rpcRecHitAnalysis", "rpcRecHitAnalysis/Bx0", "rpcRecHitAnalysis/BxX"]

dirList = []
histList = []
histNameList = []
for dir in valDir:
  for k in inputs[0][-1].Get(dir).GetListOfKeys():
    if "h" in k.GetName():
      dirList.append(dir)
      h = inputs[0][-1].Get(dir+"/"+k.GetName())
      histList.append(h.GetName())
      if "Bx0" in dir:
        h.SetName("Bx0-%s" % h.GetName())
      if "BxX" in dir:
        h.SetName("BxX-%s" % h.GetName())
      histNameList.append(h.GetName())

gROOT.cd()
allObjs = []
print histList

for j,histName in enumerate(histList):
  c = TCanvas("c%s" % histNameList[j], "c%s" % histNameList[j], 500, 800)
  pad1 = TPad("p1%s" % histName,"p1%s" % histNameList[j], 0, 0.3, 1, 1)
  pad2 = TPad("p2%s" % histName,"p2%s" % histNameList[j], 0, 0, 1, 0.3)

  padList = [pad1, pad2]
  leg = TLegend(0.6, 0.65, 0.9, 0.85,"")

  for pad in padList:
    pad.SetMargin(0.05, 0.05, 0.05, 0.05)
    pad.SetGridy()
    pad.Draw() 
  padList[0].cd() 
  isFirst = True
  for item in inputs:
    opt = "samee"
    h = item[-1].Get(dirList[j]+"/"+histName)
    h.SetTitle(histNameList[j])
    h.SetLineColor(item[1])
    h.SetLineStyle(item[2])
    h.SetLineWidth(3)
    leg.AddEntry(h, item[0])
    if isFirst:
      h.Draw("e")
      hR = h.Clone()
      hR.Sumw2()
      hR.SetLineColor(kBlack)
      hR.SetMarkerStyle(20)
      hR.SetTitle("")
      hIntegral = h.Integral()
      isFirst = False
    if not isFirst:
      h.Scale(hIntegral/h.Integral())
      h.Draw("samee")
      hR.Divide(h)
  leg.Draw()
  hR.SetMaximum(1.8)   
  hR.SetMinimum(0.2) 
  padList[1].cd()
  hR.Draw("ep")
  c.SaveAs("%s.png" % c.GetName())

  allObjs.extend([h,c])
           
