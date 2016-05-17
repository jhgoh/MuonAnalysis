#!/usr/bin/env python

modeDef = {
    "proton":("lamb", "Proton"),
    "pion":("ks", "Pion"),
    "Kp":("phi", "K^{+}"),
    "Km":("phi", "K^{-}"),
}

import sys, os
if len(sys.argv) < 3 or \
    sys.argv[1] not in ("MC", "RD") or \
    sys.argv[2] not in modeDef:
    print "python draw.py [RD,MC] [proton,pion,Kp,Km]"
    os.exit(1)

dataType = sys.argv[1]
submod = sys.argv[2]

from math import *
from ROOT import *
from SKKU.CommonTools.tdrStyle import *
tdrStyle.SetTitleSize(0.055, "XYZ")
tdrStyle.SetPadTopMargin(0.1)
tdrStyle.SetTitleXOffset(1.15)
tdrStyle.SetTitleYOffset(1.15)

mode, xtitle = modeDef[submod]

varNames = {
    "pt":("%s transverse momentum p_{T} (GeV)" % xtitle),
    "abseta":("%s pseudorapidity |#eta|" % xtitle)
}
idSets = [
    ("RPC","RPC Loose ID", kRed),
    ("TRK","Tracker muons (no arbitration)", kBlue),
    ("Loose","Loose muon ID", kGreen+1),
    ("Medium","Medium muon ID", kMagenta),
    ("Tight","Tight muon ID", kBlack),
    ("Soft","Soft muon ID", kAzure+1),
    #("OneLoose","TMOneStationLoose", kBlue),
    #("OneTight","TMOneStationTight", kGreen+1),
    #("TMLastLoose","TMLastStationLoose", kGreen+1),
    #("TMLastTight","TMLastStationTight", kMagenta),
    #("TM2DLoose","TM2DLoose" kBlack),
    #("TM2DTight","TM2DTight", kAzure+1),
]

frames = {}
graphs = {}

maxPt = 20
maxY = 0
f = TFile("%s/fit_%s.root" % (dataType, mode))
#for idName in [x.GetName() for x in f.GetDirectory(submod).GetListOfKeys()]:
for i, (idName, idTitle, color) in enumerate(idSets):
    idDir = f.GetDirectory("%s/%s" % (submod, idName))
    if idDir == None: continue
    for varName in [x.GetName() for x in idDir.GetListOfKeys()]:
        varDir = idDir.GetDirectory(varName)
        if varDir == None: continue
        if varName not in frames:
            frame = varDir.Get("hFrame")
            frame.SetTitle(";%s;Misidentification probability (%%)" % varNames[varName])
            frames[varName] = frame
            graphs[varName] = []
        g = varDir.Get("gRatio")
        g.SetTitle(idTitle)
        g.SetMarkerColor(color)
        g.SetLineColor(color)
        for p in xrange(g.GetN()):
            g.SetPoint(p, g.GetX()[p], g.GetY()[p]*100)
            g.SetPointError(p, g.GetEXlow()[p], g.GetEXhigh()[p], g.GetEYlow()[p]*100, g.GetEYhigh()[p]*100)
            maxYatP =  g.GetY()[p]+g.GetEYhigh()[p]
            maxXatP = g.GetX()[p]+g.GetEXhigh()[p]
            if varName != "pt" or maxXatP < maxPt*1.001:
                if (varName != "pt" and p < g.GetN()-1) or \
                   (varName == "pt" and p < frame.FindBin(maxPt*0.999)-1):
                    maxYatP /= 0.85-0.04*len(idSets)
                maxY = max(maxY, maxYatP)
        graphs[varName].append(g)
maxY = ceil(maxY*10)/10

objs = []
for varName in frames:
    c = TCanvas("c%s" % varName, varName, 500, 500)
    c.cd()

    leg = TLegend(0.2, 0.85-0.04*len(idSets), 0.7, 0.85)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    frame = frames[varName]

    ranges = [frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax(), 0, maxY]
    if varName == "pt": ranges[1] = maxPt
    frame.GetXaxis().SetRangeUser(ranges[0],ranges[1])
    frame.SetMinimum(ranges[2])
    frame.SetMaximum(ranges[3])
    frame.Draw()

    for g in graphs[varName]:
        leg.AddEntry(g, g.GetTitle(), "lp")
        g.Draw("P")
    leg.Draw()

    c.cd()
    label = TPaveText(ranges[0], ranges[3]*0.95, ranges[1], ranges[3]*1.1, "NB")
    label.SetFillStyle(0)
    label.SetBorderSize(0)
    label.SetTextFont(43)
    label.SetTextSize(16)
    label.SetTextAlign(11)
    label.SetMargin(0)
    label.AddText("CMS 2015 work in progress #sqrt{s} = 13 TeV, L=2.XY fb^{-1}")
    #label.Draw()

    c.Print("c_%s_%s_%s.png" % (submod, varName, dataType))

    objs.extend([c, leg, label])

