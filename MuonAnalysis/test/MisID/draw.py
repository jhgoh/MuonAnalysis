#!/usr/bin/env python

from ROOT import *
from SKKU.CommonTools.tdrStyle import *
tdrStyle.SetTitleSize(0.055, "XYZ")
tdrStyle.SetPadTopMargin(0.1)
tdrStyle.SetTitleXOffset(1.15)
tdrStyle.SetTitleYOffset(1.15)

#mode, submod, xtitle = "lamb", "proton", "Proton"
mode, submod, xtitle = "ks", "pion", "Pion"
#mode, submod, xtitle = "phi", "Kp", "K^{+}"
#mode, submod, xtitle = "phi", "Km", "K^{-}"
dataType = "RD"
varNames = {
    "pt":("%s transverse momentum p_{T} (GeV)" % xtitle),
    "abseta":("%s pseudorapidity |#eta|" % xtitle)
}
idNames = [
    ("RPC","RPC Loose ID"),
    ("TRK","Tracker muons (no arbitration)"),
    ("Loose","Loose muon ID"),
    ("Medium","Medium muon ID"),
    ("Tight","Tight muon ID"),
    ("Soft","Soft muon ID"),
]
colors = [kRed, kBlue, kGreen+1, kMagenta, kBlack, kAzure+1]

frames = {}
graphs = {}

f = TFile("%s/fit_%s.root" % (dataType, mode))
#for idName in [x.GetName() for x in f.GetDirectory(submod).GetListOfKeys()]:
for i, (idName, idTitle) in enumerate(idNames):
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
        g.SetMarkerColor(colors[i])
        g.SetLineColor(colors[i])
        for p in xrange(g.GetN()):
            g.SetPoint(p, g.GetX()[p], g.GetY()[p]*100)
            g.SetPointError(p, g.GetEXlow()[p], g.GetEXhigh()[p], g.GetEYlow()[p]*100, g.GetEYhigh()[p]*100)
        graphs[varName].append(g)

objs = []
for varName in frames:
    c = TCanvas("c%s" % varName, varName, 500, 500)
    c.cd()

    leg = TLegend(0.2, 0.85-0.04*len(idNames), 0.7, 0.85)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    frame = frames[varName]

    ranges = [frame.GetXaxis().GetXmin(), frame.GetXaxis().GetXmax(), 0, 4.0]
    if submod in ("Kp", "Km"): ranges[-1] = 7.0
    if varName == "pt": ranges[1] = 20
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
    label.AddText("CMS work in progress #sqrt{s} = 13 TeV, L=1.XY fb^{-1})")
    label.Draw()

    c.Print("c_%s_%s_%s.png" % (submod, varName, dataType))

    objs.extend([c, leg, label])

