#!/usr/bin/env python

from ROOT import *
from array import array

from MuonAnalysis.MuonIdentification.tdrStyle import *

gStyle.SetOptTitle(0)
gStyle.SetOptStat(0)

wps = [
    ("Loose", TFile("tp_fit_SingleMuon_Run2015D_Loose.root"), kRed),
    ("LooseWithRPC", TFile("tp_fit_SingleMuon_Run2015D_LooseWithRPC.root"), kBlue),
    ("LooseWithRPCNoTM", TFile("tp_fit_SingleMuon_Run2015D_LooseWithRPCNoTM.root"), kBlue),
#    ("Tight", TFile("tp_fit_SingleMuon_Run2015D_Tight2012.root"), kBlack),
#    ("RPCLoose", TFile("tp_fit_SingleMuon_Run2015D_RPCLoose.root"), kGreen+1),
#    ("RPCTight", TFile("tp_fit_SingleMuon_Run2015D_RPCTight.root"), kRed+2),
#    ("RPCTwoStLoose", TFile("tp_fit_SingleMuon_Run2015D_RPCTStLoose.root"), kBlue),
#    ("RPCLastStLoose", TFile("tp_fit_SingleMuon_Run2015D_RPCLStLoose.root"), kAzure+1),
#    ("RPCNotFirstStLoose", TFile("tp_fit_SingleMuon_Run2015D_RPCSStLoose.root"), kGreen+1),
]

leg = TLegend(0.5, 0.20, 0.9, 0.4)
leg.SetFillStyle(0)
leg.SetBorderSize(0)

grps = []
for name, f, color in wps:
    canvas = f.Get("tpTree/RPC_voigtExpo/fit_eff_plots/abseta_PLOT")
    srcgrp = canvas.FindObject("hxy_fit_eff")

    gROOT.cd()
    grp = TGraphAsymmErrors()
    for i in range(srcgrp.GetN()):
        grp.SetPoint(i, srcgrp.GetX()[i], srcgrp.GetY()[i])
        grp.SetPointError(i, srcgrp.GetEXlow()[i], srcgrp.GetEXhigh()[i], srcgrp.GetEYlow()[i], srcgrp.GetEYhigh()[i])

    grp.SetLineColor(color)
    grp.SetMarkerColor(color)
    leg.AddEntry(grp, name, "lp")
    grps.append(grp)
    grp = None

grp = grps[0]
nPoint = grp.GetN()
xbins = array('d', [grp.GetX()[i]-grp.GetEXlow()[i] for i in range(nPoint-1)]
                  +[grp.GetX()[nPoint-1]+grp.GetEXhigh()[nPoint-1]])

gROOT.cd()
c = TCanvas("c", "c", 500, 500)
frm = TH1F("hFrame", ";Muon pseudorapidity |#eta|;Efficiency", len(xbins)-1, xbins)
frm.SetMinimum(0.8)
frm.SetMaximum(1.01)
frm.GetYaxis().SetNdivisions(505)
frm.GetYaxis().SetTitleSize(0.05)
frm.GetXaxis().SetTitleSize(0.05)
frm.GetXaxis().SetTitleOffset(1.1)
frm.Draw()

for g in grps: g.Draw("P")
leg.Draw()


