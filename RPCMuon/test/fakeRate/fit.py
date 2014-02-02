#!/usr/bin/env python

import sys, os
from ROOT import *
from urllib import urlretrieve
if not os.path.exists('rootlogon.C'):
    urlretrieve('https://raw.github.com/cms-top-kr/tools/master/rootlogon.C', 'rootlogon.C')
gROOT.ProcessLine(".x rootlogon.C")

if len(sys.argv) < 2:
    print "fit.py : calculate fakerate by fitting"
    print "Usage : python -i fit.py MODE"
    print "        read hist_MODE.root and writes fit_MODE.root and image_MODE"
    sys.exit()

mode = sys.argv[1]
inputFileName = "hist_%s.root" % mode
outputFileName = "fit_%s.root" % mode
imageDirName = "image_%s" % mode

objs = []

def fit(hA, hB, c = None):
    nA = hA.Integral()
    nB = hB.Integral()
    nTotal = nA+nB

    massMin = hA.GetXaxis().GetXmin()
    massMax = hA.GetXaxis().GetXmax()

    mass = RooRealVar("mass", "mass", massMin, massMax, "GeV/c^{2}")
    mass.setBinning(RooBinning(hA.GetNbinsX(), massMin, massMax))
    hDataA = RooDataHist("hDataA", "mass", RooArgList(mass), hA)
    hDataB = RooDataHist("hDataB", "mass", RooArgList(mass), hB)

    ws = RooWorkspace("ws")
    getattr(ws, 'import')(mass)
    getattr(ws, 'import')(hDataA)
    getattr(ws, 'import')(hDataB)

    ws.factory("m0[%f, %f, %f]" % ((massMax+massMin)/2, massMin, massMax))
    if 'Kshort' in mode: ws.factory("Voigtian::sigA(mass, m0, w0[1e-2, 5e-3, 2e-2], sigmaA[1e-2, 1e-3, 1e-1])") # Kshort
    elif 'Phi'  in mode: ws.factory("Voigtian::sigA(mass, m0, w0[5e-3, 1e-3, 2e-2], sigmaA[2e-3, 1e-3, 5e-3])") # Phi
    elif 'Jpsi' in mode: ws.factory("Voigtian::sigA(mass, m0, w0[1e-2, 1e-3, 2e-2], sigmaA[2e-2, 5e-3, 5e-2])") # Jpsi
    ws.factory("Voigtian::sigB(mass, m0, w0, sigmaA)")
    ws.factory("Chebychev::bkgA(mass, {p0A[0, -5, 5]})")#, p1A[0, -5, 5]})")
    ws.factory("Chebychev::bkgB(mass, {p0B[0, -5, 5]})")#, p1B[0, -5, 5]})")
    if 'Jpsi' in mode: ws.factory("ratio[0.999, 0, 1]")
    else: ws.factory("ratio[0.003, 0, 1]")
    ws.factory("EXPR::nSigA('nSig*ratio', nSig[%f, 0, %f], ratio)" % (0.5*nTotal, 1.1*nTotal))
    ws.factory("EXPR::nSigB('nSig*(1-ratio)', nSig, ratio)")
    ws.factory("SUM::pdfA(nSigA*sigA, nBkgA[%f, 0, %f]*bkgA)" % (0.5*nA, 1.1*nA))
    ws.factory("SUM::pdfB(nSigB*sigB, nBkgB[%f, 0, %f]*bkgB)" % (0.5*nB, 1.1*nB))
    ws.factory("weight[1, 0, 1e9]")
    ws.factory("index[A,B]")
    ws.exportToCint()

    ws.index = ws.cat('index')
    weight = ws.var('weight')
    ws.pdfA = ws.pdf('pdfA')
    ws.pdfB = ws.pdf('pdfB')

    simPdf = RooSimultaneous("simPdf", "simPdf", ws.index);
    simPdf.addPdf(ws.pdfA, "A");
    simPdf.addPdf(ws.pdfB, "B");

    dataA = RooDataSet("dataA", "mass", RooArgSet(mass, weight), RooFit.WeightVar("weight"))
    dataB = RooDataSet("dataB", "mass", RooArgSet(mass, weight), RooFit.WeightVar("weight"))
    for i in range(hDataA.numEntries()):
        dataA.add(hDataA.get(i), hDataA.weight())
        dataB.add(hDataB.get(i), hDataB.weight())
    dataSim = RooDataSet("dataSim", "mass", RooArgSet(mass, weight), RooFit.Index(ws.index),
                        RooFit.Import("A", dataA), RooFit.Import("B", dataB), RooFit.WeightVar(weight))

    #simPdf.fitTo(dataSim)
    #simPdf.fitTo(dataSim, RooFit.Extended())

    RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
    simNLL = simPdf.createNLL(dataSim, RooFit.Extended(True))
    scanner = RooMinimizer(simNLL)

    nll = RooProfileLL("simPdfNLL", "", simNLL, RooArgSet(ws.var("ratio")))
    scanner.minimize("Minuit2","Scan")

    nll.getVal()
    profMinuit = nll.minimizer()
    profMinuit.setProfile(True)
    profMinuit.setStrategy(2)
    profMinuit.setPrintLevel(1)
    profMinuit.migrad()
    profMinuit.migrad()
    profMinuit.hesse()
    profMinuit.minos(RooArgSet(ws.var("ratio")))
    result = profMinuit.save()

    RooMsgService.instance().setGlobalKillBelow(RooFit.ERROR)
    #result = simPdf.fitTo(dataSim, RooFit.Save(), RooFit.Extended(), RooFit.Minos())

    if c != None:
        frameA = mass.frame()
        dataSim.plotOn(frameA, RooFit.Cut("index==index::A"))
        proj = RooFit.ProjWData(RooArgSet(ws.index), dataSim)
        s = RooFit.Slice(ws.index, "A")
        simPdf.plotOn(frameA, s, proj, RooFit.LineColor(kGreen))
        simPdf.plotOn(frameA, s, proj, RooFit.LineColor(kGreen), RooFit.Components("bkgA"), RooFit.LineStyle(kDashed));

        frameB = mass.frame()
        s = RooFit.Slice(ws.index, "B")
        dataSim.plotOn(frameB, RooFit.Cut("index==index::B"))
        simPdf.plotOn(frameB, s, proj, RooFit.LineColor(kRed))
        simPdf.plotOn(frameB, s, proj, RooFit.LineColor(kRed), RooFit.Components("bkgB"), RooFit.LineStyle(kDashed))

        c.Divide(2,2)
        c.cd(1)
        frameA.Draw()
        c.cd(2)
        frameB.Draw()
        c.cd(3)
        frameNLL = ws.var("ratio").frame(RooFit.Range(0,2e-2))
        nll.plotOn(frameNLL, RooFit.Range(0,2e-2))
        frameNLL.Draw()

        c.cd(4)
        l = TPaveText(0,0,1,1)
        l.SetTextAlign(11)
        l.SetFillStyle(0)
        l.AddText("Fit results")
        pars = result.floatParsFinal()
        for i in xrange(pars.getSize()):
            par = pars[i]
            if par.hasAsymError():
                l.AddText(" %s = %f +%f - %f" % (par.GetName(), par.getVal(), par.getErrorHi(), par.getErrorLo()))
            else:
                l.AddText(" %s = %f +- %f" % (par.GetName(), par.getVal(), par.getError()))
        l.Draw()
        objs.append(l)

    #result.Print("v");

    #cA->Print(Form("cA_%s_eta_RPC_%s.png", category, region));
    #cB->Print(Form("cB_%s_eta_RPC_%s.png", category, region));

    ratio = result.floatParsFinal().find('ratio')
    return ratio;

if not os.path.isdir(imageDirName) :os.mkdir(imageDirName)
histFile = TFile(inputFileName)
fitFile  = TFile(outputFileName, "RECREATE")
for catName in [x.GetName() for x in histFile.GetListOfKeys()]:
    catDir = histFile.GetDirectory(catName)
    if catDir == None: continue
    outCatDir = fitFile.mkdir(catName)
    if not os.path.isdir("%s/%s" % (imageDirName, catName)): os.mkdir("%s/%s" % (imageDirName, catName))

    for varName in [x.GetName() for x in catDir.GetListOfKeys()]:
        varDir = catDir.GetDirectory(varName)
        if varDir == None: continue
        hFrame = varDir.Get("hFrame")
        if hFrame == None: continue
        outVarDir = outCatDir.mkdir(varName)
        outVarDir.cd()
        hFrame = hFrame.Clone()

        if not os.path.isdir("%s/%s/%s" % (imageDirName, catName, varName)): os.mkdir("%s/%s/%s" % (imageDirName, catName, varName))

        hFrame.SetMinimum(0)
        hFrame.SetMaximum(1)
        hFrame.GetYaxis().SetTitle("Fake rate (%)")
        c = TCanvas("c_%s_%s" % (catName, varName), "%s %s" % (catName, varName), 500, 500)
        grp = TGraphAsymmErrors()
        grp.SetName("ratio")
        grp.SetTitle(catName)
        for bin in range(hFrame.GetNbinsX()):
            binName = 'bin_%d' % bin
            binDir = varDir.GetDirectory(binName)
            if binDir == None: continue

            hM_pass = binDir.Get("hM_pass")
            hM_fail = binDir.Get("hM_fail")
            cFitCanvas = TCanvas("cFit_%s_%s_%s" % (catName, varName, binName), "c %s %s %s" % (catName, varName, binName), 600, 600)
            ratio = fit(hM_pass, hM_fail, cFitCanvas)
            cFitCanvas.Write()
            objs.append(cFitCanvas)
            cFitCanvas.Print("%s/%s/%s/%s.png" % (imageDirName, catName, varName, cFitCanvas.GetName()))
            cFitCanvas.Print("%s/%s/%s/%s.pdf" % (imageDirName, catName, varName, cFitCanvas.GetName()))

            x    = hFrame.GetBinCenter(bin+1)
            dx   = hFrame.GetBinWidth(bin+1)/2
            y    = 100*ratio.getVal()
            #dy   = abs(100*ratio.getError())
            dyLo  = abs(100*ratio.getErrorLo())
            dyHi  = abs(100*ratio.getErrorHi())
            grp.SetPoint(bin, x, y)
            grp.SetPointError(bin, dx, dx, dyLo, dyHi)
            if y > hFrame.GetMaximum(): hFrame.SetMaximum(y*1.1)

        c.cd()
        hFrame.Draw()
        grp.Draw("P")

        c.Write()
        hFrame.Write()
        grp.Write()
        objs.extend([c, hFrame, grp])

        c.Print("%s/%s/%s/%s.png" % (imageDirName, catName, varName, c.GetName()))
        c.Print("%s/%s/%s/%s.pdf" % (imageDirName, catName, varName, c.GetName()))

