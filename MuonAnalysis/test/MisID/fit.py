#!/usr/bin/env python

import sys, os
from ROOT import *
from SKKU.CommonTools.tdrStyle import *

objs = []

def makedirs(d, path):
    for p in path.split('/'):
        if p == '': continue
        dd = d.GetDirectory(p)
        if dd == None: dd = d.mkdir(p)
        d = dd
    return d

def fit(hA, hB, c = None):
    nA = hA.Integral()
    nB = hB.Integral()
    nTotal = nA+nB

    massMin = hA.GetXaxis().GetXmin()
    massMax = hB.GetXaxis().GetXmax()

    ws = RooWorkspace("ws")

    mass = RooRealVar("mass", "mass", massMin, massMax, "GeV/c^{2}")
    mass.setBinning(RooBinning(hA.GetNbinsX(), massMin, massMax))
    hDataA = RooDataHist("hDataA", "mass", RooArgList(mass), hA)
    hDataB = RooDataHist("hDataB", "mass", RooArgList(mass), hB)

    getattr(ws, 'import')(mass, RooCmdArg())
    getattr(ws, 'import')(hDataA, RooCmdArg())
    getattr(ws, 'import')(hDataB, RooCmdArg())

    ws.factory("m0[%f,%f,%f]" % ((massMax+massMin)/2, massMin, massMax))
    if "ks" in mode: ws.factory("Voigtian::sigA(mass, m0, w0[1e-2, 5e-3, 2e-2], sigmaA[1e-2, 1e-3, 1e-1])")
    elif "phi" in mode: ws.factory("Voigtian::sigA(mass, m0, w0[5e-3, 1e-3, 2e-2], sigmaA[2e-3, 1e-3, 5e-3])")
    elif "lamb" in mode: ws.factory("Voigtian::sigA(mass, m0, w0[1e-2, 1e-3, 2e-2], sigmaA[2e-2, 5e-3, 5e-2])")
    ws.factory("Voigtian::sigB(mass, m0, w0, sigmaA)")
    ws.factory("Chebychev::bkgA(mass, {p0A[0, -5, 5], p1A[0, -5, 5]})")
    ws.factory("Chebychev::bkgB(mass, {p0B[0, -5, 5], p1B[0, -5, 5]})")
    ws.factory("ratio[0.003, 0, 1]")
    ws.factory("EXPR::nSigA('nSig*ratio', nSig[%f,0,%f], ratio)" % (0.5*nTotal, 1.1*nTotal))
    ws.factory("EXPR::nSigB('nSig*(1-ratio)', nSig, ratio)")
    ws.factory("SUM::pdfA(nSigA*sigA, nBkgA[%f,0,%f]*bkgA)" % (0.5*nA, 1.1*nA))
    ws.factory("SUM::pdfB(nSigB*sigB, nBkgB[%f,0,%f]*bkgB)" % (0.5*nB, 1.1*nB))
    ws.factory("weight[1,0,1e12]")
    ws.factory("index[A,B]")
#    ws.exportToCint()

    ws.index = ws.cat('index')
    weight = ws.var('weight')
    ws.pdfA, ws.pdfB = ws.pdf('pdfA'), ws.pdf('pdfB')

    simPdf = RooSimultaneous("simPdf", "simPdf", ws.index)
    simPdf.addPdf(ws.pdfA, "A")
    simPdf.addPdf(ws.pdfB, "B")

    dataA = RooDataSet("dataA", "mass", RooArgSet(mass, weight), RooFit.WeightVar("weight"))
    dataB = RooDataSet("dataB", "mass", RooArgSet(mass, weight), RooFit.WeightVar("weight"))
    for i in range(hDataA.numEntries()):
        dataA.add(hDataA.get(i), hDataA.weight())
        dataB.add(hDataB.get(i), hDataB.weight())
    dataSim = RooDataSet("dataSim", "mass", RooArgSet(mass, weight), RooFit.Index(ws.index),
                          RooFit.Import("A", dataA), RooFit.Import("B", dataB),
                          RooFit.WeightVar(weight))

    RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
    simNLL = simPdf.createNLL(dataSim, RooFit.Extended(True))
    scanner = RooMinimizer(simNLL)

    nll = RooProfileLL("simPdfNLL", "", simNLL, RooArgSet(ws.var("ratio")))
    scanner.minimize("Minuit2", "scan")

    nll.getVal()
    profMinuit = nll.minimizer()
    profMinuit.setProfile(True)
    profMinuit.setStrategy(2)
    profMinuit.setPrintLevel(1)
    profMinuit.migrad()
    profMinuit.minos(RooArgSet(ws.var("ratio")))
    result = profMinuit.save()

    RooMsgService.instance().setGlobalKillBelow(RooFit.ERROR)

    if c != None:
        frameA = mass.frame()
        dataSim.plotOn(frameA, RooFit.Cut("index==index::A"))
        proj = RooFit.ProjWData(RooArgSet(ws.index), dataSim)
        s = RooFit.Slice(ws.index, "A")
        simPdf.plotOn(frameA, s, proj, RooFit.LineColor(kGreen))
        simPdf.plotOn(frameA, s, proj, RooFit.LineColor(kGreen), RooFit.Components("bkgA"), RooFit.LineStyle(kDashed))

        frameB = mass.frame()
        dataSim.plotOn(frameB, RooFit.Cut("index==index::B"))
        proj = RooFit.ProjWData(RooArgSet(ws.index), dataSim)
        s = RooFit.Slice(ws.index, "B")
        simPdf.plotOn(frameB, s, proj, RooFit.LineColor(kRed))
        simPdf.plotOn(frameB, s, proj, RooFit.LineColor(kRed), RooFit.Components("bkgB"), RooFit.LineStyle(kDashed))


        c.Divide(2,2)
        c.cd(1)
        frameA.Draw()
        c.cd(2)
        frameB.Draw()
        c.cd(3)
        frameNLL = ws.var("ratio").frame(RooFit.Range(0, 2e-2))
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
                l.AddText(" %s = %f + %f - %f" % (par.GetName(), par.getVal(), par.getErrorHi(), par.getErrorLo()))
            else:
                l.AddText(" %s = %f +- %f" % (par.GetName(), par.getVal(), par.getError()))
        l.Draw()
        objs.append(l)

    ws = None
    ratio = result.floatParsFinal().find('ratio')
    return ratio

if __name__ == '__main__':
    gROOT.SetBatch(1)

    mode = sys.argv[1]
    fNameIn, fNameOut = sys.argv[2], sys.argv[3]

    fIn = TFile(fNameIn)
    fOut = TFile(fNameOut, "RECREATE")

    modeDir = fIn.GetDirectory(mode)
    if modeDir == None: os.exit(1)
    modeDirOut = fOut.mkdir(mode)

    for idName in set([x.GetName()[:-5] for x in modeDir.GetListOfKeys()]):
        idDir1 = modeDir.GetDirectory("%s_leg1" % idName)
        idDir2 = modeDir.GetDirectory("%s_leg2" % idName)
        for varName in [x.GetName() for x in idDir1.GetListOfKeys()]:
            varDir1 = idDir1.Get(varName)
            varDir2 = idDir2.Get(varName)
            hFrame = varDir1.Get("hFrame")

            if mode == 'lamb':
                varDirOut1 = makedirs(fOut, '/'.join(['proton', idName, varName]))
                varDirOut1.cd()
                hFrame.Clone().Write()
                gRatio1 = TGraphAsymmErrors()
                gRatio1.SetName("gRatio")
            elif mode == 'ks':
                varDirOut1 = makedirs(fOut, '/'.join(['pion', idName, varName]))
                varDirOut1.cd()
                gRatio1 = TGraphAsymmErrors()
                gRatio1.SetName("gRatio")
            elif mode == 'phi':
                varDirOut1 = makedirs(fOut, '/'.join(['Kp', idName, varName]))
                varDirOut1.cd()
                gRatio1 = TGraphAsymmErrors()
                gRatio1.SetName("gRatio")
                varDirOut2 = makedirs(fOut, '/'.join(['Km', idName, varName]))
                varDirOut2.cd()
                gRatio2 = TGraphAsymmErrors()
                gRatio2.SetName("gRatio")

            for b in range(hFrame.GetNbinsX()):
                bb = b+1
                hA1 = varDir1.Get("bin%d/hPass" % (bb))
                hB1 = varDir1.Get("bin%d/hFail" % (bb))
                hA2 = varDir2.Get("bin%d/hPass" % (bb))
                hB2 = varDir2.Get("bin%d/hFail" % (bb))

                x  = hFrame.GetXaxis().GetBinCenter(bb)
                ex = hFrame.GetXaxis().GetBinWidth(bb)/2

                if mode == 'lamb':
                    varDirOut1.cd()
                    c = TCanvas("c_bin%d" % bb, "c_bin%d" % bb, 500, 500)
                    res = fit(hA1, hB1, c)
                    c.Write()

                    y = res.getVal()
                    if res.hasAsymError(): ey1, ey2 = res.getErrorHi(), res.getErrorLo()
                    else:  ey1, ey2 = res.getError(), res.getError()

                    gRatio1.SetPoint(b, x, y)
                    gRatio1.SetPointError(b, ex, ex, eyLo, eyHi)

                elif mode == 'ks':
                    varDirOut1.cd()
                    c = TCanvas("c", "c", 500, 500)
                    hA1.Add(hA2)
                    hB1.Add(hB2)
                    res = fit(hA1, hB1, c)
                    c.Write()

                    y = res.getVal()
                    if res.hasAsymError(): eyHi, eyLo = res.getErrorHi(), res.getErrorLo()
                    else:  eyHi, eyLo = res.getError(), res.getError()

                    gRatio1.SetPoint(b, x, y)
                    gRatio1.SetPointError(b, ex, ex, eyLo, eyHi)

                elif mode == 'phi':
                    varDirOut1.cd()
                    c = TCanvas("c", "c", 500, 500)
                    res = fit(hA1, hB1, c)
                    c.Write()

                    y = res.getVal()
                    if res.hasAsymError(): eyHi, eyLo = res.getErrorHi(), res.getErrorLo()
                    else:  eyHi, eyLo = res.getError(), res.getError()

                    gRatio1.SetPoint(b, x, y)
                    gRatio1.SetPointError(b, ex, ex, eyLo, eyHi)

                    varDirOut2.cd()
                    c = TCanvas("c", "c", 500, 500)
                    res = fit(hA2, hB2, c)
                    c.Write()

                    y = res.getVal()
                    if res.hasAsymError(): eyHi, eyLo = res.getErrorHi(), res.getErrorLo()
                    else:  eyHi, eyLo = res.getError(), res.getError()

                    gRatio2.SetPoint(b, x, y)
                    gRatio2.SetPointError(b, ex, ex, eyLo, eyHi)
            gRatio1.Write()
            if gRatio2 in locals(): gRatio2.Write()
            break
        break


