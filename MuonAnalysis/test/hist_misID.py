#!/usr/bin/env python
import os
from ROOT import *

modeSet = {
    "ks":{"massbin":(100,0.45, 0.55), "ptbins":[4, 5, 6, 8, 10, 15, 20, 30, 50, 200], "aetabins":[0, 0.9, 1.2, 1.6, 2.4],},
    "phi":{"massbin":(100,0.99,1.06), "ptbins":[4, 5, 6, 8, 10, 15, 20, 30, 50, 200], "aetabins":[0, 0.9, 1.2, 1.6, 2.4],},
    "lamb":{"massbin":(100,1.08,1.22), "ptbins":[4, 5, 6, 8, 10, 15, 20, 30, 50, 200], "aetabins":[0, 0.9, 1.2, 1.6, 2.4],},
}

precut = "1"
idSet = {
    "Tight":"isTight", "Medium":"isMedium", "Loose":"isLoose",
    "Soft":"isSoft", "HighPt":"isHighPt",
    "GLB":"isGLB", "TRK":"isTRK", "STA":"isSTA", "RPC":"isRPC",
    "GLBPT":"isGLBPT", "TMLastLoose":"isTMLastLoose", "TMLastTight":"isTMLastTight",
    "TM2DLoose":"isTM2DLoose", "TM2DTight":"isTM2DTight",
    "OneLoose":"isOneLoose", "OneTight":"isOneTight",
    "LastLowPtLoose":"isLastLowPtLoose", "LastLowPtTight":"isLastLowPtTight",
    "GMTkChi2Compat":"isGMTkChi2Compat", "GMStaChi2Compat":"isGMStaChi2Compat", "GMTkKinkTight":"isGMTkKinkTight",
    "TMLastAngLoose":"isTMLastAngLoose", "TMLastAngTight":"isTMLastAngTight", "TMOneAngLoose":"isTMOneAngLoose", "TMOneAngTight":"isTMOneAngTight",
}

def makedirs(d, path):
    for p in path.split('/'):
        dd = d.GetDirectory(p)
        if dd == None: dd = d.mkdir(p)
        d = dd
    return d

def project(dirName, mode, fName):
    print "@@ Processing", fName

    ptbins = modeSet[mode]["ptbins"]
    aetabins = modeSet[mode]["aetabins"]
    nbins, minMass, maxMass = modeSet[mode]["massbin"]
    binW = 1000*(maxMass-minMass)/nbins ## in MeV unit

    if fName.startswith("root://"): f = TNetXNGFile(fName)
    else: f = TFile(fName)
    if f == None or f.IsZombie(): return
    tree = f.Get("%s/tree" % mode)

    fout = TFile("%s/%s/%s" % (dirName, mode, os.path.basename(fName)), "RECREATE")
    fout.cd()

    eventList = TEventList("eventList");
    eventList.SetDirectory(0)
    for ptbin in range(len(ptbins)-1):
        minPt, maxPt = ptbins[ptbin], ptbins[ptbin+1]
        for leg in range(2):
            cutBin = "trk_pt[%d] >= %f && trk_pt[%d] < %f" % (leg, minPt, leg, maxPt)
            print "@@@@ Building event list", os.path.basename(fName), "ptbin", ptbin, "...",
            tree.Draw(">>eventList", "(%s) && (%s)" % (precut, cutBin))
            print "done"

            for idName in idSet:
                cutID = "mu_dR[%d]<0.01 && mu_%s[%d]" % (leg, idSet[idName], leg)

                outdir = makedirs(fout, "%s/%s/leg%d_ptbin%d" % (mode, idName, leg+1, ptbin))
                outdir.cd()

                cutStrObjPass = TObjString("cutPass")
                cutStrObjFail = TObjString("cutFail")
                cutStrObjPass.SetString("(%s) && (%s) &&  (%s)" % (precut, cutBin, cutID))
                cutStrObjFail.SetString("(%s) && (%s) && !(%s)" % (precut, cutBin, cutID))
                cutStrObjPass.Write()
                cutStrObjFail.Write()

                hPass = TH1D("hPass", "Passing;Mass (GeV);Candidates per %f MeV" % binW, nbins, minMass, maxMass)
                hFail = TH1D("hFail", "Failing;Mass (GeV);Candidates per %f MeV" % binW, nbins, minMass, maxMass)
                tree.Draw("vtx_mass>>hPass", " (%s)" % cutID, "goff")
                tree.Draw("vtx_mass>>hFail", "!(%s)" % cutID, "goff")
                hPass.Write()
                hFail.Write()

    for aetabin in range(len(aetabins)-1):
        minAeta, maxAeta = aetabins[aetabin], aetabins[aetabin+1]
        for leg in range(2):
            print "@@@@ Building event list", os.path.basename(fName), "aetabin", aetabin, "...",
            cutBin = "fabs(trk_eta[%d]) >= %f && fabs(trk_eta[%d]) < %f" % (leg, minAeta, leg, maxAeta)
            tree.Draw(">>eventList", "(%s) && (%s)" % (precut, cutBin))
            print "done"

            for idName in idSet:
                cutID = "mu_dR[%d]<0.01 && mu_%s[%d]" % (leg, idSet[idName], leg)

                outdir = makedirs(fout, "%s/%s/leg%d_aetabin%d" % (mode, idName, leg+1, aetabin))
                outdir.cd()

                cutStrObjPass = TObjString("cutPass")
                cutStrObjFail = TObjString("cutFail")
                cutStrObjPass.SetString("(%s) && (%s) &&  (%s)" % (precut, cutBin, cutID))
                cutStrObjFail.SetString("(%s) && (%s) && !(%s)" % (precut, cutBin, cutID))
                cutStrObjPass.Write()
                cutStrObjFail.Write()

                hPass = TH1D("hPass", "Passing;Mass (GeV);Candidates per %f MeV" % binW, nbins, minMass, maxMass)
                hFail = TH1D("hFail", "Failing;Mass (GeV);Candidates per %f MeV" % binW, nbins, minMass, maxMass)
                tree.Draw("vtx_mass>>hPass", " (%s)" % cutID, "goff")
                tree.Draw("vtx_mass>>hFail", "!(%s)" % cutID, "goff")
                hPass.Write()
                hFail.Write()

    fout.Write()
    print "@@ Finished", fName

if __name__ == '__main__':
    import sys

    eosBase = "root://eoscms//eos/cms/store/user/jhgoh/MuonMisID/20160412_1"
    #eosBase = "/afs/cern.ch/user/j/jhgoh/eos/cms/store/user/jhgoh/MuonMisID/20160412_1"

    if len(sys.argv) == 1:
        from multiprocessing import Pool
        ## To run on a multicore host
        filesRD = ['%s/JetHT_2015D/JetHT/crab_20160412_152954/160412_133018/0000/ntuple_%d.root' % (eosBase, i) for i in range(1,162)]
        filesMC = ['%s/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160412_153031/160412_133052/0000/ntuple_%d.root' % (eosBase, i) for i in range(1,629)]

        p = Pool(processes = 8)
        for mode in modeSet:
            if not os.path.exists("RD/%s" % mode): os.makedirs("RD/%s" % mode)
            if not os.path.exists("MC/%s" % mode): os.makedirs("MC/%s" % mode)

            for fName in filesRD: p.apply_async(project, ["RD", mode, fName])
            for fName in filesMC: p.apply_async(project, ["MC", mode, fName])
        p.close()
        p.join()
    else:
        ## To run on single core with one file
        dataType, mode, fileIndex = sys.argv[1:]
        fileIndex = int(fileIndex)

        if not os.path.exists("RD/%s" % mode): os.makedirs("RD/%s" % mode)
        if not os.path.exists("MC/%s" % mode): os.makedirs("MC/%s" % mode)

        if dataType == 'RD':
            f = '%s/JetHT_2015D/JetHT/crab_20160412_152954/160412_133018/0000/ntuple_%d.root' % (eosBase, fileIndex)
            project(dataType, mode, f)
        elif dataType == 'MC':
            f = '%s/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160412_153031/160412_133052/0000/ntuple_%d.root' % (eosBase, fileIndex)
            project(dataType, mode, f)

