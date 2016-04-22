#!/usr/bin/env python
import os
from ROOT import *
from array import array

modeSet = {
    "ks":{"massbin":(100,0.45, 0.55), "vars":{
        "pt":{'expr':"trk_pt[%d]",'title':"p_{T} (GeV)",'bins':[4, 5, 6, 8, 10, 15, 20, 30, 50, 200]},
        "abseta":{'expr':'fabs(trk_eta[%d])','title':'|#eta|','bins':[0, 0.9, 1.2, 1.6, 2.4]},
    }},
    "phi":{"massbin":(100,0.99,1.06), "vars":{
        "pt":{'expr':"trk_pt[%d]",'title':"p_{T} (GeV)",'bins':[4, 5, 6, 8, 10, 15, 20, 30, 50, 200]},
        "abseta":{'expr':'fabs(trk_eta[%d])','title':'|#eta|','bins':[0, 0.9, 1.2, 1.6, 2.4]},
    }},
    "lamb":{"massbin":(100,1.08,1.22), "vars":{
        "pt":{'expr':"trk_pt[%d]",'title':"p_{T} (GeV)",'bins':[4, 5, 6, 8, 10, 15, 20, 30, 50, 200]},
        "abseta":{'expr':'fabs(trk_eta[%d])','title':'|#eta|','bins':[0, 0.9, 1.2, 1.6, 2.4]},
    }},
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

def project(dirName, mode, fName):
    print "@@ Processing", fName

    nbins, minMass, maxMass = modeSet[mode]["massbin"]
    binW = 1000*(maxMass-minMass)/nbins ## in MeV unit
    xyTitles = "Mass (GeV);Candidates per %f MeV" % binW

    if fName.startswith("root://"): f = TNetXNGFile(fName)
    else: f = TFile(fName)
    if f == None or f.IsZombie(): return
    tree = f.Get("%s/tree" % mode)

    fout = TFile("%s/%s/%s" % (dirName, mode, os.path.basename(fName)), "RECREATE")
    modedir = fout.mkdir(mode)
    modedir.cd()

    varSet = modeSet[mode]["vars"]

    for idName in idSet:
        for leg in range(2):
            cutID = "mu_dR[%d]<0.01 && mu_%s[%d]" % (leg, idSet[idName], leg)
            idDir = modedir.mkdir("%s_leg%d" % (idName, leg+1))
            for varName in varSet:
                print "%s/%s/%s" % (idName, leg, varName)
                varDir = idDir.mkdir(varName)

                bins = varSet[varName]['bins']
                title = varSet[varName]['title']
                expr = varSet[varName]['expr']

                varDir.cd()
                hFrame = TH1D("hFrame", "%s;%s" % (varName, title), len(bins)-1, array('d', bins))
                hFrame.Write()

                for b in range(len(bins)-1):
                    minX, maxX = bins[b], bins[b+1]
                    cutBin = "%s >= %f && %s < %f" % ((expr%leg), minX, (expr%leg), maxX)
                    cutPass = "(%s) && (%s) &&  (%s)" % (precut, cutBin, cutID)
                    cutFail = "(%s) && (%s) && !(%s)" % (precut, cutBin, cutID)

                    binDir = varDir.mkdir("bin%d" % (b+1))
                    binDir.cd()

                    hPass = TH1D("hPass", "Passing=%s;%s" % (cutPass, xyTitles), nbins, minMass, maxMass)
                    hFail = TH1D("hFail", "Failing=%s;%s" % (cutFail, xyTitles), nbins, minMass, maxMass)
                    #tree.Draw("vtx_mass>>hPass", cutPass, "goff")
                    #tree.Draw("vtx_mass>>hFail", cutFail, "goff")
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
        filesRD = ['%s/JetHT_2015D/JetHT/crab_20160412_152954/160412_133018/0000/ntuple_%d.root' % (eosBase, i) for i in range(1,2)]#162)]
        filesMC = ['%s/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160412_153031/160412_133052/0000/ntuple_%d.root' % (eosBase, i) for i in range(1,2)]#629)]

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

