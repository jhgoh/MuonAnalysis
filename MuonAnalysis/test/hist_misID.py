#!/usr/bin/env python
import sys, os
from ROOT import *
from multiprocessing import Pool

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
    "TMLastAngLoose":"isTMLastAngLoose", "TMLastAngTight":"TMLastAngTight", "TMOneAngLoose":"isTMOneAngLoose", "TMOneAngTight":"isTMOneAngTight",
}

def project(*args):
    dirName, fName = args[0]
    print "@@ Processing", fName
    if fName.startswith("root://"): f = TNetXNGFile(fName)
    else: f = TFile(fName)
  
    if f == None or f.IsZombie(): return

    fout = TFile("%s/%s" % (dirName, os.path.basename(fName)), "RECREATE")
    fout.cd()

    for mode in modeSet:
        ptbins = modeSet[mode]["ptbins"]
        aetabins = modeSet[mode]["aetabins"]
        nbins, minMass, maxMass = modeSet[mode]["massbin"]

        tree = f.Get("%s/tree" % mode)
        modeDir = fout.mkdir(mode)
        for idName in idSet:
            print "@@ Processing", mode, idName
            idDir = modeDir.mkdir(idName)
            cutID1 = "mu_dR[0]<0.01 && mu_%s[0]" % idSet[idName]
            cutID2 = "mu_dR[1]<0.01 && mu_%s[1]" % idSet[idName]

            for i in range(len(ptbins)-1):
                minPt, maxPt = ptbins[i], ptbins[i+1]
                for j in range(len(aetabins)-1):
                    binDir1 = idDir.mkdir("leg1_bin_pt%d_aeta%d" % (i, j))
                    binDir2 = idDir.mkdir("leg2_bin_pt%d_aeta%d" % (i, j))

                    minAeta, maxAeta = aetabins[j], aetabins[j+1]
                    cutBin1 = "trk_pt[0] >= %f && trk_pt[0] < %f && fabs(trk_eta[0]) >= %f && fabs(trk_eta[0]) < %f" % (minPt, maxPt, minAeta, maxAeta)
                    cutBin2 = "trk_pt[1] >= %f && trk_pt[1] < %f && fabs(trk_eta[1]) >= %f && fabs(trk_eta[1]) < %f" % (minPt, maxPt, minAeta, maxAeta)

                    binDir1.cd()
                    hPass1 = TH1D("hPass", "Passing;Mass (GeV);Candidates per %f MeV", nbins, minMass, maxMass)
                    hFail1 = TH1D("hFail", "Failing;Mass (GeV);Candidates per %f MeV", nbins, minMass, maxMass)
                    tree.Draw("vtx_mass>>hPass", "(%s) && (%s) &&  (%s)" % (precut, cutBin1, cutID1), "goff")
                    tree.Draw("vtx_mass>>hFail", "(%s) && (%s) && !(%s)" % (precut, cutBin1, cutID1), "goff")
                    hPass1.Write()
                    hFail1.Write()

                    binDir2.cd()
                    hPass2 = TH1D("hPass", "Passing;Mass (GeV);Candidates per %f MeV", nbins, minMass, maxMass)
                    hFail2 = TH1D("hFail", "Failing;Mass (GeV);Candidates per %f MeV", nbins, minMass, maxMass)
                    tree.Draw("vtx_mass>>hPass", "(%s) && (%s) &&  (%s)" % (precut, cutBin2, cutID2), "goff")
                    tree.Draw("vtx_mass>>hFail", "(%s) && (%s) && !(%s)" % (precut, cutBin2, cutID2), "goff")
                    hPass2.Write()
                    hFail2.Write()
    fout.Write()
    print "@@ Finished", fName

if __name__ == '__main__':
    if not os.path.exists("RD"): os.mkdir("RD")
    if not os.path.exists("MC"): os.mkdir("MC")

    eosBase = "root://eoscms//eos/cms/store/user/jhgoh/MuonMisID/20160412_1"
    filesRD = ['%s/JetHT_2015D/JetHT/crab_20160412_152954/160412_133018/0000/hist_%d.root' % (eosBase, i) for i in range(1,162)]
    filesMC = ['%s/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160412_153031/160412_133052/0000/hist_%d.root' % (eosBase, i) for i in range(1, 629)]

    p = Pool(processes = 8)
    p.map(project, zip(["RD"]*len(filesRD), filesRD))

    p.map(project, zip(["MC"]*len(filesMC), filesMC))

