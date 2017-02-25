import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from SKKU.RPCMuonAnalysis.tp_fit_template_cff import *
process.tnp = fitTemplate.clone()
process.tnp.InputFileNames = cms.vstring()
process.tnp.Efficiencies.RPC_voigtExpo.UnbinnedVariables = cms.vstring("mass")
process.tnp.OutputFileName = cms.string("tp_fit_SingleMuon_Run2015D.root")

process.p = cms.Path(process.tnp)

import sys
if len(sys.argv) > 2:
    #wp can be RPC, or combinations of (RPC, RPCLSt, RPCSSt, RPCTSt) x (Loose, Tight)
    wp = sys.argv[2]
    cutwp = 'cut'+wp

    if hasattr(process.tnp.Categories, wp): cutwp = wp

    if not hasattr(process.tnp.Cuts, cutwp):
        setattr(process.tnp.Variables, wp, cms.vstring(wp, "-1", "2", "")) ## Register the flag variable, this must be defined in the tree
        setattr(process.tnp.Cuts, cutwp, cms.vstring(cutwp, wp, "0.5"))

    for key in process.tnp.Efficiencies.parameters_().keys():
        if not hasattr(getattr(process.tnp.Efficiencies, key), 'EfficiencyCategoryAndState'): continue
        if hasattr(process.tnp.Categories, cutwp):
            getattr(process.tnp.Efficiencies, key).EfficiencyCategoryAndState = [cutwp, 'pass']
        elif hasattr(process.tnp.Cuts, cutwp):
            getattr(process.tnp.Efficiencies, key).EfficiencyCategoryAndState = [cutwp, 'above']

    process.tnp.OutputFileName = "tp_fit_SingleMuon_Run2015D_%s.root" % wp

import os
home = os.environ["HOME"]
for i in range(1, 162):
    #process.tnp.InputFileNames.append('root://eoscms//eos/cms/store/user/jhgoh/RPCMuonEff/20160524_1/SingleMuon_2015D/SingleMuon/crab_20160524_161748/160524_141823/0000/tnp_%d.root' % i)
    process.tnp.InputFileNames.append(home+'/eos/cms/store/user/jhgoh/RPCMuonEff/20160524_1/SingleMuon_2015D/SingleMuon/crab_20160524_161748/160524_141823/0000/tnp_%d.root' % i)
