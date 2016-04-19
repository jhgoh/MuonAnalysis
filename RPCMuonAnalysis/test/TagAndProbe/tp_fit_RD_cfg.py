import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from SKKU.RPCMuonAnalysis.tp_fit_template_cff import *
process.tnp_RPCMuon = fitTemplate.clone()
process.tnp_RPCMuon.InputFileNames = cms.vstring()
process.tnp_RPCMuon.Efficiencies.RPC_voigtExpo.UnbinnedVariables = cms.vstring("mass")
process.tnp_RPCMuon.OutputFileName = cms.string("tp_fit_SingleMuon_Run2015D.root")

process.p = cms.Path(process.tnp_RPCMuon)

import sys
if len(sys.argv) > 2:
    #wp can be RPC, or combinations of (RPC, RPCLSt, RPCSSt, RPCTSt) x (Loose, Tight)
    wp = sys.argv[2]
    for key in process.tnp_RPCMuon.Efficiencies.parameters_().keys():
        getattr(process.tnp_RPCMuon.Efficiencies, key).EfficiencyCategoryAndState = [wp, 'pass']
    process.tnp_RPCMuon.OutputFileName = "tp_fit_SingleMuon_Run2015D_%s.root" % wp

for i in range(1, 162):
    if i < 100: continue
    if i in (2,3,4,19,20,30,46,49,60,67,69,78,81,93): continue
    if i in (101,103,106,109,115,121,130, 137, 138, 146, 148, 152, 154): continue
    process.tnp_RPCMuon.InputFileNames.append('root://eoscms//eos/cms/store/user/jhgoh/RPCMuonEff/20150415_1/SingleMuon_2015D/SingleMuon/crab_20160415_143600/160415_123642/0000/tnp_%d.root' % i)
