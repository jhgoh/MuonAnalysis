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
    setattr(process.tnp.Categories, wp, cms.vstring(wp, "dummy[pass=1,fail=0]"))
    for key in process.tnp.Efficiencies.parameters_().keys():
        if not hasattr(getattr(process.tnp.Efficiencies, key), 'EfficiencyCategoryAndState'): continue
        getattr(process.tnp.Efficiencies, key).EfficiencyCategoryAndState = [wp, 'pass']
    process.tnp.OutputFileName = "tp_fit_SingleMuon_Run2015D_%s.root" % wp

for i in range(1, 162):
    if i in (2,3,4,19,20,30,46,49,60,67,69,78,81,93): continue
    if i in (101,103,106,109,115,121,130, 137, 138, 146, 148, 152, 154): continue
    process.tnp.InputFileNames.append('root://eoscms//eos/cms/store/user/jhgoh/RPCMuonEff/20150415_1/SingleMuon_2015D/SingleMuon/crab_20160415_143600/160415_123642/0000/tnp_%d.root' % i)
