import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from SKKU.RPCMuonAnalysis.tp_fit_template_cff import *
process.tnp_RPCMuon = fitTemplate.clone()
process.tnp_RPCMuon.Efficiencies.RPC_voigtExpo.UnbinnedVariables = cms.vstring("mass")#, "weight")
process.tnp_RPCMuon.OutputFileName = cms.string("tp_fit_DYJets_MG.root")

process.p = cms.Path(process.tnp_RPCMuon)

import sys
if len(sys.argv) > 2:
    #wp can be RPC, or combinations of (RPC, RPCLSt, RPCSSt, RPCTSt) x (Loose, Tight)
    wp = sys.argv[2]
    for key in process.tnp_RPCMuon.Efficiencies.parameters_().keys():
        getattr(process.tnp_RPCMuon.Efficiencies, key).EfficiencyCategoryAndState = [wp, 'pass']
    process.tnp_RPCMuon.OutputFileName = "tp_fit_DYJets_MG_%s.root" % wp

for i in range(1, 25):
    if i in (15, 18): continue
    process.tnp_RPCMuon.InputFileNames.append('root://eoscms//eos/cms/store/user/jhgoh/RPCMuonEff/20150415_1/DYJets_MG/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_20160415_143449/160415_123541/0000/tnp_%d.root' %i)
