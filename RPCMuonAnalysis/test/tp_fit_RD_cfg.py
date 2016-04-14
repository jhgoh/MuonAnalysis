import FWCore.ParameterSet.Config as cms
process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

from SKKU.RPCMuonAnalysis.tp_fit_template_cff import *
process.tnp_RPCMuon = fitTemplate.clone()
process.tnp_RPCMuon.InputFileNames = cms.vstring("RPCMuonEff/tp_prod_SingleMuon.root",)
process.tnp_RPCMuon.Efficiencies.RPC_voigtExpo.UnbinnedVariables = cms.vstring("mass")
process.tnp_RPCMuon.OutputFileName = cms.string("tp_fit_SingleMuon.root")

process.p = cms.Path(process.tnp_RPCMuon)
