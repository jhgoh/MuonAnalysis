import FWCore.ParameterSet.Config as cms

process = cms.Process("Analyze")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.GlobalTag.globaltag = ""
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012D/SingleMu/RECO/PromptReco-v1/000/203/777/24F0F702-110B-E211-99C9-001D09F248F8.root',
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RPCAnalysis.root'),
)

process.RPC = cms.EDAnalyzer("MuonRPCRecHitAnalyzer",
    muon = cms.InputTag("muons"),
)

process.p = cms.Path(
    process.RPC
)
