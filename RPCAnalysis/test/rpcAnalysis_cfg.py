import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCA")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.autoCond_condDBv2 import autoCond
process.GlobalTag.globaltag = autoCond['run2_mc']

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source.fileNames = [
    '/store/mc/RunIIFall15DR76Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/premixPU50_deterministic_76X_mcRun2_asymptotic_v12_ext1-v2/00000/002236D0-08F6-E511-A030-0CC47A4D7646.root'
]

process.rpcRecHitAnalysis = cms.EDAnalyzer("RPCRecHitAnalysis",
    recHit = cms.InputTag("rpcRecHits"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.rpcRecHitAnalysis
)


