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
    '/store/mc/RunIIFall15DR76/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PU25nsPoisson50_76X_mcRun2_asymptotic_v12_ext1-v1/20000/0069F61C-CBF3-E511-929C-02163E01769E.root'
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


