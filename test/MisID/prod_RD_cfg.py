import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load('Configuration/StandardSequences/Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.autoCond_condDBv2 import autoCond
process.GlobalTag.globaltag = autoCond['run2_data']

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/data/Run2015D/JetHT/AOD/16Dec2015-v1/00000/0A2C6696-AEAF-E511-8551-0026189438EB.root'
]

process.load("SKKU.MuonAnalysis.muonMisIDNtupleMaker_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

process.p = cms.Path(
    process.misIDSeq
)

