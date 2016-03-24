import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load('Configuration/StandardSequences/Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.autoCond_condDBv2 import autoCond
process.GlobalTag.globaltag = autoCond['run2_mc']
#process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/mc/RunIIFall15DR76/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/D0188731-42BC-E511-B9C0-02163E00B00F.root',
]

process.load("SKKU.MuonAnalysis.muonMisIDNtupleMaker_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.misIDSeq
)

