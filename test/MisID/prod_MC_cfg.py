import FWCore.ParameterSet.Config as cms
process = cms.Process("CATeX")

process.load('Configuration/StandardSequences/Services_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.autoCond_condDBv2 import autoCond
#process.GlobalTag.globaltag = autoCond['run2_mc']
process.GlobalTag.globaltag = "76X_mcRun2_asymptotic_v12"

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    allowUnscheduled = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/mc/RunIIFall15DR76/TT_TuneCUETP8M1_13TeV-powheg-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/001194BD-31BB-E511-BFC2-0CC47A4C8E46.root',
    #'/store/mc/RunIIFall15DR76/JpsiToMuMu_OniaFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/08789E29-3DAC-E511-BD75-002590D0AFF4.root'
    #'/store/mc/RunIIFall15DR76/BuToJpsiK_BFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/72B114B2-E6A6-E511-8132-02163E013B6C.root',
    #'/store/mc/RunIIFall15DR76/InclusiveBtoJpsitoMuMu_JpsiPt8_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/80000/EC05CF1E-2AF5-E511-82F0-02163E0133DF.root',
]

process.load("MuonAnalysis.MuonIdentification.muonMisIDNtupleMaker_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

process.p = cms.Path(
    process.misIDSeq
)

