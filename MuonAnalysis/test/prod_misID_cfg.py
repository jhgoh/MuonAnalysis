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

process.ks = cms.EDAnalyzer("MuonMisIDNtupleMaker",
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    vertex = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    pfCandidates = cms.InputTag("packedCandidates"), ## For the MiniAOD
    #genParticles = cms.InputTag("packedGenParticles"), ## For the MiniAOD

    trkMinPt = cms.double(4.0),
    trkMaxEta = cms.double(2.4),
    trkChi2 = cms.double(5.),
    trkNHit = cms.int32(6),
    trkSignif = cms.double(-5),
    trkDCA = cms.double(1.),

    vtxType = cms.string("kshort"),
    vtxMinLxy = cms.double(.0),
    vtxMaxLxy = cms.double(4),
    vtxChi2 = cms.double(7.),
    vtxSignif = cms.double(.0),
)
process.phi = process.ks.clone(vtxType = cms.string("phi"))
process.lamb = process.ks.clone(vtxType = cms.string("lambda"))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("hist.root"),
)

process.p = cms.Path(
    process.ks + process.phi + process.lamb
)

