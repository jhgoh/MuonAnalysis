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
#process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9')

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

process.ks = cms.EDAnalyzer("MuonMisIDNtupleMaker",
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    vertex = cms.InputTag("offlinePrimaryVertices"),
    #genParticles = cms.InputTag("genParticles"),
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

