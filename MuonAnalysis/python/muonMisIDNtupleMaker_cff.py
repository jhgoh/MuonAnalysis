import FWCore.ParameterSet.Config as cms

ks = cms.EDAnalyzer("MuonMisIDNtupleMaker",
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    pfCandidates = cms.InputTag("packedCandidates"), ## For the MiniAOD
    vertex = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    #vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    #genParticles = cms.InputTag("packedGenParticles"), ## For the MiniAOD

    applyGenFilter = cms.untracked.bool(False),
    useBeamSpot = cms.untracked.bool(True),

    minTrkPt = cms.untracked.double(4.0),
    maxTrkEta = cms.untracked.double(2.5),
    maxTrkChi2 = cms.untracked.double(10.),
    minTrkNHit = cms.untracked.int32(7),
    minTrkSigXY = cms.untracked.double(2),
    minTrkSigZ = cms.untracked.double(-1),

    maxVtxDCA = cms.untracked.double(1.),

    vtxType = cms.untracked.string("kshort"),
    minVtxLxy = cms.untracked.double(-4),
    maxVtxLxy = cms.untracked.double(40),
    maxVtxChi2 = cms.untracked.double(7.),
    minVtxSignif = cms.untracked.double(-5),
)

phi = ks.clone(vtxType = cms.untracked.string("phi"), minTrkSigXY = cms.untracked.double(-999),)
lamb = ks.clone(vtxType = cms.untracked.string("lambda"))

