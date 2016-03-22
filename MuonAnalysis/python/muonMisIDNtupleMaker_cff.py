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

    trkMinPt = cms.untracked.double(4.0),
    trkMaxEta = cms.untracked.double(2.5),
    trkChi2 = cms.untracked.double(10.),
    trkNHit = cms.untracked.int32(7),
    trkSigXY = cms.untracked.double(2),
    trkSigZ = cms.untracked.double(-1),

    vtxDCA = cms.untracked.double(1.),

    vtxType = cms.untracked.string("kshort"),
    vtxMinLxy = cms.untracked.double(-4),
    vtxMaxLxy = cms.untracked.double(40),
    vtxChi2 = cms.untracked.double(7.),
    vtxSignif = cms.untracked.double(-5),
)

phi = ks.clone(vtxType = cms.untracked.string("phi"), trkSigXY = cms.untracked.double(0))
lamb = ks.clone(vtxType = cms.untracked.string("lambda"))

