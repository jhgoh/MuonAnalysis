import FWCore.ParameterSet.Config as cms

ks = cms.EDAnalyzer("MuonMisIDNtupleMaker",
    muons = cms.InputTag("muons"),
    tracks = cms.InputTag("generalTracks"),
    vertex = cms.InputTag("offlinePrimaryVertices"),
    genParticles = cms.InputTag("genParticles"),
    pfCandidates = cms.InputTag("packedCandidates"), ## For the MiniAOD
    #genParticles = cms.InputTag("packedGenParticles"), ## For the MiniAOD

    applyGenFilter = cms.untracked.bool(False),
    useBeamSpot = cms.untracked.bool(True),

    trkMinPt = cms.untracked.double(4.0),
    trkMaxEta = cms.untracked.double(2.5),
    trkChi2 = cms.untracked.double(5.),
    trkNHit = cms.untracked.int32(6),
    trkSignif = cms.untracked.double(-5),
    trkDCA = cms.untracked.double(1.),

    vtxType = cms.untracked.string("kshort"),
    vtxMinLxy = cms.untracked.double(-4),
    vtxMaxLxy = cms.untracked.double(40),
    vtxChi2 = cms.untracked.double(7.),
    vtxSignif = cms.untracked.double(-5),
)

phi = ks.clone(vtxType = cms.untracked.string("phi"))
lamb = ks.clone(vtxType = cms.untracked.string("lambda"))

