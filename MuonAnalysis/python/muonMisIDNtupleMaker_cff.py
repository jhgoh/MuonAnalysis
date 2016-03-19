import FWCore.ParameterSet.Config as cms

ks = cms.EDAnalyzer("MuonMisIDNtupleMaker",
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

phi = ks.clone(vtxType = cms.string("phi"))
lamb = ks.clone(vtxType = cms.string("lambda"))

