import FWCore.ParameterSet.Config as cms

from SKKU.RPCMuonAnalysis.rpcMuonIds_cfi import *
rpcMuonIds.src = "muons"

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

    idMaps = cms.VPSet(
        cms.PSet(name=cms.untracked.string("RPCLoose"), src=cms.InputTag("rpcMuonIds:Loose")),
        cms.PSet(name=cms.untracked.string("RPCTight"), src=cms.InputTag("rpcMuonIds:Tight")),
        cms.PSet(name=cms.untracked.string("TStLoose"), src=cms.InputTag("rpcMuonIds:TwoStationLoose")),
        cms.PSet(name=cms.untracked.string("TStLoose"), src=cms.InputTag("rpcMuonIds:TwoStationTight")),
        cms.PSet(name=cms.untracked.string("LStLoose"), src=cms.InputTag("rpcMuonIds:LastStationLoose")),
        cms.PSet(name=cms.untracked.string("LStLoose"), src=cms.InputTag("rpcMuonIds:LastStationTight")),
        cms.PSet(name=cms.untracked.string("SStLoose"), src=cms.InputTag("rpcMuonIds:SecondStationLoose")),
        cms.PSet(name=cms.untracked.string("SStLoose"), src=cms.InputTag("rpcMuonIds:SecondStationTight")),
    ),
)

phi = ks.clone(vtxType = cms.untracked.string("phi"), minTrkSigXY = cms.untracked.double(-999),)
lamb = ks.clone(vtxType = cms.untracked.string("lambda"))
D0 = ks.clone(vtxType = cms.untracked.string("D0"))
Dp = ks.clone(vtxType = cms.untracked.string("D+"))
Bp = ks.clone(vtxType = cms.untracked.string("B+"))
jpsi = ks.clone(vtxType = cms.untracked.string("jpsi"))

#misIDSeq = cms.Sequence(ks + phi + lamb + D0 + Dp + Bp + jpsi)
misIDSeq = cms.Sequence(rpcMuonIds * ks + phi + lamb)

