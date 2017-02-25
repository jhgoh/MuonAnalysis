import FWCore.ParameterSet.Config as cms

rpcMuonIds = cms.EDProducer("RPCMuonIdMapProducer",
#    src = cms.InputTag("muons"),
    src = cms.InputTag("patMuons"),
)
