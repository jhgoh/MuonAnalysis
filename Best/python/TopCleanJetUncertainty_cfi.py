import FWCore.ParameterSet.Config as cms

cleanJets = cms.EDFilter("TopJetUncertainty",
    doFilter = cms.bool(False),
    debug = cms.untracked.bool(False),
    jet = cms.InputTag("loosePatJetsPF"),
    met = cms.InputTag("patMETsPF"),
    selection = cms.PSet(
        cut = cms.string(""),
        minPt = cms.double(30),
        maxEta = cms.double(2.5),
    ),
    cleaning = cms.PSet(
        overlapCands = cms.VInputTag(),
        overlapDeltaR = cms.double(0.5),
        #cleanMethod = cms.untracked.string("subtract"),
        cleanMethod = cms.string(""),
        #cleanMethod = cms.untracked.string(""),
    ),
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999),
)

