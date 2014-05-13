import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond["com10"]

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.throw = False
process.hltHighLevel.HLTPaths = ["HLT_IsoMu24_v*", "HLT_IsoMu24_eta2p1_v*",]

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.suppressError = cms.untracked.vstring('patTriggerFull')
process.MessageLogger.suppressWarning = cms.untracked.vstring('patTriggerFull')
process.MessageLogger.suppressInfo = cms.untracked.vstring('patTriggerFull')
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012D/SingleMu/RECO/PromptReco-v1/000/203/777/24F0F702-110B-E211-99C9-001D09F248F8.root',
    ),
    inputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_generalTracks_*_*",
        "keep *_siPixel*_*_*", "keep *_siStrip*_*_*",
        "keep *_dt*_*_*", "keep *_csc*_*_*", "keep *_rpc*_*_*",
        "keep *EcalRecHit*_*_*_*", "keep *H*RecHit*_h*reco_*_*",
        "keep *CaloTower*_*_*_*", "keep *_*CaloJets_*_*",
        "keep *_offlineBeamSpot_*_*",
        "keep *_offlinePrimaryVertices_*_*",
        "keep *_TriggerResults_*_*",
        "keep *_hltTriggerSummaryAOD_*_*",
        "keep *_l1extra*_*_*",
    ),
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
)

## Muon rereconstruction
if not hasattr(process, "muidRPCMuLoose"):
    process.muidRPCMuLoose = cms.EDProducer("MuonSelectionTypeValueMapProducer",
        inputMuonCollection = cms.InputTag("muons1stStep"),
        selectionType = cms.string("RPCMuLoose"),
    )
    process.muonSelectionTypeSequence += process.muidRPCMuLoose

## Disable PFMuon for temporary solution not to crash in 53X
process.muons.FillPFIsolation = False
process.muons.FillPFMomentumAndAssociation = False
process.muons.PFCandidates = "particleFlow"

process.muonRereco = cms.Sequence(
    process.muonrecoComplete
  * process.muoncosmicreco * process.regionalCosmicTracksSeq
  * process.muoncosmichighlevelreco #* process.muonshighlevelreco
  * process.muons
)

## Define common filters
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

# Good vertex requirement
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True),
)

process.commonFilters = cms.Sequence(
    process.noscraping + process.goodOfflinePrimaryVertices
  + process.hltHighLevel
)

## Setup Tag and probes
process.selectedMuons = cms.EDFilter("MuonSelector",
    src = cms.InputTag("muons"),
    cut = cms.string("pt > 20 && abs(eta) < 1.8"),
)

process.promptMuons = cms.EDFilter("PromptMuonSelector",
    src = cms.InputTag("selectedMuons"),
    maxDxy = cms.untracked.double(0.2),
    maxDz = cms.untracked.double(0.5),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
)

process.selectLooseMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(
        " (isGlobalMuon || isTrackerMuon)"# && isPFMuon"
    ),
)

process.selectMediumMuons = process.selectLooseMuons.clone(
    cut = cms.string(
        " isGlobalMuon"
        " && globalTrack().normalizedChi2()<10.0"
        " && innerTrack().hitPattern().numberOfValidPixelHits() > 0"
        " && track().hitPattern().trackerLayersWithMeasurement() > 5"
    ),
)

process.selectTightMuons = process.selectLooseMuons.clone(
    cut = cms.string(
        #" isGlobalMuon && isPFMuon"
        " isGlobalMuon"
        #" && isolationR03().sumPt<0.05*pt"
        " && innerTrack().hitPattern().numberOfValidPixelHits() > 0"
        " && track().hitPattern().trackerLayersWithMeasurement() > 5"
        " && globalTrack().normalizedChi2()<10.0"
        " && globalTrack().hitPattern().numberOfValidMuonHits()>0"
        " && numberOfMatchedStations>1"
    ),
)

process.selectLooseRPCMuons  = process.selectLooseMuons.clone(cut = cms.string("isRPCMuon"))
process.selectMediumRPCMuons = process.selectLooseMuons.clone(cut = cms.string("isRPCMuon && numberOfMatchedStations('RPCHitAndTrackArbitration') >= 2"))
process.selectTightRPCMuons  = process.selectLooseMuons.clone(cut = cms.string("isRPCMuon && numberOfMatchedStations('RPCHitAndTrackArbitration') >= 2 && numberOfMatchedRPCLayers >= 3"))

## Trigger matching with PAT
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "promptMuons")
#useExtendedL1Match(process)
#addHLTL1Passthrough(process)

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(
        #" isGlobalMuon && isPFMuon"
        " isGlobalMuon"
        #" && isolationR03().sumPt<0.05*pt"
        " && innerTrack().hitPattern().numberOfValidPixelHits() > 0"
        " && track().hitPattern().trackerLayersWithMeasurement() > 5"
        " && globalTrack().normalizedChi2()<10.0"
        " && globalTrack().hitPattern().numberOfValidMuonHits()>0"
        " && numberOfMatchedStations>1"
        " && (!triggerObjectMatchesByPath('HLT_IsoMu24_v*').empty()"
        "  || !triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*').empty() )"
    ),
    filter = cms.bool(True),
)

process.tagMuonFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("tagMuons"),
    minNumber = cms.uint32(1),
)

process.goodTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string("pt>20 && abs(eta)<1.8 && numberOfValidHits >= 10"),
)

process.trackCands = cms.EDProducer("ConcreteChargedCandidateProducer",
    src  = cms.InputTag("goodTracks"),
    particleType = cms.string("mu+"),
)

process.promptTrackCands = cms.EDFilter("PromptTrackCandSelector",
    src = cms.InputTag("trackCands"),
    maxDxy = cms.untracked.double(0.2),
    maxDz = cms.untracked.double(0.5),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
)

process.trackProbes = cms.EDFilter("CandViewRefSelector",
    src = cms.InputTag("promptTrackCands"),
    cut = cms.string("abs(eta)<1.8"),
)

## Matchings
process.matchLooseMuons = cms.EDProducer("MatcherUsingTracks",
    src     = cms.InputTag("promptTrackCands"),
    matched = cms.InputTag("selectLooseMuons"),
    algorithm = cms.string("byDirectComparison"),
    srcTrack = cms.string("tracker"),
    srcState = cms.string("atVertex"),
    matchedTrack = cms.string("tracker"),
    matchedState = cms.string("atVertex"),
    maxDeltaR        = cms.double(0.01),
    maxDeltaLocalPos = cms.double(0.01),
    maxDeltaPtRel    = cms.double(0.01),
    sortBy           = cms.string("deltaR"),
)

## Passing probe candidates
process.passLooseMuons = cms.EDProducer("MatchedCandidateSelector",
    src = cms.InputTag("trackProbes"),
    match = cms.InputTag("matchLooseMuons"),
)

process.muonSelectionSequence = cms.Sequence(
    process.selectedMuons * process.promptMuons
  * process.patMuonsWithTriggerSequence * process.tagMuons * process.tagMuonFilter
  + process.goodTracks * process.trackCands * process.promptTrackCands * process.trackProbes
  + process.selectLooseMuons * process.matchLooseMuons * process.passLooseMuons
)
for muonId in ["Medium", "Tight", "LooseRPC", "MediumRPC", "TightRPC"]:
    setattr(process, "match%sMuons" % muonId, process.matchLooseMuons.clone(
        matched = cms.InputTag("select%sMuons" % muonId)
    ))

    setattr(process, "pass%sMuons" % muonId, process.passLooseMuons.clone(
        match = cms.InputTag("match%sMuons" % muonId)
    ))

    process.muonSelectionSequence += getattr(process, "select%sMuons" % muonId)
    process.muonSelectionSequence *= getattr(process, "match%sMuons" % muonId)
    process.muonSelectionSequence *= getattr(process, "pass%sMuons" % muonId)

## Build Tag-Probe pair and tree
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ trackProbes@-"),
    cut   = cms.string("60 < mass < 120"),
)

process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    variables = cms.PSet(
        eta = cms.string("eta"),
        phi = cms.string("phi"),
        pt  = cms.string("pt"),
        abseta    = cms.string("abs(eta)"),
    ),
    flags = cms.PSet(
        probe = cms.InputTag("trackProbes"),
        looseMuons     = cms.InputTag("passLooseMuons"   ),
        mediumMuons    = cms.InputTag("passMediumMuons"  ),
        tightMuons     = cms.InputTag("passTightMuons"   ),
        looseRPCMuons  = cms.InputTag("passLooseRPCMuons"),
        mediumRPCMuons = cms.InputTag("passMediumRPCMuons"),
        tightRPCMuons  = cms.InputTag("passTightRPCMuons"),
    ),
    isMC = cms.bool(False),
)
    
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("tnpTree.root"),
)

process.p = cms.Path(
    process.commonFilters
  + process.muonRereco
  * process.muonSelectionSequence
  * process.tpPairs * process.muonEffs
)

