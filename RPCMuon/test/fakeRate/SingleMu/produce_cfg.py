import FWCore.ParameterSet.Config as cms
import sys, os

process = cms.Process("FakeRate")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond["com10"]
#process.GlobalTag.globaltag = autoCond["startup"]

process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.throw = False
process.hltHighLevel.HLTPaths = ["HLT_IsoMu24_v*", "HLT_IsoMu24_eta2p1_v*"]

process.MessageLogger.suppressError = cms.untracked.vstring('patTriggerFull',)
process.MessageLogger.suppressWarning = cms.untracked.vstring('patTriggerFull',)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/data/Run2012D/SingleMu/RECO/PromptReco-v1/000/203/853/F814E722-C70B-E211-8C51-BCAEC5329727.root",
#        "/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/GEN-SIM-RECO/v2/00000/52149A23-1E2A-E211-8BD3-003048678B14.root",
#        "/store/relval/CMSSW_5_3_6-START53_V14/RelValProdTTbar/GEN-SIM-RECO/v2/00000/B86B2DE8-122A-E211-AD41-003048678B84.root",
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
if not hasattr(process, 'muidRPCMuLoose'):
    process.muidRPCMuLoose = cms.EDProducer("MuonSelectionTypeValueMapProducer",
        inputMuonCollection = cms.InputTag("muons1stStep"),
        selectionType = cms.string('RPCMuLoose'),
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

## Trigger matching with PAT
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string(
        " (!triggerObjectMatchesByPath('HLT_IsoMu24_v*').empty()"
        "  || !triggerObjectMatchesByPath('HLT_IsoMu24_eta2p1_v*').empty() )"
    ),
    filter = cms.bool(True),
)

process.muonSelectionSeq = cms.Sequence(
    process.patMuonsWithTriggerSequence * process.tagMuons
)

## User analysis
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("fakeTree.root"),
)

process.load("SKKU.RPCMuon.VertexCandProducer_cfi")
process.load("SKKU.RPCMuon.fakeMuonAnalyzer_cfi")

process.fakeKshort.muon = "patMuonsWithTrigger"

process.commonFilters = cms.Sequence(
    process.noscraping + process.goodOfflinePrimaryVertices
  + process.hltHighLevel
)

process.pKs = cms.Path(
    process.commonFilters
  + process.kshortVertex
  + process.muonRereco * process.muonSelectionSeq
  * process.fakeKshort
)

process.pLambda = cms.Path(
    process.commonFilters
  + process.lambdaVertex
  + process.muonRereco * process.muonSelectionSeq
  #* process.fakeLambda2
  * process.fakeLambda
)

process.pPhi = cms.Path(
    process.commonFilters
  + process.phiVertex
  + process.muonRereco * process.muonSelectionSeq
  * process.fakePhi
)

process.pJpsi = cms.Path(
    process.commonFilters
  + process.jpsiVertex
  + process.muonRereco * process.muonSelectionSeq
  * process.fakeJpsi
)

