import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
#    allowUnscheduled = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.goodVertexFilter = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 25 && position.Rho <= 2"),
    filter = cms.bool(True),
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)
process.fastFilter = cms.Sequence(process.goodVertexFilter + process.noScraping)

from RecoMuon.MuonIdentification.calomuons_cfi import calomuons;
process.mergedMuons = cms.EDProducer("CaloMuonMerger",
    mergeTracks = cms.bool(True),
    mergeCaloMuons = cms.bool(False), # AOD
    muons     = cms.InputTag("muons"),
    caloMuons = cms.InputTag("calomuons"),
    tracks    = cms.InputTag("generalTracks"),
    minCaloCompatibility = calomuons.minCaloCompatibility,
    ## Apply some minimal pt cut
    muonsCut     = cms.string("pt > 3 && track.isNonnull"),
    caloMuonsCut = cms.string("pt > 3"),
    tracksCut    = cms.string("pt > 3"),
)

process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
process.muonMatchHLTL2.maxDeltaR = 0.3 # Zoltan tuning - it was 0.5
process.muonMatchHLTL3.maxDeltaR = 0.1
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
changeRecoMuonInput(process, "mergedMuons")

from MuonAnalysis.TagAndProbe.common_variables_cff import *
process.load("MuonAnalysis.TagAndProbe.common_modules_cff")

process.tagMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("pt > 21 && "+MuonIDFlags.Tight2012.value()+
                     " && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"+
                     " && !triggerObjectMatchesByPath('HLT_IsoMu20_v*').empty()"+
                     " && pfIsolationR04().sumChargedHadronPt/pt < 0.2"),
)
process.oneTag  = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tagMuons"), minNumber = cms.uint32(1))

process.probeMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string("track.isNonnull"),  # no real cut now
)
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    #cut = cms.string('60 < mass < 140 && abs(daughter(0).vz - daughter(1).vz) < 4'),
    cut = cms.string('60 < mass && abs(daughter(0).vz - daughter(1).vz) < 4'),
    decay = cms.string('tagMuons@+ probeMuons@-')
)
process.onePair = cms.EDFilter("CandViewCountFilter", src = cms.InputTag("tpPairs"), minNumber = cms.uint32(1))

process.tpTree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("None"),
    variables = cms.PSet(
        abseta = cms.string('abs(eta)'),
        charge = cms.string('charge'),
        phi = cms.string('phi'),
        pt = cms.string('pt'),
        eta = cms.string('eta'),
    ),
    flags = cms.PSet(
        Glb = cms.string('isGlobalMuon'),
        GlbPT = cms.string('muonID("GlobalMuonPromptTight")'),
        Loose = cms.string('isLooseMuon()'),
        Medium = cms.string('isPFMuon && innerTrack.validFraction >= 0.8 && ( isGlobalMuon && globalTrack.normalizedChi2 < 3 && combinedQuality.chi2LocalPosition < 12 && combinedQuality.trkKink < 20 && segmentCompatibility >= 0.303 || segmentCompatibility >= 0.451 )'),
        PF = cms.string('isPFMuon()'),
        TM = cms.string('isTrackerMuon()'),
        TMA = cms.string('muonID("TrackerMuonArbitrated")'),
        TMLSAT = cms.string('muonID("TMLastStationAngTight")'),
        TMLST = cms.string('muonID("TMLastStationTight")'),
        TMOSL = cms.string('muonID("TMOneStationLoose")'),
        TMOST = cms.string('muonID("TMOneStationTight")'),
        TMOSTQual = cms.string('muonID("TMOneStationTight") && track.numberOfValidHits > 10 && track.normalizedChi2()<1.8 && track.hitPattern.pixelLayersWithMeasurement>1'),
        Tight2012 = cms.string('isPFMuon && numberOfMatchedStations > 1 && muonID("GlobalMuonPromptTight") && abs(dB) < 0.2 && track.hitPattern.trackerLayersWithMeasurement > 5 && track.hitPattern.numberOfValidPixelHits > 0'),
        RPC = cms.string('isRPCMuon()'),
        #RPCLoose = cms.string('muonID("RPCMuLoose")'),
        #RPCMedium = cms.string('muonID("RPCMuMedium")'),
        #RPCTight = cms.string('muonID("RPCMuTight")'),
    ),
    tagVariables = cms.PSet(
        nVertices   = cms.InputTag("nverticesModule"),
        abseta = cms.string('abs(eta)'),
        eta = cms.string('eta'),
        phi = cms.string('phi'),
        pt = cms.string('pt'),
    ),
    tagFlags = cms.PSet(),
    pairVariables = cms.PSet(
        deltaR = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)'),
        dz = cms.string('daughter(0).vz - daughter(1).vz'),
        nJets30 = cms.InputTag("njets30Module"),
        pt = cms.string('pt'),
        rapidity = cms.string('rapidity')
    ),
    pairFlags = cms.PSet(),
    tagMatches       = cms.InputTag("tagMuonsMCMatch"),
    probeMatches     = cms.InputTag("probeMuonsMCMatch"),
    motherPdgId      = cms.vint32(22, 23),
    makeMCUnbiasTree       = cms.bool(False),
    checkMotherInUnbiasEff = cms.bool(False),
    allProbes              = cms.InputTag("probeMuons"),
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tnp.root"))

