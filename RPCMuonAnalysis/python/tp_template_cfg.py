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
    fileNames = cms.untracked.vstring(
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/00A3E567-75A8-E511-AD0D-0CC47A4D769E.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/06CC1B3A-FDA7-E511-B02B-00259073E388.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/0A9FEDA2-6DA8-E511-A451-002590596490.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/0AEF074D-EBA7-E511-B229-0002C94CDAF4.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/12998942-7BA8-E511-B1AA-003048FFCB84.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/145E4DB2-EFA7-E511-8E21-00266CF3DFE0.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/148E0F6C-EEA7-E511-A70E-0090FAA588B4.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/149A16F7-6DA8-E511-8A40-003048FFCC0A.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/18D542EB-FAA7-E511-A011-00259073E4E8.root',
#            '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/24537A2D-0BA8-E511-8D7C-20CF300E9ECF.root',
'/store/data/Run2015D/SingleMuon/RECO/PromptReco-v4/000/259/890/00000/00AE40B3-2C7C-E511-9A51-02163E014162.root',
#  489EAD68-277C-E511-9AC0-02163E014729.root  7EA5B246-2F7C-E511-893A-02163E011B73.root  D87A16C4-467C-E511-AAF5-02163E0145F8.root
#00E81978-227C-E511-9866-02163E0144E7.root  4AF8EC23-357C-E511-8F1F-02163E01455F.root  863E8BB6-2C7C-E511-BA33-02163E0139B6.root  DE7E0D35-2F7C-E511-9FB0-02163E012878.root
#04B77360-307C-E511-8429-02163E01458E.root  4E1E26D4-357C-E511-97BE-02163E011CD2.root  8A0C3C8B-617C-E511-BA71-02163E014162.root  E225E5EB-2C7C-E511-8FEB-02163E013757.root
#04D0A662-2D7C-E511-8BF3-02163E0133E5.root  527B50B4-427C-E511-8860-02163E014475.root  8E07D753-277C-E511-937F-02163E0145FE.root  E463C4DC-3F7C-E511-8910-02163E01186E.root
#06273F77-227C-E511-A418-02163E0141CE.root  56877E9D-2F7C-E511-AC10-02163E0142BC.root  9C3CA3D8-297C-E511-BD27-02163E0142BC.root  E802E7DA-3F7C-E511-80C5-02163E0145C4.root
#08930E39-287C-E511-B74E-02163E0145A7.root  5699ECB5-2D7C-E511-B5B1-02163E0141EA.root  A21665E5-347C-E511-B85B-02163E014367.root  E807034B-277C-E511-922C-02163E0139BE.root
#0EFD31D0-2F7C-E511-A90C-02163E01444B.root  58C79FCE-277C-E511-B8B1-02163E0141A2.root  AEEBC7E9-2D7C-E511-97AF-02163E0141A2.root  ECD6480E-287C-E511-A2DB-02163E01430F.root
#145FC319-427C-E511-9204-02163E0142BC.root  5A3AC284-2F7C-E511-9F52-02163E01418C.root  B0FDD32A-277C-E511-BB79-02163E014541.root  EEDDDF34-347C-E511-BDCD-02163E014162.root
#14EB3D63-3D7C-E511-8D19-02163E014358.root  64FA3E68-2D7C-E511-ABE0-02163E0146C7.root  BC5E4F74-227C-E511-AB1E-02163E011D7C.root  F0CFFE6F-227C-E511-BF4E-02163E0146C7.root
#243F13B1-297C-E511-A5A1-02163E0138F7.root  68E6312D-947C-E511-B6DB-02163E011AB1.root  C0226C02-377C-E511-9FB0-02163E01455F.root  F6F772A8-297C-E511-87F5-02163E01389F.root
#246C8FDC-2C7C-E511-B3B7-02163E011D7C.root  6ABB1490-227C-E511-B717-02163E013958.root  C6499958-277C-E511-BF8B-02163E014370.root  FA0C3CC6-297C-E511-9003-02163E013975.root
#28F56503-2A7C-E511-886F-02163E014553.root  6C9616FE-437C-E511-9EB7-02163E01438F.root  C66648C1-297C-E511-9E5B-02163E0134B1.root  FCF35905-2A7C-E511-93A7-02163E01446B.root
#323B47F0-297C-E511-ADF4-02163E014553.root  6CA22AF0-2C7C-E511-A682-02163E011D7C.root  C6AD3B60-2D7C-E511-A917-02163E0133E5.root  FE7DB0E6-2E7C-E511-A975-02163E01453A.root
#46951863-227C-E511-A93F-02163E0144A8.root  720A2342-517C-E511-AD8C-02163E0125D6.root  D2A02A3E-2A7C-E511-8EA5-02163E0144CE.root
#469BC7F9-4E7C-E511-B141-02163E0142BB.root  7E7241D6-8E7C-E511-B124-02163E0118C4.root  D8341A01-2E7C-E511-9F6A-02163E0142B3.root

    ),
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = '76X_dataRun2_v15'
#process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

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
    cut = cms.string("pt > 15 && "+MuonIDFlags.Tight2012.value()+
                     " && !triggerObjectMatchesByCollection('hltL3MuonCandidates').empty()"+
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

def finalizeTnP(process, isMC=True):
    process.tpTree.isMC = cms.bool(False)#cms.bool(isMC)
    if isMC:
        process.tagMuonsMCMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
            src     = cms.InputTag("tagMuons"), # RECO objects to match
            matched = cms.InputTag("goodGenMuons"),   # mc-truth particle collection
            mcPdgId     = cms.vint32(13),  # one or more PDG ID (13 = muon); absolute values (see below)
            checkCharge = cms.bool(False), # True = require RECO and MC objects to have the same charge
            mcStatus = cms.vint32(1),      # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
            maxDeltaR = cms.double(0.3),   # Minimum deltaR for the match
            maxDPtRel = cms.double(0.5),   # Minimum deltaPt/Pt for the match
            resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
            resolveByMatchQuality = cms.bool(True), # False = just match input in order; True = pick lowest deltaR pair first
        )
        process.probeMuonsMCMatch = process.tagMuonsMCMatch.clone(
            src = "probeMuons", maxDeltaR = 0.3, maxDPtRel = 1.0,
            resolveAmbiguities = False,  resolveByMatchQuality = False
        )

        process.tpTree.pairVariables.genWeight = cms.InputTag("genWeightInfo", "genWeight")
        process.tpTree.addRunLumiInfo = cms.bool(False)

        process.p = cms.Path(
            process.fastFilter
          + process.mergedMuons + process.patMuonsWithTriggerSequence
          * process.tagMuons * process.oneTag
          + process.probeMuons * process.tpPairs * process.onePair
          * process.nverticesModule + process.njets30Module
          + process.genWeightInfo
          * process.tpTree
        )

    else:
        process.tpTree.tagVariables.instLumi = cms.InputTag("addEventInfo", "instLumi")
        process.tpTree.addRunLumiInfo = cms.bool(True)

        process.p = cms.Path(
            process.fastFilter
          + process.mergedMuons + process.patMuonsWithTriggerSequence
          * process.tagMuons * process.oneTag
          + process.probeMuons * process.tpPairs * process.onePair
          * process.nverticesModule + process.njets30Module
          + process.addEventInfo
          * process.tpTree
        )

