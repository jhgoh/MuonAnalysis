from MuonAnalysis.MuonIdentification.tp_prod_template_cfg import *

process.source.fileNames = [
    '/store/mc/RunIIFall15DR76/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/B05BAC03-ACA6-E511-A7A3-02163E0169D9.root',
]
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

process.tpTree.isMC = cms.bool(False)#cms.bool(isMC)

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
  + process.rpcMuonIds
  * process.tpTree
)

