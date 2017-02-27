from MuonAnalysis.MuonIdentification.tp_prod_template_cfg import *

process.source.fileNames = [
#   '/store/data/Run2015D/SingleMuon/AOD/16Dec2015-v1/10000/00A3E567-75A8-E511-AD0D-0CC47A4D769E.root',
    '/store/data/Run2015D/SingleMuon/RECO/PromptReco-v4/000/259/890/00000/00AE40B3-2C7C-E511-9A51-02163E014162.root',
]
process.GlobalTag.globaltag = '76X_dataRun2_v15'

process.tpTree.isMC = cms.bool(False)#cms.bool(isMC)

process.tpTree.tagVariables.instLumi = cms.InputTag("addEventInfo", "instLumi")
process.tpTree.addRunLumiInfo = cms.bool(True)

process.p = cms.Path(
    process.fastFilter
  + process.mergedMuons + process.patMuonsWithTriggerSequence
  * process.tagMuons * process.oneTag
  + process.probeMuons * process.tpPairs * process.onePair
  * process.nverticesModule + process.njets30Module
  + process.addEventInfo
  + process.rpcMuonIds
  * process.tpTree
)

