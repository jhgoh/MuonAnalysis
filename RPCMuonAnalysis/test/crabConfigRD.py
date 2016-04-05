from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True
config.General.workArea = "RPCMuonEff/SingleMuon_2015D"

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'prod_tp_RD_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.inputDataset = '/SingleMuon/Run2015D-16Dec2015-v1/AOD'
config.Data.splitting = "LumiBased"
config.Data.unitsPerJob = 200
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
#crab checkwrite --site=T3_KR_KISTI --lfn=/store/group/CAT/
config.Data.outLFNDirBase = '/store/user/jhgoh/RPCMuonEff/20150404_1/SingleMuon_2015D'

