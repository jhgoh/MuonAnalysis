from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'prod_RD_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/AOD'
config.Data.splitting = "LumiBased"

config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
config.Data.unitsPerJob = 200

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'

dateTrial = "20170306_1"
config.General.workArea = "MuonMisID/JetHT_2015D"
config.Data.outLFNDirBase = '/store/user/jhgoh/MuonMisID/%s/JetHT_2015D' % dateTrial
config.General.requestName = "MuonMisID_JetHT_2015D_%s" % dateTrial

#config.Data.lumiMask = 'notFinishedLumis.json'
#config.Data.unitsPerJob = 1
