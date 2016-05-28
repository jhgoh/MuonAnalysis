from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True
config.General.workArea = "RPCMuonEff/DYJets_MG"

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'tp_prod_MC_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.inputDataset = "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/AODSIM"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 20

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
#crab checkwrite --site=T3_KR_KISTI --lfn=/store/group/CAT/
config.Data.outLFNDirBase = '/store/user/jhgoh/RPCMuonEff/20150524_1/DYJets_MG'

