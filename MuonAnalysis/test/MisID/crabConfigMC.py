from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True
config.General.workArea = "MuonMisID/TT_powheg"

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'prod_misID_MC_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.inputDataset = "/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/AODSIM"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 10

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'
#crab checkwrite --site=T3_KR_KISTI --lfn=/store/group/CAT/
config.Data.outLFNDirBase = '/store/user/jhgoh/MuonMisID/20160416_1/TT_powheg'

