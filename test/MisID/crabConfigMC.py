from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.transferLogs    = False
config.General.transferOutputs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'prod_MC_cfg.py'

config.section_("Data")
config.Data.publication  = False
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 25

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'

## Configurations that can change
dateTrial = "20170306_2"
config.Data.inputDataset = "/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/AODSIM"
import os
if 'DATASET' in os.environ:
    config.Data.inputDataset = os.environ['DATASET']
pd = config.Data.inputDataset.split('/')[1].replace('_TuneCUETP8M1_13TeV_pythia8', '')

config.General.workArea = "%s/%s" % (dateTrial, pd)
config.General.requestName = "MuonMisID_%s_%s" % (dateTrial, pd)
config.Data.outLFNDirBase = '/store/user/jhgoh/MuonMisID/%s/%s' % (dateTrial, pd)

