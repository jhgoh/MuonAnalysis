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
config.Data.splitting = "LumiBased"

config.Data.unitsPerJob = 200

config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite = 'T2_CH_CERN'

## Configurations that can change
dateTrial = "20170306_2"
config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/AOD'
import os
if 'DATASET' in os.environ:
    config.Data.inputDataset = os.environ['DATASET']
pd, sd = config.Data.inputDataset.split('/')[1:3]
lumiMaskBase = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/'
if '2015' in sd:
    config.Data.lumiMask = lumiMaskBase+'Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
elif '2016' in sd:
     config.Data.lumiMask = lumiMaskBase+'Collisions16/13TeV/Final/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'

config.General.workArea = "%s/%s_%s" % (dateTrial, pd, sd)
config.General.requestName = "MuonMisID_%s_%s_%s" % (dateTrial, pd, sd)
config.Data.outLFNDirBase = '/store/user/jhgoh/MuonMisID/%s/%s_%s' % (dateTrial, pd, sd)

#config.Data.lumiMask = 'notFinishedLumis.json'
#config.Data.unitsPerJob = 1
