from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = '1_NANOREG_2018'
config.General.workArea = '/afs/cern.ch/work/n/nchernya/crab_projects/'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.numCores = 4
config.JobType.maxMemoryMB = 6000
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'myNanoProdMc2018_NANO.py'

config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'

config.Data.inputDBS = 'global'
#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 50000
#config.Data.totalUnits = 50000*2
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.totalUnits = 1
config.Data.outLFNDirBase = '/store/user/nchernya/NanoAOD/NanoReg18/' 
config.Data.publication = True
config.Data.outputDatasetTag = 'NANO18_102_RegInputs'

config.Site.storageSite = 'T3_CH_PSI'
config.Site.whitelist = ['T2_CH_CERN', 'T1_IT_CNAF', 'T2_IT_Rome', 'T2_CH_CSCS']
config.Site.blacklist = ['T2_US_Nebraska']
