from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'VHBB_V25_interference'
config.General.workArea = '/afs/cern.ch/work/n/nchernya/CMSSW_8_0_25_log//crab_projects/'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'heppy_crab_fake_pset.py'
config.JobType.scriptExe = 'heppy_crab_script.sh'
import os
#os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")

os.system("tar czf python.tar.gz --directory $CMSSW_BASE python `find $CMSSW_BASE/src -name python | perl -pe s#$CMSSW_BASE/## `")
#onfig.JobType.sendPythonFolder = True
config.JobType.maxMemoryMB = 2500
config.JobType.maxJobRuntimeMin = 2000
config.JobType.inputFiles = ['heppy_config.py',
                             'heppy_crab_script.py',
                             'python.tar.gz',
                             #'MVAJetTags_620SLHCX_Phase1And2Upgrade.db',
                             'combined_cmssw.py',
                             '../vhbb.py',
                              '../vhbb_combined.py',
                             'TMVAClassification_BDT.weights.xml',
                             'puData.root',
                             'puDataMinus.root',
                             'puDataPlus.root',
                             '../triggerEmulation.root',
                             'puMC.root',
                              'json.txt',
                              #"../Zll-spring15.weights.xml",
                              #"../Wln-spring15.weights.xml",
                              #"../Znn-spring15.weights.xml",
                              #"../VBF-spring15.weights.xml",
                              #"../ttbar-fall15_TargetGenOverPt_GenPtCut0.weights.xml",
			      #'../ttbar-spring16-80X.weights.xml',
#			      '../ttbar-pumoriond17--500k-13d-300t.weights.xml',
                              '../ttbar-G25-500k-13d-300t.weights.xml',	
			      '../TMVA_blikelihood_vbf_cmssw76_h21trained.weights.xml'
]
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/LLJJ_INT_SM_5f_LO_13TeV_madgraph-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
#config.Data.totalUnits = 8000
config.Data.allowNonValidInputDataset = True # to run on datasets in PRODUCTION
config.Data.outLFNDirBase = '/store/group/phys_higgs/vbfHbb/V25/'
config.Data.publication = True
config.Data.outputDatasetTag = 'VHBB_HEPPY_V25'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
#config.Site.storageSite = "T2_IT_Pisa"
#config.Site.storageSite = "T3_CH_PSI"

#config.Data.ignoreLocality = True
