from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'ggHHbbgg_resonant700'
#config.General.requestName = 'ZHbbll'
config.General.requestName = 'V25_energyRings_TTbar_JECv10_ptd'
config.General.workArea = '/afs/cern.ch/work/n/nchernya/CMSSW_8_0_25_log//crab_projects'
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
                 #             '../ttbar-G25-500k-13d-300t.weights.xml',	
                              '../gravall-v25.weights.xml',	
			      '../TMVA_blikelihood_vbf_cmssw76_h21trained.weights.xml'
]
#config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/GluGluToRadionToHHTo2B2G_M-500_narrow_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/GluGluToRadionToHHTo2B2G_M-700_narrow_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset = '/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#config.Data.inputDataset ='/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
#config.Data.inputDataset = ''
#config.Data.inputDataset = ''
#config.Data.inputDataset = ''
#config.Data.inputDataset = ''
#config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 50000
#config.Data.unitsPerJob = 10000   # for HT bin 400-600 ext
#config.Data.totalUnits = 80000000
#config.Data.unitsPerJob = 10000
#config.Data.totalUnits = 500000
config.Data.allowNonValidInputDataset = True # to run on datasets in PRODUCTION
config.Data.outLFNDirBase = '/store/group/phys_higgs/vbfHbb/V25/'
config.Data.publication = False
config.Data.outputDatasetTag = 'VHBB_HEPPY_V25_energyRings_JECv10_ptd'

config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"
#config.Site.storageSite = "T2_IT_Pisa"
#config.Site.storageSite = "T3_CH_PSI"

#config.Data.ignoreLocality = True
#config.Data.ignoreLocality = True
#config.Data.ignoreLocality = True
#config.Data.ignoreLocality = True
#config.Data.ignoreLocality = True
