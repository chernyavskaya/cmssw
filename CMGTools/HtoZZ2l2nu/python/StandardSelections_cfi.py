import FWCore.ParameterSet.Config as cms

from CMGTools.External.puJetIDAlgo_cff import PhilV1, PuJetIdMinMVA, PuJetIdOptMVA
from CMGTools.HtoZZ2l2nu.TriggerSequences_cff import getTriggerPaths
DoubleElectronTrigs, DoubleMuTrigs, MuEGTrigs, PhotonTrigs, SingleMuTrigs, mcTrigs = getTriggerPaths(version=2012)

# base values for trigger event
BaseTriggerSelection = cms.PSet( source = cms.InputTag("TriggerResults::HLT"),
                                 triggerPaths = cms.PSet( gamma=cms.vstring(PhotonTrigs),
                                                          ee=cms.vstring(DoubleElectronTrigs),
                                                          mumu=cms.vstring(DoubleMuTrigs),
                                                          emu=cms.vstring(MuEGTrigs),
                                                          singleMu=cms.vstring(SingleMuTrigs)
                                                          )
                                 )

# base values for the vertex selection ------------------------------------------
BaseGeneratorSelection = cms.PSet( source = cms.InputTag("prunedGen"),
                                   filterId = cms.int32(25),
                                   genJets=cms.InputTag("selectedPatJetsPFlow"),
                                   puReweight=cms.InputTag("puWeights:puWeight")
                                   )


# base values for the vertex selection ------------------------------------------
BaseVertexSelection = cms.PSet( source = cms.InputTag("offlinePrimaryVertices"),
                                beamSpot = cms.InputTag("offlineBeamSpot"),
                                maxZ = cms.double(24),
                                maxRho = cms.double(2.0),
                                minNDOF = cms.int32(4)
                                )

# base values for muon selection ----------------------------------------------
Base2011MuonsSelection = cms.PSet( source = cms.InputTag("selectedPatMuons"), #PFlow"),
                                   minPt = cms.double(20),
                                   maxEta = cms.double(2.4),
                                   requireGlobal = cms.bool(True),
                                   minValidMuonHits=cms.int32(1),
                                   minMatchingMuonStations = cms.int32(2),
                                   minValidTrackerHits = cms.int32(11),
                                   minPixelHits = cms.int32(1),
                                   maxTrackChi2 = cms.double(10),
                                   maxRelPtUncertainty = cms.double(0.1),
                                   maxD0=cms.double(0.02),
                                   maxDz=cms.double(0.1),
                                   id = cms.string(""),
                                   maxRelIso = cms.double(999999.),#0.15),
                                   applySoftMuonIsolationVeto=cms.bool(False),
                                   usePFIso = cms.bool(False),
                                   doDeltaBetaCorrection = cms.bool(False)
                               )

BaseMuonsSelection = cms.PSet( source = cms.InputTag("selectedPatMuons"), #PFlow"),
                               minPt = cms.double(20),
                               maxEta = cms.double(2.4),
                               requireGlobal = cms.bool(False),
                               minValidMuonHits=cms.int32(0),
                               minMatchingMuonStations = cms.int32(0),
                               minValidTrackerHits = cms.int32(0),
                               minPixelHits = cms.int32(0),
                               maxTrackChi2 = cms.double(9999.),
                               maxRelPtUncertainty = cms.double(99999.),
                               maxD0=cms.double(1.0),
                               maxDz=cms.double(1.0),
                               id = cms.string(""),
                               maxRelIso = cms.double(999999.),#0.15),
                               applySoftMuonIsolationVeto=cms.bool(False),
                               usePFIso = cms.bool(False),
                               doDeltaBetaCorrection = cms.bool(False)
                               )

# base values for loose muon selection ----------------------------------------------
BaseLooseMuonsSelection = BaseMuonsSelection.clone( minPt = cms.double(10) )
BaseSoftMuonsSelection  = BaseMuonsSelection.clone( minPt = cms.double(3),
                                                    requireGlobal = cms.bool(False),
                                                    minValidMuonHits=cms.int32(0),
                                                    minMatchingMuonStations = cms.int32(0),
                                                    minValidTrackerHits = cms.int32(11),
                                                    minPixelHits = cms.int32(0),
                                                    maxTrackChi2 = cms.double(9999.),
                                                    maxRelPtUncertainty = cms.double(9999.),
                                                    maxD0=cms.double(0.2),
                                                    maxDz=cms.double(0.2),
                                                    id = cms.string("TMLastStationAngTight"),
                                                    maxRelIso = cms.double(999999.),
                                                    applySoftMuonIsolationVeto=cms.bool(True),
                                                    usePFIso = cms.bool(False),
                                                    doDeltaBetaCorrection = cms.bool(False) )
                                                    

# base values for photon selection ----------------------------------------------
BasePhotonsSelection = cms.PSet( source = cms.InputTag("photons"),
                                 conversions=cms.InputTag("allConversions"),
                                 gsfElectrons = cms.InputTag("gsfElectrons"),
                                 #cf. https://twiki.cern.ch/twiki/bin/view/CMS/RegressionSCCorrections
                                 scCorrector = cms.string("${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/PhoEnRegress.root"),
                                 rho25 = cms.InputTag("kt6PFJetsForIso:rho"),
                                 ebrechits = cms.InputTag("reducedEcalRecHitsEB"),
                                 eerechits = cms.InputTag("reducedEcalRecHitsEE"),
                                 minEt = cms.double(5), 
                                 maxEta = cms.double(2.5),
                                 maxHoE = cms.double(0.05),
                                 minSipipEB = cms.double(0.001),
                                 minSihihEB = cms.double(0.001),
                                 maxSihihEB = cms.double(0.011),
                                 maxSihihEE = cms.double(0.03),
                                 trkIsoCoeffsEB = cms.vdouble(2.0,  0.001,  0.0167),
                                 trkIsoCoeffsEE = cms.vdouble(2.0,  0.001,  0.0032),
                                 ecalIsoCoeffsEB = cms.vdouble(4.2, 0.006,  0.183),
                                 ecalIsoCoeffsEE = cms.vdouble(4.2, 0.006,  0.090),
                                 hcalIsoCoeffsEB = cms.vdouble(2.2, 0.0025, 0.062),
                                 hcalIsoCoeffsEE = cms.vdouble(2.2, 0.025,  0.180),
                                 trackSource = cms.InputTag("generalTracks"),
                                 gsfTrackSource = cms.InputTag("gsfElectronTracks")
                                 )

# base values for electron selection ----------------------------------------------
BaseElectronsSelection = cms.PSet( source = cms.InputTag("selectedPatElectrons"), #PFlow"),
                                   #cf. https://twiki.cern.ch/twiki/bin/view/CMS/RegressionSCCorrections
                                   scCorrector = cms.string("${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/EleEnRegress.root"),
                                   minPt = cms.double(20),
                                   maxEta = cms.double(2.5),
                                   vetoTransitionElectrons = cms.bool(True),
                                   #these cuts correspond to VBTF80 https://twiki.cern.ch/twiki/bin/view/CMS/VbtfEleID2011 with tigher hoe for EE
                                   vbtf2011 = cms.PSet( maxSihih     = cms.vdouble(0.01,  0.03),
                                                        maxDphiTrack = cms.vdouble(0.06,  0.04),
                                                        maxDetaTrack = cms.vdouble(0.004, 0.007),
                                                        maxHoE       = cms.vdouble(0.04,  0.1),
                                                        maxD0        = cms.vdouble(0.02,  0.02),
                                                        maxDz        = cms.vdouble(0.1,   0.1),
                                                        maxTrackLostHits = cms.vint32(0,  0),
                                                        applyConversionVetoFrom = cms.string("simpleEleId80relIso")
                                                        ),
                                   maxRelIso    = cms.double(999999.),#0.1),
                                   minDeltaRtoMuons = cms.double(0.1),
                                   usePFIso = cms.bool(False),
                                   doDeltaBetaCorrection = cms.bool(False)
                                   )

# base values for electron selection ----------------------------------------------
BaseLooseElectronsSelection = BaseElectronsSelection.clone(minPt = cms.double(10))

#my base values for jet selection -------------------------------------------------
BaseJetSelection = cms.PSet( source = cms.InputTag("selectedPatJetsPFlow"),
                             rho = cms.InputTag("kt6PFJetsPFlow:rho"),
                             minPt = cms.double(10),
                             maxEta = cms.double(5.0),
                             minDeltaRtoLepton = cms.double(0.4),
                             puJetIds = cms.VPSet( PhilV1,
                                                   PuJetIdMinMVA,
                                                   PuJetIdOptMVA
                                                   )
                             )
AssocJetSelection = BaseJetSelection.clone(source = cms.InputTag("selectedPatJetsPFlowNoPuSub") )

# base values for the dilepton selection ------------------------------------------
BaseDileptonSelection = cms.PSet( minDileptonMass = cms.double(0),
                                  maxDileptonMass = cms.double(7000),
                                  maxDz = cms.double(1.0)
                                  )

# base values for met selection -----------------------------------------------------
BaseMetSelection = cms.PSet( source = cms.InputTag("patMETsPFlow"),
                             trksource = cms.InputTag("trackMetProducer"),
                             hzzmetSources = cms.VInputTag("ClusteredPFMetProducer:assoc",                #1
                                                           "ClusteredPFMetProducer:standard",             #2  
                                                           "ClusteredPFMetProducer:central",              #3
                                                           "ClusteredPFMetProducer:cleaned",              #4 
                                                           "ClusteredPFMetProducer:assocCharged",         #5 
                                                           "ClusteredPFMetProducer:assocWithFwd",         #6
                                                           "ClusteredPFMetProducer:assoc",             #7  //to be replaced by something else
                                                           "ClusteredPFMetProducer:assocWithFwd",      #8  //to be replaced by something else
                                                           "ClusteredPFMetProducer:assoc",             #9  //to be replaced by something else
                                                           "ClusteredPFMetProducer:assocWithFwd",      #10 //to be replaced by something else
                                                           "ClusteredPFMetProducer:assocBeta",            #11
                                                           "ClusteredPFMetProducer:assocWithFwdBeta"),    #12
                             pfCands = cms.InputTag("particleFlow"),
                             pvAssocCandidatesSource = cms.InputTag("ClusteredPFMetProducer:pvAssocCandidates"),
                             sumEtSources = cms.InputTag("ClusteredPFMetProducer:globalPfMetSums")
                             )

