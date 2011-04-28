import FWCore.ParameterSet.Config as cms
from CMGTools.Common.factories.cmgTau_cfi import cmgTau
from CMGTools.Common.factories.cmgDiObject_cfi import diObjectFactory

ditauFactory = diObjectFactory.clone(leg1Collection = cms.InputTag("cmgTau"),
                                     leg2Collection = cms.InputTag("cmgTau"),
                                     metCollection = cms.InputTag("")
                                     )


cmgDiTau = cms.EDFilter("DiTauPOProducer",
                        cfg = ditauFactory.clone(),
                        cuts =  cms.PSet( pt = cms.string("pt() > 0"),
                                          mass = cms.string(' 0 < mass() && mass() < 500'),             
                                          )
                        )

