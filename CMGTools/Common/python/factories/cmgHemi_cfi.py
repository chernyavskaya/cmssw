######################################################################

import FWCore.ParameterSet.Config as cms

from CMGTools.Common.selections.kinematics_cfi import kinematics

hemiFactory = cms.PSet(
    inputCollection = cms.InputTag("dummy"),
    maxCand = cms.uint32(1)
    )

cmgHemi = cms.EDFilter("HemispherePOProducer",
                       cfg = hemiFactory.clone(
    inputCollection = cms.InputTag("cmgPFJetSel")
    ),
                       cuts = cms.PSet(kinematics = kinematics.clone())
                       )

######################################################################
