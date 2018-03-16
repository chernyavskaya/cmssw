#! /usr/bin/env python
fastObjects=True

#Switch to True to produce x1, x2, id1, id2, pdf scale
doPDFVars = False

import ROOT
from DataFormats.FWLite import *
import PhysicsTools.HeppyCore.framework.config as cfg
from VHbbAnalysis.Heppy.vhbbobj import *
from PhysicsTools.HeppyCore.utils.deltar import deltaPhi
from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import * 
from VHbbAnalysis.Heppy.btagSF import btagSFhandle, get_event_SF

import logging
logging.basicConfig(level=logging.WARNING)

import os

cfg.Analyzer.nosubdir = True

treeProducer= cfg.Analyzer(
	class_object=AutoFillTreeProducer, 
	defaultFloatType = "F",
	verbose=False, 
	vectorTree = True,
        globalVariables	= [
                 NTupleVariable("puWeightUp", lambda ev : getattr(ev,"puWeightPlus",1.), help="Pileup up variation",mcOnly=True),
                 NTupleVariable("puWeightDown", lambda ev : getattr(ev,"puWeightMinus",1.), help="Pileup down variation",mcOnly=True),
                 NTupleVariable("json", lambda ev : getattr(ev,"json",True), help="Passing json selection"),
                 NTupleVariable("json_silver", lambda ev : getattr(ev,"json_silver",True), help="Passing silver json selection"),
                 NTupleVariable("nPU0", lambda ev : [bx.nPU() for bx in  ev.pileUpInfo if bx.getBunchCrossing()==0][0], help="nPU in BX=0",mcOnly=True),
                 NTupleVariable("nPVs", lambda ev : len(ev.goodVertices), help="total number of good PVs"),
		 NTupleVariable("Vtype", lambda ev : ev.Vtype, help="Event classification"),
		 NTupleVariable("VtypeSim", lambda ev : ev.VtypeSim, help="Event classification",mcOnly=True),
		 NTupleVariable("VMt", lambda ev : ev.V.goodMt, help="Transverse mass of the vector boson"),
		 NTupleVariable("HVdPhi", lambda ev : deltaPhi(ev.V.phi(),ev.H.phi()), help="Delta phi between Higgs and Z/W"),
		 NTupleVariable("fakeMET_sumet", lambda ev : ev.fakeMET.sumet, help="Fake SumET from Zmumu events removing muons"),
    	         NTupleVariable("bx",  lambda ev: ev.input.eventAuxiliary().bunchCrossing(), int, help="bunch crossing number"),
		 NTupleVariable("caloMetPt",  lambda ev: ev.met.caloMETPt(), float, help="calo met pt "),
		 NTupleVariable("caloMetPhi",  lambda ev: ev.met.caloMETPhi(), float, help="calo met phi"),
		 NTupleVariable("rho",  lambda ev: ev.rho, float, help="kt6PFJets rho"),
		 NTupleVariable("rhoN",  lambda ev: ev.rhoN, float, help="rho with neutrals only"),
		 NTupleVariable("rhoCHPU",  lambda ev: ev.rhoCHPU, float, help="rho with charged pileup particles"),
		 NTupleVariable("rhoCentral",  lambda ev: ev.rho, float, help="rho central"),
          
		 NTupleVariable("deltaR_jj",  lambda ev: deltaR(ev.hJets[0].eta(),ev.hJets[0].phi(),ev.hJets[1].eta(),ev.hJets[1].phi()) if len(ev.hJets) > 1 else -1, float, help="deltaR higgsJets"),
		 NTupleVariable("lheNj",  lambda ev: ev.lheNj, float,mcOnly=True, help="number of jets at LHE level"),
                 NTupleVariable("lheNb",  lambda ev: ev.lheNb, float,mcOnly=True, help="number of b-jets at LHE level"),
                 NTupleVariable("lheNc",  lambda ev: ev.lheNc, float,mcOnly=True, help="number of c-jets at LHE level"),
                 NTupleVariable("lheNg",  lambda ev: ev.lheNg, float,mcOnly=True, help="number of gluon jets at LHE level"),
                 NTupleVariable("lheNl",  lambda ev: ev.lheNl, float,mcOnly=True, help="number of light(uds) jets at LHE level"),
		 NTupleVariable("lheV_pt",  lambda ev: ev.lheV_pt, float,mcOnly=True, help="Vector pT at LHE level"),
                 NTupleVariable("lheHT",  lambda ev: ev.lheHT, float,mcOnly=True, help="HT at LHE level"),
                 #NTupleVariable("LHE_originalWeight",  lambda ev: ev.LHE_originalWeight, float, mcOnly=True, help="LHE original weight (for normalisation of LHE_weights)"),
                 NTupleVariable("genTTHtoTauTauDecayMode", lambda ev: ev.genTTHtoTauTauDecayMode, int,mcOnly=True, help="gen level ttH, H -> tautau decay mode"),        		 
		#Soft Activity vars
#                 NTupleVariable("totSoftActivityJets2", lambda ev: len([ x for x in ev.softActivityJets if x.pt()> 2 ] ), int, help="number of jets from soft activity with pt>2Gev"),
 #                NTupleVariable("totSoftActivityJets5", lambda ev: len([ x for x in ev.softActivityJets if x.pt()> 5 ] ), int, help="number of jets from soft activity with pt>5Gev"),
  #               NTupleVariable("totSoftActivityJets10", lambda ev: len([ x for x in ev.softActivityJets if x.pt()> 10 ] ), int, help="number of jets from soft activity with pt>10Gev"),
                 NTupleVariable("ttCls",  lambda ev: getattr(ev, "ttbarCls", -1), float,mcOnly=True, help="ttbar classification via GenHFHadronMatcher"),
                 NTupleVariable("heavyFlavourCategory",  lambda ev: getattr(ev, "ttbarCategory", -1), float,mcOnly=True, help="ttbar classification via GenHFHadronMatcher and the official GenTtbarClassifier code"),
		 NTupleVariable("mhtJet30",  lambda ev : ev.mhtJet30, help="mht with jets30"),
		 NTupleVariable("mhtPhiJet30",  lambda ev : ev.mhtPhiJet30, help="mht phi with jets30"),
		 NTupleVariable("htJet30",  lambda ev : ev.htJet30, help="ht  with jets30"),
                 NTupleVariable("met_sig",  lambda ev : ev.met.significance(), help="met significance from MET::significance() method"),
                 NTupleVariable("met_covXX",  lambda ev : ev.met.getSignificanceMatrix().At(0,0), help="xx element of met covariance matrix"),
		 NTupleVariable("met_covXY",  lambda ev : ev.met.getSignificanceMatrix().At(0,1), help="xy element of met covariance matrix"),
		 NTupleVariable("met_covYY",  lambda ev : ev.met.getSignificanceMatrix().At(1,1), help="yy element of met covariance matrix"),
		 NTupleVariable("met_rawpt",  lambda ev : ev.met.uncorPt(), help="raw met"),
		 NTupleVariable("metPuppi_pt",  lambda ev : ev.metPuppi.pt(), help="met from Puppi"),
		 NTupleVariable("metPuppi_phi",  lambda ev : ev.metPuppi.phi(), help="met phi from Puppi"),
		 NTupleVariable("metPuppi_rawpt",  lambda ev : ev.metPuppi.uncorPt(), help="raw met from Puppi"),
		 NTupleVariable("metType1p2_pt",  lambda ev : ev.met.shiftedPt(12,2), help="type1type2met"),
		 NTupleVariable("tkMet_pt",  lambda ev : ev.quickTkMET.pt(), help="pt of quick tk E_{T}^{miss}"),
		 NTupleVariable("tkMet_phi",  lambda ev : ev.quickTkMET.phi(), help="phi of quick tk E_{T}^{miss}"),
#		 NTupleVariable("tkMet_pt",  lambda ev : ev.tkMet.pt(), help="E_{T}^{miss} from tracks"),
#		 NTupleVariable("tkMet_phi",  lambda ev : ev.tkMet.phi(), help="phi of E_{T}^{miss} from tracks"),
#		 NTupleVariable("tkMetPVchs_pt",  lambda ev : ev.tkMetPVchs.pt(), help="E_{T}^{miss} from tracks"),
#		 NTupleVariable("tkMetPVchs_phi",  lambda ev : ev.tkMetPVchs.phi(), help="phi of E_{T}^{miss} from tracks"),
		 NTupleVariable("isrJetVH",  lambda ev : ev.isrJetVH, help="Index of ISR jet in VH"),
#		 NTupleVariable("Flag_hbheIsoFilter",  lambda ev : ev.hbheFilterIso, help="hbheFilterIso, after rerun"),
#		 NTupleVariable("Flag_hbheFilterNew",  lambda ev : ev.hbheFilterNew, help="hbheFilterIso, after rerun"),
		 NTupleVariable("simPrimaryVertex_z", lambda ev: ev.genvertex, float,mcOnly=True, help="z coordinate of the simulated primary vertex"),
		 NTupleVariable("genHiggsDecayMode", lambda ev: ev.genHiggsDecayMode, float, mcOnly=True, help="decay mode of the Higgs boson"),
		 NTupleVariable("triggerEmulationWeight", lambda ev: ev.triggerEmulationWeight, float, mcOnly=True, help="Emulates SL/DL triggers"),
	],
	globalObjects = {
          "met"    : NTupleObject("met",     metType, help="PF E_{T}^{miss}, after default type 1 corrections"),
          "fakeMET"    : NTupleObject("fakeMET", fourVectorType, help="fake MET in Zmumu event obtained removing the muons"),
          "H"    : NTupleObject("H", fourVectorType, help="higgs"),
          "H_reg"    : NTupleObject("H_reg", fourVectorType, help="regressed higgs"),
          "H_reg_corrJECUp"    : NTupleObject("H_reg_corrJECUp", fourVectorType, help="regressed higgs for JECUp"),
          "H_reg_corrJECDown"    : NTupleObject("H_reg_corrJECDown", fourVectorType, help="regressed higgs for JECDown"),
          "H_reg_corrJERUp"    : NTupleObject("H_reg_corrJERUp", fourVectorType, help="regressed higgs for JERUp"),
          "H_reg_corrJERDown"    : NTupleObject("H_reg_corrJERDown", fourVectorType, help="regressed higgs for JERDown"),
	  "HCSV"    : NTupleObject("HCSV", fourVectorType, help="higgs CSV selection"),
	  "HCSV_reg"    : NTupleObject("HCSV_reg", fourVectorType, help="regresses higgs CSV selection"),
          "HCSV_reg_corrJECUp"    : NTupleObject("HCSV_reg_corrJECUp", fourVectorType, help="reg_corrresses higgs for JECUp CSV selection"),
          "HCSV_reg_corrJECDown"    : NTupleObject("HCSV_reg_corrJECDown", fourVectorType, help="reg_corrresses higgs for JECDown CSV selection"),
          "HCSV_reg_corrJERUp"    : NTupleObject("HCSV_reg_corrJERUp", fourVectorType, help="reg_corrresses higgs for JERUp CSV selection"),
          "HCSV_reg_corrJERDown"    : NTupleObject("HCSV_reg_corrJERDown", fourVectorType, help="regressed higgs for JERDown CSV selection"),
	  "HCMVAV2"    : NTupleObject("HCMVAV2", fourVectorType, help="higgs CMVAV2 selection"),
	  "HCMVAV2_reg"    : NTupleObject("HCMVAV2_reg", fourVectorType, help="regresses higgs CMVAV2 selection"),
          "HCMVAV2_reg_corrJECUp"    : NTupleObject("HCMVAV2_reg_corrJECUp", fourVectorType, help="reg_corrresses higgs for JECUp CMVAV2 selection"),
          "HCMVAV2_reg_corrJECDown"    : NTupleObject("HCMVAV2_reg_corrJECDown", fourVectorType, help="reg_corrresses higgs for JECDown CMVAV2 selection"),
          "HCMVAV2_reg_corrJERUp"    : NTupleObject("HCMVAV2_reg_corrJERUp", fourVectorType, help="reg_corrresses higgs for JERUp CMVAV2 selection"),
          "HCMVAV2_reg_corrJERDown"    : NTupleObject("HCMVAV2_reg_corrJERDown", fourVectorType, help="regressed higgs for JERDown CMVAV2 selection"),
          "HaddJetsdR08"    : NTupleObject("HaddJetsdR08", fourVectorType, help="higgs with cen jets added if dR<0.8 from hJetsCSV selection"),
          "V"    : NTupleObject("V", fourVectorType, help="z or w"),
          "softActivityJets"    : NTupleObject("softActivity", softActivityType, help="VBF soft activity variables"),
          "softActivityVHJets"    : NTupleObject("softActivityVH", softActivityType, help="VH soft activity variables"),
          "softActivityEWKJets"    : NTupleObject("softActivityEWK", softActivityType, help="EWK soft activity variables"),
          "l1MET"       : NTupleObject("l1MET",   twoVectorType , help="Stage-2 L1 trigger MET", mcOnly=False),        
       #   "l1MET2"       : NTupleObject("l1MET2",   twoVectorType , help="Stage-2 L1 trigger MET", mcOnly=False),   #l1MET2 is defined in "l1t::EtSum" but it is empty
          "l1MHT"       : NTupleObject("l1MHT",   twoVectorType , help="Stage-2 L1 trigger MHT", mcOnly=False),        
          "l1ET"       : NTupleObject("l1ET",   twoVectorType , help="Stage-2 L1 trigger ET", mcOnly=False),        
          "l1HT"       : NTupleObject("l1HT",   twoVectorType , help="Stage-2 L1 trigger HT", mcOnly=False),      
        },
	collections = {
		#standard dumping of objects
   	        "selectedLeptons" : NTupleCollection("selLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
#   	        "inclusiveLeptons" : NTupleCollection("incLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
		#old style stuff, can be removed at some point
   	        "vLeptons" : NTupleCollection("vLeptons", leptonTypeVHbb, 8, help="Leptons after the preselection"),
   	        "aLeptons" : NTupleCollection("aLeptons", leptonTypeVHbb, 8, help="Additional leptons, not passing the preselection"),
# now store only indices, this lines are left commented for possible debugging
#	        "hJets"       : NTupleCollection("hJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Higgs jets"),
#	        "aJets"       : NTupleCollection("aJets",     jetTypeVHbb, 8, sortDescendingBy = lambda jet : jet.btag('combinedSecondaryVertexBJetTags'),help="Additional jets"),
                "hjidx"       : NTupleCollection("hJidx",    objectInt, 2,help="Higgs jet indices"),
                "hjidxDiJetPtByCSV"       : NTupleCollection("hJidx_sortcsv",    objectInt, 2,help="Higgs jet indices within hJets with CSV sorting "),
                "ajidx"       : NTupleCollection("aJidx",    objectInt, 8,help="additional jet indices"),
                "hjidxCSV"       : NTupleCollection("hJCidx",    objectInt, 2,help="Higgs jet indices CSV"),
                "ajidxCSV"       : NTupleCollection("aJCidx",    objectInt, 8,help="additional jet indices CSV"),
                "hjidxCMVAV2"       : NTupleCollection("hJCMVAV2idx",    objectInt, 2,help="Higgs jet indices CMVAV2"),
                "ajidxCMVAV2"       : NTupleCollection("aJCMVAV2idx",    objectInt, 8,help="additional jet indices CMVAV2"),
                "cleanJetsAll"       : NTupleCollection("Jet",     jetTypeVHbb, 25, help="Cental+fwd jets after full selection and cleaning, sorted by b-tag"),
                "hjidxaddJetsdR08"       : NTupleCollection("hjidxaddJetsdR08",    objectInt, 5,help="Higgs jet indices with Higgs formed adding cen jets if dR<0.8 from hJetsCSV"),
                "ajidxaddJetsdR08"       : NTupleCollection("ajidxaddJetsdR08",    objectInt, 8,help="additional jet indices with Higgs formed adding cen jets if dR<0.8 from hJetsCSV"),
		"dRaddJetsdR08"       : NTupleCollection("dRaddJetsdR08",    objectFloat, 5,help="dR of add jet with Higgs formed adding cen jets if dR<0.8 from hJetsCSV"),        
#                "discardedJets"       : NTupleCollection("DiscardedJet",     jetTypeVHbb, 15, help="jets that were discarded"),
                "inclusiveTaus"  : NTupleCollection("TauGood", tauTypeVHbb, 25, help="Taus after the preselection"),
                "softActivityJets"    : NTupleCollection("softActivityJets", fourVectorType, 5, help="jets made for soft activity"),
                "softActivityVHJets"    : NTupleCollection("softActivityVHJets", fourVectorType, 5, help="jets made for soft activity VH version"),
                "softActivityEWKJets"    : NTupleCollection("softActivityEWKJets", fourVectorType, 5, help="jets made for soft activity EWK version"),
                "goodVertices"    : NTupleCollection("primaryVertices", primaryVertexType, 4, help="first four PVs"), 

		#dump of gen objects
                #"generatorSummary"    : NTupleCollection("GenSummary", genParticleWithLinksType, 30, help="Generator summary, see description in Heppy GeneratorAnalyzer",mcOnly=True),
                "genJets"    : NTupleCollection("GenJet",   genJetType, 15, help="Generated jets with hadron matching, sorted by pt descending",filter=lambda x: x.pt() > 20,mcOnly=True),
                "vh_genHiggsSisters"    : NTupleCollection("GenHiggsSisters",     genParticleType, 4, help="Sisters of the Higgs bosons",mcOnly=True),
                "vh_gentopquarks"    : NTupleCollection("GenTop",     genTopType, 4, help="Generated top quarks from hard scattering",mcOnly=True),
                "vh_genallstatus2bhadrons"    : NTupleCollection("GenStatus2bHad",     genParticleType, 15, help="Generated Status 2 b Hadrons",mcOnly=True),
                "gennusFromTop"    : NTupleCollection("GenNuFromTop",     genParticleType, 4, help="Generated neutrino from t->W decay",mcOnly=True),
                "vh_genbquarksFromH"      : NTupleCollection("GenBQuarkFromH",  genParticleType, 4, help="Generated bottom quarks from Higgs decays",mcOnly=True),
                "vh_genbquarksFromTop"      : NTupleCollection("GenBQuarkFromTop",  genParticleType, 4, help="Generated bottom quarks from top decays",mcOnly=True),
                "vh_genbquarksFromHafterISR"      : NTupleCollection("GenBQuarkFromHafterISR",  genParticleType, 4, help="Generated bottom quarks from Higgs decays",mcOnly=True),
                "vh_gengluonfromb"      : NTupleCollection("GenGluonFromB",  genParticleType, 4, help="Generated gluons from b-quarks",mcOnly=True),
                "vh_gengluonfromt"      : NTupleCollection("GenGluonFromTop",  genParticleType, 4, help="Generated gluons from top quarks",mcOnly=True),
                "vh_genwzquarks"     : NTupleCollection("GenWZQuark",   genParticleType, 6, help="Generated quarks from W/Z decays",mcOnly=True),
                "vh_genleps"         : NTupleCollection("GenLep",     genParticleWithAncestryType, 6, help="Generated leptons from W/Z/Higgs decays",mcOnly=True),
		"gennus"         : NTupleCollection("GenNu",     genParticleWithAncestryType, 6, help="Generated neutrino from W/Z decays",mcOnly=True),
                "gentaus"         : NTupleCollection("GenTaus",     genParticleWithAncestryType, 6, help="Generated taus",mcOnly=True),
		"gennusFromTau" : NTupleCollection("GenNuFromTau",     genParticleType, 8, help="Generated neutrino from tau decay",mcOnly=True),
		
		"vh_genlepsRecovered"         : NTupleCollection("GenLepRecovered",     genParticleType, 4, help="Generated leptons from recovered W/Z decays",mcOnly=True),
                "genlepsFromTop"         : NTupleCollection("GenLepFromTop",     genParticleType, 4, help="Generated leptons from t->W decays",mcOnly=True),
		"gentauleps"      : NTupleCollection("GenLepFromTau", genParticleType, 6, help="Generated leptons from decays of taus from W/Z/Higgs decays",mcOnly=True),
		"vh_genHiggsBosons"   : NTupleCollection("GenHiggsBoson", genParticleType, 4, help="Generated Higgs boson ",mcOnly=True),
		#"genZbosonsToLL"  : NTupleCollection("GenZbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		#"genWbosonsToLL"  : NTupleCollection("GenWbosonsToLL", genParticleType, 6, help="Generated W or Z bosons decaying to LL"),
		"vh_genvbosons"       : NTupleCollection("GenVbosons", genParticleType, 6, help="Generated W or Z bosons, mass > 30",mcOnly=True),
                "vh_genvbosonsRecovered"  : NTupleCollection("GenVbosonsRecovered", genParticleType, 6, help="Generated W or Z bosons recovered from daughters, mass > 30",mcOnly=True),
		"pileUpVertex_z"       : NTupleCollection("pileUpVertex_z",    objectFloat, 5,help="z position of hardest pile-up collisions",mcOnly=True),        
		"pileUpVertex_ptHat"   : NTupleCollection("pileUpVertex_ptHat",    objectFloat, 5,help="z position of hardest pile-up collisions",mcOnly=True),        
		"LHE_weights_scale"       : NTupleCollection("LHE_weights_scale",   weightsInfoType , 6 ,help="LHE weights for scale variation", mcOnly=True),        
		"LHE_weights_pdf"       : NTupleCollection("LHE_weights_pdf",   weightsInfoType , 103 ,help="LHE weights for pdf variation (NNPDF)", mcOnly=True),        
                "HTXSRivetProducer_cat0"       : NTupleCollection("HTXSRivetProducer_cat0",   objectInt , 1 ,help="HTXSRivetProducer event category step 0", mcOnly=True),
                "HTXSRivetProducer_cat1"       : NTupleCollection("HTXSRivetProducer_cat1",   objectInt , 1 ,help="HTXSRivetProducer event category step 1", mcOnly=True),
                "l1Jets"       : NTupleCollection("l1Jets",   l1CandidateType , 20 ,help="Stage-2 L1 trigger jets", mcOnly=False),
		"l1Taus"       : NTupleCollection("l1Taus",   l1CandidateType , 20 ,help="Stage-2 L1 trigger taus", mcOnly=False),        
		"l1Muons"       : NTupleCollection("l1Muons",   l1CandidateType , 20 ,help="Stage-2 L1 trigger muons", mcOnly=False),        
		"l1EGammas"       : NTupleCollection("l1EGammas",   l1CandidateType , 20 ,help="Stage-2 L1 trigger EGammas", mcOnly=False),        
	}
	)

#Create shifted MET Ntuples
metNames={y.__get__(ROOT.pat.MET):x for x,y in ROOT.pat.MET.__dict__.items() if  (x[-2:]=="Up" or x[-4:]=="Down")}
print "met Names", metNames
shifted_met_keys = ["met_shifted_{0}".format(n) for n in range(12)] #we do not need noShift I gueess
shifted_met_names = ["met_shifted_%s"%metNames[n] for n in range(12)] #we do not need noShift I gueess
shifted_mets = {mk: NTupleObject(nm, shiftedMetType, help="PF E_{T}^{miss}, after default type 1 corrections, shifted with %s" %mk) for mk,nm in zip(shifted_met_keys,shifted_met_names)}
treeProducer.globalObjects.update(shifted_mets)

btag_weights = {}
for algo in ["CSV", "CMVAV2"]:
	for syst in ["central", "up_jes", "down_jes", "up_lf", "down_lf", "up_hf", "down_hf", "up_hfstats1", "down_hfstats1", "up_hfstats2", "down_hfstats2", "up_lfstats1", "down_lfstats1", "up_lfstats2", "down_lfstats2", "up_cferr1", "down_cferr1", "up_cferr2", "down_cferr2"]:
		syst_name = "" if syst=="central" else ("_"+syst) 
		btag_weights["btagWeight"+algo+syst_name] = NTupleVariable("btagWeight"+algo+syst_name,
									   lambda ev, get_event_SF=get_event_SF, syst=syst, algo=algo, btagSFhandle=btagSFhandle : get_event_SF(ev.cleanJetsAll, syst, algo, btagSFhandle)
									   , float, mcOnly=True, help="b-tag "+algo+"continuous  weight, variating "+syst
									   )
treeProducer.globalVariables += list(btag_weights.values())

ZllKinFitResults = {}
for analysis in ["","corrJECUp", "corrJECDown", "corrJERUp", "corrJERDown"]:
	name = "ZllKinFit"+("" if analysis=="" else "_")+analysis
	ZllKinFitResults[name+"_mass"] =  NTupleVariable(name+"_mass", lambda ev, name=name : getattr(ev, name+"_mass", -99.), float, mcOnly=False, help="Zll kin fit mass (=mass of the HCSV jets after kinematic fit) for analysis "+analysis) 
	ZllKinFitResults[name+"_njet"] =  NTupleVariable(name+"_njet", lambda ev, name=name : getattr(ev, name+"_njet", -99), int, mcOnly=False, help="Zll kin fit njet (=number of jets used for kinematic fit) for analysis "+analysis) 
	ZllKinFitResults[name+"_status"] =  NTupleVariable(name+"_status", lambda ev, name=name : getattr(ev, name+"_status", -99), int, mcOnly=False, help="Zll kin fit status (=0 if success, =1 otherwise) for analysis "+analysis) 
treeProducer.globalVariables += list(ZllKinFitResults.values())

'''

for syst in ["JES", "LF", "HF", "HFStats1", "HFStats2", "LFStats1", "LFStats2", "cErr1", "cErr2"]:
	for sdir in ["Up", "Down"]:
		name = "bTagWeight"+syst+sdir
		btag_weights[name] = NTupleVariable("bTagWeight_" + syst + sdir,
			lambda ev, sname=syst+sdir: bweightcalc.calcEventWeight(
				ev.cleanJetsAll, kind="final", systematic=sname
			), float, mcOnly=True, help="b-tag CSV weight, variating "+syst+" "+sdir
		)
btag_weights["bTagWeight"] = NTupleVariable("bTagWeight",
	lambda ev: bweightcalc.calcEventWeight(
		ev.cleanJetsAll, kind="final", systematic="nominal"
	), float ,mcOnly=True, help="b-tag CSV weight, nominal"
)
#print list(btag_weights.values())
treeProducer.globalVariables += list(btag_weights.values())
'''

# Lepton Analyzer, take its default config and fix loose iso consistent with tight definition
from PhysicsTools.Heppy.analyzers.objects.LeptonAnalyzer import LeptonAnalyzer
LepAna = LeptonAnalyzer.defaultConfig
#LepAna.doElectronScaleCorrections = {
#	'data' : 'EgammaAnalysis/ElectronTools/data/76X_16DecRereco_2015', 
#	'GBRForest' : [os.environ['CMSSW_BASE']+'/src/VHbbAnalysis/Heppy/data/egamma/egamma_epComb_GBRForest_76X.root','gedelectron_p4combination_25ns'], 
#	'isSync' : False
#	}
LepAna.mu_isoCorr  = "deltaBeta"
#LepAna.loose_muon_isoCut = lambda muon : muon.relIso04 < 0.25 
LepAna.loose_muon_isoCut = lambda muon : (muon.relIso04 < 0.4 or muon.miniRelIso < 0.4)
LepAna.loose_muon_pt = 4.5
LepAna.ele_isoCorr = "rhoArea"
LepAna.loose_electron_isoCut = lambda electron : (electron.relIso03 < 0.4 or electron.miniRelIso < 0.4) 
#LepAna.loose_electron_isoCut = lambda electron : ele_mvaEleID_Trig_preselection(electron)
LepAna.loose_electron_pt = 6.5
LepAna.loose_electron_eta = 2.5
#LepAna.loose_electron_id = "mvaEleID-Spring15-25ns-Trig-V1-wp90"
LepAna.doMiniIsolation = True

from PhysicsTools.Heppy.analyzers.objects.VertexAnalyzer import VertexAnalyzer
VertexAna = VertexAnalyzer.defaultConfig
VertexAna.keepFailingEvents = True
VertexAna.scores = 'offlineSlimmedPrimaryVertices'

from PhysicsTools.Heppy.analyzers.objects.PhotonAnalyzer import PhotonAnalyzer
PhoAna = PhotonAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.TauAnalyzer import TauAnalyzer
TauAna = TauAnalyzer.defaultConfig
TauAna.inclusive_ptMin = 18.
TauAna.inclusive_etaMax = 2.5
TauAna.inclusive_dxyMax = 1000.
TauAna.inclusive_dzMax = 0.4
TauAna.inclusive_vetoLeptons = False
TauAna.inclusive_leptonVetoDR = 0.4
TauAna.inclusive_decayModeID = "decayModeFindingNewDMs"
TauAna.inclusive_tauID = "decayModeFindingNewDMs"
TauAna.inclusive_vetoLeptonsPOG = False
TauAna.inclusive_tauAntiMuonID = ""
TauAna.inclusive_tauAntiElectronID = ""
TauAna.inclusive_tauLooseID = "decayModeFindingNewDMs"

from PhysicsTools.Heppy.analyzers.objects.JetAnalyzer import JetAnalyzer
JetAna = JetAnalyzer.defaultConfig
JetAna.calculateSeparateCorrections = True # CV: needed for ttH prompt lepton MVA
JetAna.lepSelCut = lambda lep : (abs(lep.pdgId()) == 11 and lep.relIso03 < 0.4) or (abs(lep.pdgId()) == 13 and lep.relIso04 < 0.4)

from PhysicsTools.Heppy.analyzers.gen.LHEAnalyzer import LHEAnalyzer 
LHEAna = LHEAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.gen.LHEWeightAnalyzer import LHEWeightAnalyzer 
LHEWeightAna = LHEWeightAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.gen.GeneratorAnalyzer import GeneratorAnalyzer 
GenAna = GeneratorAnalyzer.defaultConfig
GenAna.allGenTaus = True

from VHbbAnalysis.Heppy.VHGeneratorAnalyzer import GeneratorAnalyzer as  VHGeneratorAnalyzer
VHGenAna = VHGeneratorAnalyzer.defaultConfig

from PhysicsTools.Heppy.analyzers.objects.METAnalyzer import METAnalyzer
METAna = METAnalyzer.defaultConfig
##METAna.metCollection = "slimmedMETs::EX"
METAna.metCollection = "slimmedMETs"
METAna.recalibrate = False
METAna.applyJetSmearing = False
METAna.doTkMet = False
METAna.doMetNoPU = False
METAna.doTkGenMet = False
METAna.includeTkMetPVLoose = False
METAna.includeTkMetPVTight = False

METPuppiAna = copy.copy(METAna)
METPuppiAna.metCollection = "slimmedMETsPuppi"
METPuppiAna.doMetNoPU = False
METPuppiAna.recalibrate = False
METPuppiAna.collectionPostFix = "Puppi"
METPuppiAna.copyMETsByValue = True
METPuppiAna.doTkMet = False
METPuppiAna.doMetNoPU = False

from PhysicsTools.Heppy.analyzers.core.PileUpAnalyzer import PileUpAnalyzer
PUAna = PileUpAnalyzer.defaultConfig

from VHbbAnalysis.Heppy.VHbbAnalyzer import VHbbAnalyzer
JetAna.jetPt = 15
JetAna.jetEta = 4.7
JetAna.doQG=True
JetAna.QGpath=os.environ['CMSSW_BASE']+"/src/PhysicsTools/Heppy/data/pdfQG_AK4chs_13TeV_cmssw8020_v2.root"
JetAna.recalibrateJets=True
JetAna.jecPath=os.environ['CMSSW_BASE']+"/src/VHbbAnalysis/Heppy/data/jec"
#JetAna.mcGT="Fall15_25nsV2_MC"
#JetAna.dataGT = "Fall15_25nsV2_DATA"

JetAna.mcGT="Summer16_23Sep2016V3_MC"
JetAna.dataGT="Summer16_23Sep2016GV3_DATA"
JetAna.addJECShifts=True
JetAna.addJERShifts=True

factorizedJetCorrections = [
        "AbsoluteStat",
        "AbsoluteScale",
        "AbsoluteFlavMap",
        "AbsoluteMPFBias",
        "Fragmentation",
        "SinglePionECAL",
        "SinglePionHCAL",
        "FlavorQCD",
        "TimePtEta",
        "RelativeJEREC1",
        "RelativeJEREC2",
        "RelativeJERHF",
        "RelativePtBB",
        "RelativePtEC1",
        "RelativePtEC2",
        "RelativePtHF",
        "RelativeBal",
        "RelativeFSR",
        "RelativeStatFSR",
        "RelativeStatEC",
        "RelativeStatHF",
        "PileUpDataMC",
        "PileUpPtRef",
        "PileUpPtBB",
        "PileUpPtEC1",
        "PileUpPtEC2",
        "PileUpPtHF",
        "PileUpMuZero",
        "PileUpEnvelope",
        "SubTotalPileUp",
        "SubTotalRelative",
        "SubTotalPt",
        "SubTotalScale",
        "SubTotalAbsolute",
        "SubTotalMC",
        "Total",
        "TotalNoFlavor",
        "TotalNoTime",
        "TotalNoFlavorNoTime",
        "FlavorZJet",
        "FlavorPhotonJet",
        "FlavorPureGluon",
        "FlavorPureQuark",
        "FlavorPureCharm",
        "FlavorPureBottom",
        "TimeRunBCD",
        "TimeRunEF",
        "TimeRunG",
        "TimeRunH",
        "CorrelationGroupMPFInSitu",
        "CorrelationGroupIntercalibration",
        "CorrelationGroupbJES",
        "CorrelationGroupFlavor",
        "CorrelationGroupUncorrelated",
]
JetAna.factorizedJetCorrections = factorizedJetCorrections

for jet_type in [jetTypeVHbb, patSubjetType, subjetcorrType]:
        for jet_corr in factorizedJetCorrections:
            for sdir in ["Up", "Down"]:
                name = jet_corr + sdir
                jet_type.variables += [
                    NTupleVariable(
                        "corr_{0}".format(name),
                        lambda x, name = name: getattr(x, 'corr{0}'.format(name), -99),
                        float,
                        mcOnly=True,
                        help=""
                    )
                ]

        
        
# HTT Subjets
for subjet in ["sjW1", "sjW2", "sjNonW"]:
    for jet_corr in factorizedJetCorrections:
        for sdir in ["Up", "Down"]:
            name = jet_corr + sdir
            httType.variables += [
                NTupleVariable(
                    "{0}_corr_{1}".format(subjet, name),
                        lambda x, name=name, subjet=subjet : getattr( getattr(x,subjet), 'corr{0}'.format(name), -99),
                    float,
                    mcOnly=True,
                    help=""
                )
            ]



# delta-beta corrected isolation for muons:
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
def mu_pfRelIso04(mu):
	return (mu.pfIsolationR04().sumChargedHadronPt + max( mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5 * mu.pfIsolationR04().sumPUPt,0.0)) / mu.pt()

# the MVA-triggering preseletion according to:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recommended_MVA_recipes_for_2015
def ele_mvaEleID_Trig_preselection(ele) : 
	return (ele.pt()>15 and 
		( ( abs(ele.superCluster().eta()) < 1.4442 and ele.full5x5_sigmaIetaIeta() < 0.012 and ele.hcalOverEcal() < 0.09 and (ele.ecalPFClusterIso() / ele.pt()) < 0.37 and (ele.hcalPFClusterIso() / ele.pt()) < 0.25 and (ele.dr03TkSumPt() / ele.pt()) < 0.18 and abs(ele.deltaEtaSuperClusterTrackAtVtx()) < 0.0095 and abs(ele.deltaPhiSuperClusterTrackAtVtx()) < 0.065 ) or 
		  ( abs(ele.superCluster().eta()) > 1.5660 and ele.full5x5_sigmaIetaIeta() < 0.033 and ele.hcalOverEcal() <0.09 and (ele.ecalPFClusterIso() / ele.pt()) < 0.45 and (ele.hcalPFClusterIso() / ele.pt()) < 0.28 and (ele.dr03TkSumPt() / ele.pt()) < 0.18 ) ) )

from VHbbAnalysis.Heppy.ttHLeptonMVAAnalyzer import ttHLeptonMVAAnalyzer
ttHLeptonMVA = cfg.Analyzer(
    verbose = False,
    class_object = ttHLeptonMVAAnalyzer,
)

from VHbbAnalysis.Heppy.TriggerEmulation import TriggerEmulationAnalyzer
trigemu = cfg.Analyzer(
    verbose = False,
    class_object = TriggerEmulationAnalyzer,
    calibrationFile="triggerEmulation.root",
    slEleSelection = lambda x : x.pt() > 25 and getattr(x,"mvaIdSpring16GeneralPurposePOG90",False) and ele_mvaEleID_Trig_preselection(x),
    slMuSelection = lambda x : x.pt() > 25 and x.muonID("POG_ID_Tight") and mu_pfRelIso04(x) < 0.15,
    dlEleSelection = lambda x : x.pt() > 15 and getattr(x,"mvaIdSpring16GeneralPurposePOG80",False) and ele_mvaEleID_Trig_preselection(x),
    dlMuSelection = lambda x : x.pt() > 15 and x.muonID("POG_ID_Loose") and mu_pfRelIso04(x) < 0.25,
)

VHbb = cfg.Analyzer(
    verbose=False,
    class_object=VHbbAnalyzer,
    wEleSelection = lambda x : x.pt() > 25 and getattr(x,"mvaIdSpring16GeneralPurposePOG90",False)     and ele_mvaEleID_Trig_preselection(x),
    wMuSelection = lambda x : x.pt() > 25 and x.muonID("POG_ID_Tight") and mu_pfRelIso04(x) < 0.15,
    zEleSelection = lambda x : x.pt() > 15 and getattr(x,"mvaIdSpring16GeneralPurposePOG80",False)     and ele_mvaEleID_Trig_preselection(x),
    zMuSelection = lambda x : x.pt() > 15 and x.muonID("POG_ID_Loose") and mu_pfRelIso04(x) < 0.25,
    zLeadingElePt = 20,
    zLeadingMuPt = 20,
    singleJetThreshold = 200,
    sumPtThreshold = 160.,
    higgsJetsPreSelection = lambda x: (( x.puJetId() > 0 and x.jetID('POG_PFID_Loose')) or abs(x.eta())>3.0 ) and x.pt() >  20 ,
    higgsJetsPreSelectionVBF = lambda x: x.pt() >  20 ,
#    higgsJetsPreSelectionVBF = lambda x: (( x.puJetId() > 0 and x.jetID('POG_PFID_Loose')) or abs(x.eta())>3.0 ) and x.pt() >  20,
    passall=False,
    doSoftActivityVH=True,
    doZllKinematicFit=True,
    doSoftActivityEWK=True,
    doVBF=True,
    regressions = [
        {"weight":"ttbar-G25-500k-13d-300t.weights.xml", "name":"jet0Regression", "vtypes":[0,1,2,3,4,5,-1]},
    ],
    regressionVBF = [
        {"weight":"ttbar-G25-500k-13d-300t.weights.xml", "name":"jet0Regression_vbf", "vtypes":[0,1,2,3,4,5,-1]}
    ],
    VBFblikelihood = {"weight":"TMVA_blikelihood_vbf_cmssw76_h21trained.weights.xml", "name":"BDGT"}
)

from VHbbAnalysis.Heppy.TTHtoTauTauAnalyzer import TTHtoTauTauAnalyzer
TTHtoTauTau = cfg.Analyzer(
    verbose = False,
    class_object = TTHtoTauTauAnalyzer,
)
from VHbbAnalysis.Heppy.TTHtoTauTauGeneratorAnalyzer import TTHtoTauTauGeneratorAnalyzer
TTHtoTauTauGen = cfg.Analyzer(
    verbose = False,
    class_object = TTHtoTauTauGeneratorAnalyzer,
)

#from VHbbAnalysis.Heppy.HeppyShell import HeppyShell
#sh = cfg.Analyzer( class_object=HeppyShell)

from PhysicsTools.Heppy.analyzers.core.TriggerBitAnalyzer import TriggerBitAnalyzer

from VHbbAnalysis.Heppy.TriggerTable import triggerTable
from VHbbAnalysis.Heppy.TriggerTableData import triggerTable as triggerTableData

TrigAna = cfg.Analyzer(
    verbose = False,
    class_object = TriggerBitAnalyzer,
    triggerBits = triggerTable,  #default is MC, use the triggerTableData in -data.py files
   processName = 'HLT',
   fallbackProcessName = 'HLT2',
#   outprefix = 'HLT'
   )

from PhysicsTools.HeppyCore.framework.services.tfile import TFileService 
output_service = cfg.Service(
      TFileService,
      'outputfile',
      name="outputfile",
      fname='tree.root',
      option='recreate'
    )

from PhysicsTools.Heppy.analyzers.core.TriggerBitAnalyzer import TriggerBitAnalyzer
FlagsAna = TriggerBitAnalyzer.defaultEventFlagsConfig
FlagsAna.triggerBits.update( { "chargedHadronTrackResolutionFilter" : ["Flag_chargedHadronTrackResolutionFilter"], "muonBadTrackFilter" : ["Flag_muonBadTrackFilter"] , "GlobalTightHalo2016Filter" : ["Flag_globalTightHalo2016Filter"] } )


#from VHbbAnalysis.Heppy.hbheAnalyzer import *
#hbheAna = hbheAnalyzer.defaultConfig


### add trigger objects ####
from PhysicsTools.Heppy.analyzers.core.TriggerObjectsAnalyzer import TriggerObjectsAnalyzer
from VHbbAnalysis.Heppy.TriggerObjectsList import *
TriggerObjectsAna = TriggerObjectsAnalyzer.defaultConfig
TriggerObjectsAna.triggerObjectsCfgs = triggerObjectCollections

for collectionName in triggerObjectCollectionsFull.keys():
    treeProducer.collections["trgObjects_"+collectionName] = NTupleCollection("trgObjects_"+collectionName, triggerObjectsType, 5, help="")

for collectionName in triggerObjectCollectionsOnlyPt.keys():
    treeProducer.collections["trgObjects_"+collectionName] = NTupleCollection("trgObjects_"+collectionName, triggerObjectsOnlyPtType, 5, help="")

for collectionName in triggerObjectCollectionsOnlySize.keys():
    treeProducer.collections["trgObjects_"+collectionName] = NTupleCollection("trgObjects_"+collectionName, triggerObjectsNothingType , 5, help="")
#    treeProducer.globalVariables.append(NTupleVariable("trgObjects_"+collectionName+"_size", lambda ev : len(getattr(ev,"trgObjects_"+collectionName,[])), int, help="trigger objects size"))

### add L1 trigger objects ####
from PhysicsTools.Heppy.analyzers.core.L1TriggerAnalyzer import L1TriggerAnalyzer
L1TriggerAna = cfg.Analyzer(
    class_object = L1TriggerAnalyzer,
    processName = 'HLT',
)
###




from PhysicsTools.Heppy.analyzers.gen.PDFWeightsAnalyzer import PDFWeightsAnalyzer
PdfAna = cfg.Analyzer(PDFWeightsAnalyzer,
    PDFWeights = [],
    doPDFVars = doPDFVars
)

if doPDFVars:
    treeProducer.globalVariables += [
        NTupleVariable("pdf_x1",  lambda ev: ev.pdf_x1, float,mcOnly=True, help="PDF energy fraction of first parton"),
        NTupleVariable("pdf_x2",  lambda ev: ev.pdf_x2, float,mcOnly=True, help="PDF energy fraction of second parton"),
        NTupleVariable("pdf_id1",  lambda ev: ev.pdf_id1, int,mcOnly=True, help="PDF id of first parton"),
        NTupleVariable("pdf_id2",  lambda ev: ev.pdf_id2, int,mcOnly=True, help="PDF id of second parton"),
        NTupleVariable("pdf_scale",  lambda ev: ev.pdf_scale, float,mcOnly=True, help="PDF scale"),
    ]

TrigAna.unrollbits=True

from PhysicsTools.Heppy.analyzers.core.JSONAnalyzer import JSONAnalyzer
jsonAna = cfg.Analyzer(JSONAnalyzer,
      passAll=True
      )
silverJsonAna = cfg.Analyzer(JSONAnalyzer,
      passAll=True,
      json="silver.txt",
      suffix="_silver"
      )

sequence = [
    jsonAna,LHEAna,LHEWeightAna,FlagsAna,
    #hbheAna, 
    GenAna,VHGenAna,PUAna,TrigAna,
    VertexAna,LepAna,PhoAna,TauAna,JetAna,
    ttHLeptonMVA,METAna, METPuppiAna,
    PdfAna,
    VHbb,TTHtoTauTau,TTHtoTauTauGen,TriggerObjectsAna,L1TriggerAna,trigemu,treeProducer
]

from PhysicsTools.Heppy.utils.miniAodFiles import miniAodFiles
sample = cfg.MCComponent(
	files = [
		#"root://xrootd.ba.infn.it//store/mc/RunIIFall15MiniAODv1/TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/30000/045996FE-A19D-E511-B76D-D4AE526A0B47.root" ##ttbar
		#"root://xrootd.ba.infn.it//store/mc/RunIISpring16MiniAODv1/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/0899BDA9-AE01-E611-A239-008CFA05EA2C.root"
#		"root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_8_1_0_pre9/RelValTTbar_13/MINIAODSIM/PU25ns_81X_mcRun2_asymptotic_v2_hip0p8_mtoff-v1/10000/2A7336F1-D851-E611-AA11-003048D15D48.root"
#/store/mc/RunIISpring16MiniAODv2/GluGluToBulkGravitonToHHTo4B_M-550_narrow_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v2/90000/4E40D2E2-9E3A-E611-8C5B-00259081FB18.root"
#root://stormgf1.pi.infn.it:1094//store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/E8090432-8628-E611-8713-001EC9ADFDC9.root",
		#"root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/mc/RunIISpring16MiniAODv2/WW_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/0AF21AF1-121B-E611-B652-549F35AE4F88.root"
		# "root://stormgf1.pi.infn.it:1094//store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/00000/E8090432-8628-E611-8713-001EC9ADFDC9.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/184FBB8C-45C8-E611-B5C5-001EC9ADC726.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/601C84AE-60D0-E611-8A73-0025905A48D6.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/2029BD2C-7DCB-E611-8E08-0025904A87E2.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/GluGluHToBB_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/06ADFE29-08CA-E611-943C-0025905B860C.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/GluGluHToBB_M125_13TeV_powheg_herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/BA7BFE5F-93CC-E611-94E7-0CC47A57D164.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/02D08DD9-79B9-E611-B828-00266CFFBE14.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/3E7EED14-C5B6-E611-BC04-484D7E8DF0C6.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/04643CD6-CFC2-E611-9FD8-0CC47AD98F6E.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/TT_TuneEE5C_13TeV-powheg-herwigpp/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/02656FC1-B0B5-E611-B2F9-44A842CFCA27.root",
	#	"root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/WH_HToBB_WToLNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/7C89D7C0-B7D0-E611-AB35-0CC47A4C8EE8.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/ggZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8_CUETP8M1Up/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/D466A5FE-11CC-E611-AA18-842B2B17E3BA.root",
		# "root://cmsxrootd.fnal.gov///store/mc/RunIISummer16MiniAODv2/ZH_HToBB_ZToLL_M125_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/60000/28CA6CC9-FEC2-E611-BF37-008CFA5D2758.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_0.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_1.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/Heppy/Loop_drop_1/cmsswPreProcessing.root",
"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/Heppy/Loop02/cmsswPreProcessing.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_2.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_3.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_4.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_5.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_6.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_7.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_8.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_9.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_10.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_11.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_12.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_13.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_14.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_15.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_16.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_17.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_18.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_19.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_20.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_21.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_22.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_23.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_24.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_25.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_26.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_27.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_28.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_29.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_30.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_31.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_32.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_33.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_34.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_35.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_36.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_37.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_38.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_39.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_40.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_41.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_42.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_43.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_44.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_45.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_46.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_47.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_48.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_49.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_50.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_51.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_52.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_53.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_54.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_55.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_56.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_57.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_58.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_59.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_60.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_61.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_62.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_63.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_64.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_65.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_66.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_67.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_68.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_69.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_70.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_71.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_72.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_73.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_74.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_75.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_76.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_77.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_78.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_79.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_80.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_81.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_82.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_83.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_84.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_85.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_86.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_87.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_88.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_89.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_90.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_91.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_92.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_93.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_94.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_95.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_96.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_97.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_98.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_99.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_100.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_101.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_102.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_103.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_104.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_105.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_106.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_107.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_108.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_109.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_110.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_111.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_112.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_113.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_114.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_115.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_116.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_117.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_118.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_119.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_120.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_121.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_122.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_123.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_124.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_125.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_126.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_127.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_128.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_129.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_130.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_131.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_132.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_133.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_134.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_135.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_136.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_137.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_138.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_139.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_140.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_141.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_142.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_143.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_144.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_145.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_146.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_147.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_148.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_149.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_150.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_151.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_152.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_153.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_154.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_155.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_156.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_157.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_158.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_159.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_160.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_161.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_162.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_163.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_164.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_165.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_166.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_167.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_168.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_169.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_170.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_171.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_173.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_174.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_175.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_176.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_177.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_178.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_179.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_180.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_181.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_182.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_183.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_184.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_185.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_186.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_187.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_188.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_189.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_190.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_191.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_192.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_193.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_194.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_195.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_196.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_197.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_198.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_199.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_200.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_201.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_202.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_203.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_204.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_205.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_206.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_207.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_208.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_209.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_210.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_211.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_212.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_213.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_214.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_215.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_216.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_217.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_218.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_219.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_220.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_221.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_222.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_223.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_224.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_225.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_226.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_227.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_228.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_229.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_230.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_231.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_232.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_233.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_234.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_235.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_236.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_237.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_238.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_239.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_240.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_241.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_242.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_243.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_244.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_245.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_246.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_247.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_248.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_249.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_250.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_251.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_252.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_253.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_254.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_255.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_256.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_257.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_258.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_259.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_260.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_261.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_262.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_263.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_264.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_265.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_266.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_267.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_268.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_269.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_270.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_271.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_272.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_273.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_274.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_275.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_276.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_277.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_278.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_279.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_280.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_281.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_282.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_283.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_284.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_285.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_286.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_287.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_288.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_289.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_290.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_291.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_292.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_293.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_294.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_295.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_296.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_297.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_298.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_299.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_300.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_301.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_302.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_303.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_304.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_305.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_306.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_307.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_308.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_309.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_310.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_311.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_312.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_313.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_314.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_315.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_316.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_317.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_318.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_319.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_320.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_321.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_322.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_323.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_324.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_325.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_326.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_327.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_328.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_329.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_330.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_331.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_332.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_333.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_334.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_335.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_336.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_337.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_338.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_339.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_340.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_341.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_342.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_343.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_344.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_345.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_346.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_347.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_348.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_349.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_350.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_351.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_352.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_353.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_354.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_355.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_356.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_357.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_358.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_359.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_360.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_361.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_362.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_363.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_364.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_365.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_366.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_367.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_368.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_369.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_370.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_371.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_372.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_373.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_374.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_375.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_376.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_377.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_378.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_379.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_380.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_381.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_382.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_383.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_384.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_385.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_386.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_387.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_388.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_389.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_390.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_391.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_392.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_393.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_394.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_395.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_396.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_397.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_398.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_399.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_400.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_401.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_402.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_403.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_404.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_405.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_406.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_407.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_408.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_409.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_410.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_411.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_412.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_413.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_414.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_415.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_416.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_417.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_418.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_419.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_420.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_421.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_422.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_423.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_424.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_425.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_426.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_427.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_428.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_429.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_430.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_431.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_432.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_433.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_434.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_435.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_436.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_437.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_438.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_439.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_440.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_441.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_442.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_443.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_444.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_445.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_446.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_447.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_448.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_449.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_450.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_451.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_452.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_453.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_454.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_455.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_456.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_457.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_458.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_459.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_460.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_461.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_462.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_463.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_464.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_465.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_466.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_467.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_468.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_469.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_470.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_471.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_472.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_473.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_474.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_475.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_476.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_477.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_478.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_479.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_480.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_481.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_482.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_483.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_484.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_485.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_486.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_487.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_488.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_489.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_490.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_491.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_492.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_493.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_494.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_495.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_496.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_497.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_498.root",
#"/eos/cms/store/group/phys_higgs/vbfHbb/LHE_MINIAOD/EWKinterference/MINIAODSIM/Events_499.root",
#
    ],
    # files = ["cmsswPreProcessing.root"],
   # name="ZHLL125", isEmbed=False,
    name="EWKinterference", isEmbed=False,
    puFileMC="puMC.root",
    puFileData="puData.root", 
    puFileDataPlus="puDataPlus.root", 
    puFileDataMinus="puDataMinus.root", 
    splitFactor = 5
    )
sample.isMC=True




# the following is declared in case this cfg is used in input to the heppy.py script
selectedComponents = [sample]
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
config = cfg.Config( components = selectedComponents,
                     sequence = sequence, 
		     services = [output_service],
                     events_class = Events)

class TestFilter(logging.Filter):
    def filter(self, record):
        print record

# and the following runs the process directly 
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper 
    looper = Looper( 'Loop', config, nPrint = 0, nEvents = 50000)

    import time
    import cProfile
    p = cProfile.Profile(time.clock)
    p.runcall(looper.loop)
    p.print_stats()
    looper.write()
