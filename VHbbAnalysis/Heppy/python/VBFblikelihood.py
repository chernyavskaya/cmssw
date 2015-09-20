from VHbbAnalysis.Heppy.vhbbobj import Jet_blikeVBF
from PhysicsTools.HeppyCore.utils.deltar import deltaPhi

from math import *
import ROOT
import array
from ROOT import TLorentzVector

class VBFblikelihood :
	def __init__(self,weightfile,name) :
		reader = ROOT.TMVA.Reader()
		self.Jet_pt=array.array('f',[0])
		self.Jet_eta=array.array('f',[0])
		self.Jet_btagCSV=array.array('f',[0])
		self.Jet_pt_idx=array.array('f',[0])
		self.Jet_eta_idx=array.array('f',[0])
		self.Jet_btagCSV_idx=array.array('f',[0])
		reader.AddVariable("Jet_pt",self.Jet_pt);
		reader.AddVariable("Jet_eta",self.Jet_eta);
		reader.AddVariable("Jet_btagCSV",self.Jet_btagCSV);
		reader.AddVariable("Jet_pt_idx",self.Jet_pt_idx);
		reader.AddVariable("Jet_eta_idx",self.Jet_eta_idx);
		reader.AddVariable("Jet_btagCSV_idx",self.Jet_btagCSV_idx);
      reader.BookMVA(name,weightfile)
      self.reader=reader
      self.name=name


	def evaluateBlikelihood(self, event) :
		for index,j in enumerate(event.jets,start=0) :
			j.blike_VBF[0]=-2
			if j.btagCSV>1 : j.btagCSV=1
			if j.btagCSV<0 : j.btagCSV=0
		if genWeight>0 :# how the variable is called?
			if ((j[0].pt()>92.) and (j[1].pt()>76.) and (j[2].pt()>64) and (j[3].pt()>30)):
				if (nJet>=4):# how the variable is called?
					jetsForHiggs4=[]
					jetsForHiggs4=[jet for index,jet in enumerate(event.jets) if (jet.jet_id>0) and (index<4) and (jet.pt()>20)]#how is Jet_id called here
					btag_max1=max(jetsForHiggs4.btagcsv())
					if btag_max1>0.7:
						bjet1=jetsForHiggs4[jetsForHiggs4.index.max(jetsForHiggs4.btagcsv())]
						other3jets=[]
						other3jets=[jet for index,jet in enumerate(jetsForHiggs4) if index!= btag_max1]
						maxEta1, maxEta2 = None, None
						maxDeltaEta = -1

						for i in range(len(other3jets)):
    						for j in range(i + 1, len(other3jets)):
        						d = abs(other3jets[i].eta() - other3jets[j].eta())
        						if d > maxDeltaEta:
            				maxEta1, maxEta2, maxDeltaEta = i, j, d
						if (maxEta1>=0) and (maxEta2>=0) :
							qjets=[]
							qjets.append(other3jets[maxEta1])
							qjets.append(other3jets[maxEta2])
							leftBjets=[jet for index,jet in enumerate(other3jets) if (index!=maxEta1) and (index!=maxEta2)]
							bjet2=leftBjets[0]

							Mqq=(qjets[0].p4()+qjets[1].p4()).M()
							bbDeltaPhi=abs(qjets[0].phi(),qjets[1].phi())
							qqDeltaEta=abs(qjets[0].eta()-qjets[1].eta())

							if (Mqq>460) and(qqDeltaEta>4.1)and (bbDeltaPhi<1.6) :
								if (HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v==1) : #the name of the variable

									loopMaxJet=7
									if nJet<7 : loopMaxJet=nJet  #variable name
									jetsForHiggsMax=[jet for index,jet in enumerate(event.jets) if (jet.jet_id>0) and (index<loopMaxJet) and (jet.pt()>20)]
									jetsEtaIdx=[sorted(range(len(jetsForHiggsMax)),key=lambda x:abs(jetsForHiggsMax[x].eta()))]
									jetsBtagIdx=[sorted(range(len(jetsForHiggsMax)),key=lambda x:abs(jetsForHiggsMax[x].btagCSV()))]

									for i	in range(0,loopMaxJet) :
										j=jetsForHiggsMax[i]
										if (j.jet_id<=0) : continue # the name of jet_id variable
										self.Jet_pt[0]=j.pt()# [0]?
										self.Jet_eta[0]=abs(j.eta())
										self.Jet_btagCSV[0]=j.btagcsv()
										self.Jet_pt_idx[0]=i
										self.Jet_eta_idx[0]=jetsEtaIdx[i]
										self.Jet_btagCSV_idx[0]=jetsBtagIdx[i]
										j.blike_VBF[0]=self.reader.EvaluateMVA(self.name)[0]

