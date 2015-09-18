from VHbbAnalysis.Heppy.vhbbobj import Jet_blikeVBF

from math import *
import ROOT
import array

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
		for index,j in enumerate(event.jetsForHiggs,start=0) :
			self.Jet_pt[0]=j.pt()
			self.Jet_eta[0]=abs(j.eta())
			self.Jet_btagCSV[0]=j.btagcsv()
			self.Jet_pt_idx[0]=index
			self.Jet_eta_idx[0]=
			self.Jet_btagCSV_idex=
			j.blike_VBF[0]=self.reader.EvaluateMVA(self.name)[0]
		
		
			
