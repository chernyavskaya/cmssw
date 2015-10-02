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
    reader.AddVariable("Jet_pt",self.Jet_pt)
    reader.AddVariable("Jet_eta",self.Jet_eta)
    reader.AddVariable("Jet_btagCSV",self.Jet_btagCSV)
    reader.AddVariable("Jet_pt_idx",self.Jet_pt_idx)
    reader.AddVariable("Jet_eta_idx",self.Jet_eta_idx)
    reader.AddVariable("Jet_btagCSV_idx",self.Jet_btagCSV_idx)
    reader.BookMVA(name,weightfile)
    self.reader=reader
    self.name=name


  def evaluateBlikelihood(self, event) :
    #print 'evaluate blikelihood is ok'
    nJet=len(event.cleanJetsAll)
    for index,j in enumerate(event.cleanJetsAll,start=0) :
      j.blike_VBF=-2
    if (nJet>=4):
      if ((event.cleanJetsAll[0].pt()>92.) and (event.cleanJetsAll[1].pt()>76.) and (event.cleanJetsAll[2].pt()>64) and (event.cleanJetsAll[3].pt()>30)):
        jetsForHiggs4=[jet for index,jet in enumerate(event.cleanJetsAll)  if (jet.jetID("POG_PFID_Loose")) and (index<4)]
        jetsForHiggs4btag=[jet.btag("pfCombinedInclusiveSecondaryVertexV2BJetTags") for index,jet in enumerate(event.cleanJetsAll) if (jet.jetID("POG_PFID_Loose")) and (index<4) ]
        btag_max1=max(jetsForHiggs4btag)
        btag_max1_idx=jetsForHiggs4btag.index(btag_max1)
        event.cleanJetsAll[0].VBF_btagMax=btag_max1
        event.cleanJetsAll[0].VBF_btagMaxIdx=btag_max1_idx
        event.cleanJetsAll[0].VBF_nJet=jetsForHiggs4btag[3]
        event.cleanJetsAll[0].VBF_4Len0=event.cleanJetsAll[0].btag("pfCombinedInclusiveSecondaryVertexV2BJetTags")
        event.cleanJetsAll[0].VBF_4Len1=event.cleanJetsAll[1].btag("pfCombinedInclusiveSecondaryVertexV2BJetTags")
        event.cleanJetsAll[0].VBF_4Len2=event.cleanJetsAll[2].btag("pfCombinedInclusiveSecondaryVertexV2BJetTags")
        event.cleanJetsAll[0].VBF_4Len3=event.cleanJetsAll[3].btag("pfCombinedInclusiveSecondaryVertexV2BJetTags")

        event.cleanJetsAll[0].VBF_4btagLen=len(jetsForHiggs4btag)
        if btag_max1>0.7:
          bjet1=jetsForHiggs4[btag_max1_idx]
          other3jets=[jet for index,jet in enumerate(jetsForHiggs4) if index!= btag_max1_idx]
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
            bbDeltaPhi=abs(deltaPhi(bjet1.phi(),bjet2.phi()))
            qqDeltaEta=abs(qjets[0].eta()-qjets[1].eta())
           
            event.cleanJetsAll[0].VBF_phi1=qjets[0].phi()
            event.cleanJetsAll[1].VBF_phi2=qjets[1].phi()
            event.cleanJetsAll[0].VBF_q1Idx=maxEta1
            event.cleanJetsAll[1].VBF_q2Idx=maxEta2

            if (Mqq>460) :
              if (qqDeltaEta>4.1) : 
                if (bbDeltaPhi<1.6) :
                  if (event.HLT_BIT_HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v) : 
                    loopMaxJet=7
                    if nJet<7 : loopMaxJet=nJet  
                    jetsForHiggsMax=[jet for index,jet in enumerate(event.cleanJetsAll) if (jet.jetID("POG_PFID_Loose")) and (index<loopMaxJet) and (jet.pt()>20)]
                    jetsForHiggsMaxEta=[abs(jet.eta()) for jet in jetsForHiggsMax]
                    jetsForHiggsMaxEtaSorted=sorted(range(len(jetsForHiggsMaxEta)),key=lambda x:jetsForHiggsMaxEta[x])
                    jetsEtaIdx=sorted(range(len(jetsForHiggsMaxEta)),key=lambda x:jetsForHiggsMaxEtaSorted[x])
                    jetsForHiggsMaxBtag=[jet.btag("pfCombinedInclusiveSecondaryVertexV2BJetTags") for jet in jetsForHiggsMax]
                    jetsForHiggsMaxBtagSorted=sorted(range(len(jetsForHiggsMaxBtag)),key=lambda x:jetsForHiggsMaxBtag[x], reverse=True)
                    jetsBtagIdx=sorted(range(len(jetsForHiggsMaxBtag)),key=lambda x:jetsForHiggsMaxBtagSorted[x])
                    
                    event.cleanJetsAll[0].VBF_loopJet=loopMaxJet
                    event.cleanJetsAll[0].VBF_JetMaxlen=len(jetsForHiggsMax)
                    for i  in range(len(jetsForHiggsMax)) :
                      j=jetsForHiggsMax[i]
                      self.Jet_pt[0]=j.pt()
                      self.Jet_eta[0]=abs(j.eta())
                      btag = j.btag("pfCombinedInclusiveSecondaryVertexV2BJetTags")
                      if btag < 0 : btag = 0
                      if btag > 1 : btag = 1
                      self.Jet_btagCSV[0]=btag
                      self.Jet_pt_idx[0]=i
                      self.Jet_eta_idx[0]=jetsEtaIdx[i]
                      self.Jet_btagCSV_idx[0]=jetsBtagIdx[i]
                      j.blike_VBF=self.reader.EvaluateMVA(self.name)
                      event.cleanJetsAll[i].VBFetaIdx=jetsEtaIdx[i]
                      event.cleanJetsAll[i].VBFbtagIdx=jetsBtagIdx[i]
                      event.cleanJetsAll[i].VBFptIdx=i
                      event.cleanJetsAll[i].VBFpt=j.pt()
                      event.cleanJetsAll[i].VBFbtag=btag
                      event.cleanJetsAll[i].VBFeta=abs(j.eta())
                  else :  event.cleanJetsAll[0].cut_debug=8
                else :  event.cleanJetsAll[0].cut_debug=7
              else : event.cleanJetsAll[0].cut_debug=6
            else : event.cleanJetsAll[0].cut_debug=5
          else : event.cleanJetsAll[0].cut_debug=4
        else : event.cleanJetsAll[0].cut_debug=3
      else : event.cleanJetsAll[0].cut_debug=2
    else : event.cleanJetsAll[0].cut_debug=1
    
