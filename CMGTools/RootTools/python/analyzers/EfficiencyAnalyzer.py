from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.statistics.Average import Average
from CMGTools.RootTools.statistics.Counter import Counter
from CMGTools.RootTools.statistics.Histograms import Histograms
from CMGTools.RootTools.physicsobjects.PhysicsObjects import PhysicsObject, GenParticle
from CMGTools.RootTools.utils.DeltaR import cleanObjectCollection, matchObjectCollection
from CMGTools.RootTools.utils.OOTPileUpReweighting import ootPUReweighter

from ROOT import TH1F, TH2F, TFile


class PileUpSummaryInfo( object ):
    def __init__(self, object ):
        self.object = object
        
    def __getattr__(self,name):
        '''all accessors  from cmg::DiTau are transferred to this class.'''
        return getattr(self.physObj, name)

    def nPU(self):
        return self.object.getPU_NumInteractions()
    
    def __str__(self):
        tmp = '{className} : bunchx = {bunchx}; numPU = {numpu}'.format(
            className = self.__class__.__name__,
            bunchx = self.object.getBunchCrossing(),
            numpu = self.object.getPU_NumInteractions() )
        return tmp


class PhaseSpace( object ):
    def __init__(self, name, etaRanges, ptRanges=None):
        self.name = name
        self.etaRanges = etaRanges
        self.ptRanges = ptRanges

    def etaPass(self, eta):
        for range in self.etaRanges:
            if eta>range[0] and eta<range[1]:
                return True
        return False

ECAL_barrel = PhaseSpace('EB',
                         [(-1.4, 1.4)])
ECAL = PhaseSpace('ECAL',
                  [(-2.9, 2.9)])
ECAL_endcaps = PhaseSpace('EE',
                          [(-2.9, -1.6),
                           (1.6, 2.9) ])

muon_barrel = PhaseSpace('MB',
                         [(-0.9, 0.9)])
muon_transition = PhaseSpace('MT',
                             [(-1.2, -0.9),
                              (0.9, 1.2) ])
muon_endcaps = PhaseSpace('ME',
                          [ (-2.5, -1.2),
                            (1.2, 2.5) ])
muon = PhaseSpace('MUON',
                  [(-2.5, 2.5)])

electronPhaseSpaces = [ECAL_barrel, ECAL_endcaps, ECAL]
muonPhaseSpaces = [muon_barrel, muon_transition, muon_endcaps, muon]


class EfficiencyHistograms( Histograms ):
    '''
    Define variable binning in pT
    Define detector boundaries in a generic way
    Efficiency plotter
    '''
    def __init__(self, name, space):
        self.h_pt = TH1F(name + '_h_pt', ';p_{T} (GeV);Efficiency', 50, 0, 200)
        self.h_eta = TH1F(name + '_h_eta', ';#eta;Efficiency', 50, -6, 6)
        self.h_phi = TH1F(name + '_h_phi', ';#phi (rad);Efficiency', 50, -6.3, 6.3)
        self.h_pv = TH1F(name + '_h_pv', ';# rec vertices', 40, 0, 40)
        self.h_pu = TH1F(name + '_h_pu', ';# PU interactions', 40, 0, 40)
        self.h_pup_VS_pu = TH2F(name + '_h_pup_VS_pu',
                                ';# PU interactions;# PU interactions in next bunch',
                                50, 0, 50, 50, 0, 100)
        self.space = space 
        super( EfficiencyHistograms, self).__init__(name)
 
    def ptPass(self, pt):
        ptMin = 10.
        return pt>ptMin  

    def fillParticle(self, particle, event, weight):
        eta = particle.eta()
        pt = particle.pt()

        #COLIN not good! should test only the gen particle.
        
        # temp: 
        if not self.space.etaPass(eta) or \
           not self.ptPass(pt):
            return
        
        if self.space.etaPass(eta):
            self.h_pt.Fill(pt, weight)
        if self.ptPass(pt):
                self.h_eta.Fill( eta, weight)
        if self.space.etaPass(eta) and \
           self.ptPass(pt):
            self.h_phi.Fill( particle.phi(), weight)
            self.h_pv.Fill( len(event.vertices), weight)
            self.h_pu.Fill( event.pusi[1].nPU(), weight)
            self.h_pup_VS_pu.Fill( event.pusi[1].nPU(), event.pusi[2].nPU(), weight)
            
    def fillParticles(self, particles, event, weight):
        for particle in particles:
            self.fillParticle( particle, event, weight)

        


class EfficiencyAnalyzer( Analyzer ):
    '''A simple jet analyzer for Pietro.'''

    
    def declareHandles(self):
        super(EfficiencyAnalyzer, self).declareHandles()

        instance = self.cfg_ana.instance
        type = self.cfg_ana.type
        self.handles['rec'] =  AutoHandle(
            instance,
            type 
            )

        self.handles['vertices'] =  AutoHandle(
            'offlinePrimaryVertices',
            'std::vector<reco::Vertex>'
            )

        geninstance = 'genLeptonsStatus1'
        gentype = 'std::vector<reco::GenParticle>'
        self.mchandles['gen'] =  AutoHandle(
            geninstance,
            gentype 
            )

        self.mchandles['pusi'] =  AutoHandle(
            'addPileupInfo',
            'std::vector<PileupSummaryInfo>' 
            )        
  
    def beginLoop(self):
        super(EfficiencyAnalyzer,self).beginLoop()
        self.file = TFile( '/'.join( [self.dirName, 'EfficiencyAnalyzer.root']),
                           'recreate')

        print self.cfg_ana
        self.phaseSpaces = None
        if self.cfg_ana.genPdgId==13:
            self.phaseSpaces = muonPhaseSpaces
        elif self.cfg_ana.genPdgId==11:
            self.phaseSpaces = electronPhaseSpaces

        for space in self.phaseSpaces:
            space.denomHistos = EfficiencyHistograms('_'.join([space.name,
                                                               'Denom']),
                                                     space )
            space.numHistos = EfficiencyHistograms('_'.join([space.name,
                                                             'Num']),
                                                   space)
            space.counter = Counter( space.name )
            
        # self.denomHistos = EfficiencyHistograms('Denom')
        # self.numHistos = EfficiencyHistograms('Num')
        # self.counters.addCounter('effcount') 


    def process(self, iEvent, event):
        self.readCollections( iEvent )

        event.pusi = map( PileUpSummaryInfo, self.mchandles['pusi'].product() )  
        event.rec = self.handles['rec'].product()
        event.gen = self.mchandles['gen'].product()
        event.vertices = self.handles['vertices'].product()

        refselFun = None
        if hasattr( self.cfg_ana, 'refselFun'):
            refselFun = self.cfg_ana.refselFun
        if refselFun is not None:
            event.refsel = [ PhysicsObject(obj) for obj in event.rec if refselFun(obj)]
        else:
            event.refsel = event.gen

        recselFun = None
        if hasattr( self.cfg_ana, 'recselFun'):
            recselFun = self.cfg_ana.recselFun
        if recselFun is not None:
            event.recsel = [ PhysicsObject(obj) for obj in event.rec if recselFun(obj)]
        else:
            event.recsel = event.rec

        genpdgid = self.cfg_ana.genPdgId
        event.gensel = [ GenParticle(obj) for obj in event.gen if abs(obj.pdgId())==genpdgid]

        # gen objects matched to a reference lepton
        pairs = matchObjectCollection( event.gensel, event.refsel, 0.1)
        event.genmatchedRef = [ gen for gen,rec in pairs.iteritems() if rec is not None]

        # and gen objects wich are in addition matched to a
        # selected lepton
        pairs = matchObjectCollection( event.genmatchedRef, event.recsel, 0.1)
        event.genmatched = [ gen for gen,rec in pairs.iteritems() if rec is not None]

        # reweightign OOTPU in chamonix samples to the OOTPU observed in Fall11 samples
        weight = 1
        if self.cfg_comp.name.find('Chamonix')!=-1:
            weight = ootPUReweighter.getWeight( event.pusi[1].nPU(), event.pusi[2].nPU())
        
        for space in self.phaseSpaces:        
            space.denomHistos.fillParticles( event.genmatchedRef, event, weight)
            space.numHistos.fillParticles( event.genmatched, event, weight)
            # space.counter.inc('passed')
        
        
    def write(self):
        for space in self.phaseSpaces:
            space.denomHistos.Write( self.file )
            space.numHistos.Write( self.file )
