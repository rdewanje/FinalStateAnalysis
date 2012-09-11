'''

Analyze MMT events for the WH analysis

'''

from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor
import glob
from MuTauTauTree import MuTauTauTree
import os
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import WHAnalyzerBase
import ROOT

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')
highpt_mu_fr = build_roofunctor(
    frfit_dir + '/m_wjets_pt20_pfidiso02_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
lowpt_mu_fr = build_roofunctor(
    frfit_dir + '/m_wjets_pt10_pfidiso02_muonJetPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)


tau_fr = build_roofunctor(
    frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################

# Determine MC-DATA corrections
is7TeV = bool('7TeV' in os.environ['jobid'])
print "Is 7TeV:", is7TeV

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions = glob.glob(os.path.join(
    'inputs', os.environ['jobid'], 'data_DoubleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight(
    'S6' if is7TeV else 'S7', *pu_distributions)

muon_pog_PFTight_2011 = MuonPOGCorrections.make_muon_pog_PFTight_2011()
muon_pog_PFTight_2012 = MuonPOGCorrections.make_muon_pog_PFTight_2012()

muon_pog_PFRelIsoDB02_2011 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB02_2011()
muon_pog_PFRelIsoDB02_2012 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB02_2012()

muon_pog_Mu17Mu8_Mu17_2012 = MuonPOGCorrections.make_muon_pog_Mu17Mu8_Mu17_2012()
muon_pog_Mu17Mu8_Mu8_2012 = MuonPOGCorrections.make_muon_pog_Mu17Mu8_Mu8_2012()

# takes etas of muons
muon_pog_Mu17Mu8_2011 = MuonPOGCorrections.muon_pog_Mu17Mu8_eta_eta_2011

# Get object ID and trigger corrector functions
def mc_corrector_2011(row):
    if row.run > 2:
        return 1
    pu = pu_corrector(row.nTruePU)
    #pu = 1
    mid = muon_pog_PFTight_2011(row.mPt, row.mEta)
    miso = muon_pog_PFRelIsoDB02_2011(row.mPt, row.mEta)
    trigger = muon_pog_Mu17Mu8_2011(row.m1Eta, row.m2Eta)
    return pu*m1id*m2id*m1iso*m2iso*trigger

def mc_corrector_2012(row):
    if row.run > 2:
        return 1
    pu = pu_corrector(row.nTruePU)
    mid = muon_pog_PFTight_2012(row.mPt, row.mEta)
    miso = muon_pog_PFRelIsoDB02_2012(row.mPt, row.mEta)
    mTrig = muon_pog_Mu17Mu8_Mu17_2012(row.mPt, row.mEta)
    return pu*m1id*m2id*m1iso*m2iso*mTrig

# Determine which set of corrections to use
mc_corrector = mc_corrector_2011
if not is7TeV:
    mc_corrector = mc_corrector_2012

################################################################################
#### Analysis logic ############################################################
################################################################################

class WHAnalyzeMMT(WHAnalyzerBase.WHAnalyzerBase):
    tree = 'mtt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(WHAnalyzeMTT, self).__init__(tree, outfile, MuTauTauTree, **kwargs)
        # Hack to use S6 weights for the one 7TeV sample we use in 8TeV
        target = os.environ['megatarget']
        if 'HWW3l' in target:
            print "HACK using S6 PU weights for HWW3l"
            global pu_corrector
            pu_corrector =  PileupWeight.PileupWeight('S6', *pu_distributions)

    def book_histos(self, folder):
        self.book(folder, "weight", "Event weight", 100, 0, 5)
        self.book(folder, "weight_nopu", "Event weight without PU", 100, 0, 5)
        self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        self.book(folder, "prescale", "HLT prescale", 26, -5.5, 20.5)
        self.book(folder,"t1t2Pt","Pt of 2 Leading Taus/Event",100,0,200)
        self.book(folder, "mPt", "Muon  Pt", 100, 0, 100)
        self.book(folder, "mJetPt", "Muon Jet Pt", 100, 0, 200)
        self.book(folder, "mAbsEta", "Muon AbsEta", 100, 0, 2.4)
        self.book(folder, "subMass", "subleadingMass", 200, 0, 200)
        self.book(folder, "leadMass", "leadingMass", 200, 0, 200)
        # Rank muons by less MT to MET, for WZ control region
        self.book(folder, "subMTMass", "subMTMass", 200, 0, 200)
        self.book(folder, "t1Pt", "Tau 1 Pt", 100, 0, 100)
        self.book(folder, "t2Pt", "Tau 2 Pt", 100, 0, 100)
        self.book(folder, "t1AbsEta", "Tau 1 AbsEta", 100, 0, 2.3)
        self.book(folder, "t2AbsEta", "Tau 2 AbsEta", 100, 0, 2.3)
        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)

    def fill_histos(self, histos, folder, row, weight):
        histos['/'.join(folder + ('nTruePU',))].Fill(row.nTruePU)
        def fill(name, value):
            histos['/'.join(folder + (name,))].Fill(value, weight)
        histos['/'.join(folder + ('weight',))].Fill(weight)

        fill('prescale', row.doubleMuPrescale)
        fill('rho', row.rho)
        fill('nvtx', row.nvtx)
        fill('t1t2Pt',row.t1t2Pt)
        fill('mPt', row.mPt)
        fill('mJetPt', row.mJetPt)
        fill('mAbsEta', row.mAbsEta)
        fill('subMass', row.m_t2_Mass)
        fill('leadMass', row.m_t1_Mass)
        fill('t1Pt', row.t1Pt)
        fill('t2Pt', row.t2Pt)
        fill('t1AbsEta', row.t1AbsEta)
        fill('t2AbsEta', row.t2AbsEta)
        if row.m1MtToMET > row.m2MtToMET:
            fill('subMTMass', row.m2_t_Mass)
        else:
            fill('subMTMass', row.m1_t_Mass)

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if row.t1_t2_SS:
            return False
        if row.m_t1_Mass > 80 and row.t1t2Pt > 50:
            return False
        if row.t1MuOverlap:
            return False
        if row.t2MuOverlap:
            return False
        if not row.mPFIDTight:
            return False
        if row.mRelPFIsoDB > 0.1:
            return False
        if row.mPt < 24:
            return False
        if row.t1Pt < 25:
            return False
        if row.t2Pt < 20:
            return False
        if row.mAbsEta > 2.1:
            return False
        if row.t1AbsEta > 2.3:
            return False
        if row.t2AbsEta > 2.3:
            return False

        if abs(row.mDZ) > 0.2:
            return False
        if abs(row.t1DZ) > 0.2:
            return False
        if abs(row.t2DZ) > 0.2:
            return False
                    
        if row.LT < 80:
            return False

        if row.muVetoPt5:
            return False
        if row.bjetCSVVeto:
            return False
        if row.tauVetoPt20:
            return False
        if row.eVetoCicTightIso:
            return False

        if not row.mPixHits:
            return False
        
        ## Fixme use CSV
        if row.mJetBtag > 3.3:
            return False
       
        if not row.t1AntiMuonTight:
            return False
        if not row.t2AntiMuonTight:
            return False
        if not row.t1AntiElectronLoose:
            return False
        if not row.t2AntiElectronMedium:
            return False
            
        if not self.trigger_match_m(row):
            return False
        

        return True

    @staticmethod
    def trigger_match_m(row):
        return True
        if row.mDiMuonL3p5PreFiltered8  > 0 or \
           row.mDiMuonL3PreFiltered7  > 0 or \
           row.mSingleMu13L3Filtered13  > 0 or \
           row.mSingleMu13L3Filtered17  > 0 or \
           row.mDiMuonMu17Mu8DzFiltered0p2  > 0 or \
           row.mL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17:
            return True

    

    def sign_cut(self, row):
        ''' Returns true if taus are OS '''
        return bool(row.t1_t2_OS)

    def obj1_id(self, row):
        return bool(row.mPFIDTight) and bool(row.mRelPFIsoDB < 0.2)

    def obj2_id(self, row):
        return bool(row.t1TightMVAIso)
    
    def obj3_id(self, row):
        return bool(row.t2TightMVAIso)


    def anti_wz(self, row):
        return row.t1AntiMuonTight and not row.t1MuOverlap and row.t2AntiMuonTight and not row.t2MuOverlap

    def enhance_wz(self, row):
        # Require the "tau" to be a muon, and require the third muon
        # to have M_Z +- 20
        if row.tAntiMuonTight or not row.tMuOverlap:
            return False
        # Cut on m2 PT > 20
        #if row.m2Pt < 20:
            #return False
        # Make sure any Z is from m1
        m2_good_Z = bool(71 < row.m2_t_Mass < 111)
        return not m2_good_Z

    def event_weight(self, row):
        return mc_corrector(row)

    def obj1_weight(self, row):
        return highpt_mu_fr(max(row.mJetPt, row.mPt))
        #return highpt_mu_fr(row.mPt)

    def obj2_weight(self, row):
        return tau_fr(row.t1Pt)

    def obj3_weight(self, row):
        return tau_fr(row.t2Pt)


    # For measuring charge flip probability
    # Not really used in this channel
    def obj1_obj3_SS(self, row):
        return not row.m_t1_SS

    def obj1_charge_flip(self, row):
        return 0
