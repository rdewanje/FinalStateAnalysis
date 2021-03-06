'''

Definition of 7TeV samples and data.

NB the real data samples are built automatically at the bottom.

Author: Evan K. Friis, UW Madison

'''

from datacommon import square, cube, quad, picobarns, br_w_leptons
import re

# Define a mapping between a "nice" name and a set of datasets.
# This is so you can make an update to the underlying data sample pythia/powheg
# etc.
data_name_map = {
    'Zjets' : ['Zjets_M50',],

    'QCDMu' : ['QCD_20toInf_MuPt15'],

    'Wjets' : ['WplusJets_madgraph'],

    'WW' : ['WWJetsTo2L2Nu'],
    'WZ' : ['WZJetsTo3LNu'],
    'WZ_pythia' : ['WZJetsTo3LNu_pythia'],
    'ZZ' : ['ZZJetsTo4L_pythia'],

    'ttjets': ['TTplusJets_madgraph'],

    'VGamma': ['VGjets'],

    'VH100' : ['VH_100'],
    'VH110' : ['VH_110'],
    'VH115' : ['VH_115'],
    'VH120' : ['VH_120'],
    'VH125' : ['VH_125'],
    'VH130' : ['VH_130'],
    'VH135' : ['VH_135'],
    'VH140' : ['VH_140'],
    'VH145' : ['VH_145'],
    'VH150' : ['VH_150'],
    'VH160' : ['VH_160'],

    'VH110WW' : ['WH_110_HWW3l'],
    'VH115WW' : ['WH_115_HWW3l'],
    'VH120WW' : ['WH_120_HWW3l'],
    'VH125WW' : ['WH_125_HWW3l'],
    'VH130WW' : ['WH_130_HWW3l'],
    'VH135WW' : ['WH_135_HWW3l'],
    'VH140WW' : ['WH_140_HWW3l'],
    'VH145WW' : ['WH_145_HWW3l'],
    'VH150WW' : ['WH_150_HWW3l'],
    'VH155WW' : ['WH_155_HWW3l'],
    'VH160WW' : ['WH_160_HWW3l'],

    'TTW' : ['TTWToLplus', 'TTWToLminus'],
    'TTZ' : ['TTZToLplus', 'TTZToLminus'],
    'WWW' : ['WWWTo2Lplus', 'WWWTo2Lminus'],

    'ggH_ZZ_4l_120' : ['ggH_ZZ_4l_120'],
}

datadefs = {
    ############################################################################
    #### EWK background datasets            ####################################
    ############################################################################

    'WbbToLNu_TuneZ2_7TeV-madgraph-pythia6-tauola' : {
        'datasetpath' : '/WbbToLNu_TuneZ2_7TeV-madgraph-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : 999, #NNLO
        'pu' : 'S6',
        'analyses' : ['Wbb',  'VH', 'Mu'],
        'responsible' : 'Tapas',
    },
    'Zjets_M50' : {
        'datasetpath' : '/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : 3048*picobarns, #NNLO
        'pu' : 'S6',
        'analyses' : ['HTT',  'VH', 'Tau', 'Mu'],
        'responsible' : 'Maria',
    },
    'WplusJets_madgraph' : {
        'datasetpath' : "/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 31314*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT',  'VH', 'Tau', 'Mu'],
        'responsible' : 'Maria',
    },
    'TTplusJets_madgraph' : {
        'datasetpath' : "/TTJets_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 157.5*picobarns, # NLO cross-section from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
        'pu' : 'S6',
        'analyses' : ['HTT',  'VH', 'Tau', 'Mu'],
        'responsible' : 'Ian',
    },
    # Single top samples
    'T_tW_Powheg' : {
        'datasetpath' : "/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 7.87*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Maria',
    },
    'T_t_Powheg' : {
        'datasetpath' : "/T_TuneZ2_t-channel_7TeV-powheg-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 41.92*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Maria',
    },
    'T_s_Powheg' : {
        'datasetpath' : "/T_TuneZ2_s-channel_7TeV-powheg-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 3.19*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Maria',
    },
    # Single anti-top samples
    'Tbar_tW_Powheg' : {
        'datasetpath' : "/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/Fall11-PU_S6_START42_V14B-v2/AODSIM",
        'x_sec' : 7.87*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Evan',
    },
    'Tbar_t_Powheg' : {
        'datasetpath' : "/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 22.65*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Evan',
    },
    'Tbar_s_Powheg' : {
        'datasetpath' : "/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1.44*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Evan',
    },
    'Wplus3Jets_madgraph' : {
        'datasetpath' : "/W3Jets_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 304.2*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Maria',
    },


    ############################################################################
    #### VGamma background datasets         ####################################
    ############################################################################
    'VGjets' : {
        'datasetpath' : '/GVJets_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : 56.64*picobarns,
        'pu' : 'S6',
        'analyses' : ['Mu', 'VH'],
        'responsible' : 'Evan',
    },

    ############################################################################
    #### Diboson datasets                   ####################################
    ############################################################################

    'ZZJetsTo4L_pythia' : {
        'datasetpath' : "/ZZTo4L_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.106*picobarns, # from MCFM via Ian
        'x_sec_unc' : quad(1.5, 0.2, 0.2)*0.10096*0.10096,
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'WZJetsTo3LNu' : {
        'datasetpath' : "/WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : 26.735*picobarns*3*0.03365*(0.1075+0.1057+0.1125) ,
        'x_sec_unc' : quad(2.4, 1.1, 1.0)*0.3257*0.10096,
        'analyses' : ['VH', ],
        'responsible' : 'Evan',
    },
    'WZJetsTo2L2Q' : {
        'datasetpath' : "/WZJetsTo2L2Q_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : 3.85*picobarns ,
        'analyses' : ['HTT'],
        'responsible' : 'Ian',
    },
    'WZJetsTo3LNu_pythia' : {
        'datasetpath' : "/WZTo3LNu_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        # This xsec comes from PREP and is only LO.  Different from the PREP
        # value for the madgraph sample, as Pythia does not include gamma*
        'x_sec' : 0.33*picobarns, # FIXME !!!!!!!!
        'x_sec_unc' : quad(2.4, 1.1, 1.0)*0.3257*0.10096,
        'analyses' : ['VH',  'HTT'],
        'responsible' : 'Evan',
    },
    'WWJetsTo2L2Nu' : {
        'datasetpath' : '/WWJetsTo2L2Nu_TuneZ2_7TeV-madgraph-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'pu' : 'S6',
        #'x_sec' : 3.783*picobarns, # FROM PREP
        # 55.3 +- 3.3 6.9 3.3 from EWK-11-10
        'x_sec' : 55.3*picobarns*0.3257*0.3257, # 32.57% BR to leptons
        'x_sec_unc' : quad(3.3, 6.9, 3.3)*0.3257*0.3257,
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Ian',
    },
    'WZinclusive' : {
        'datasetpath' : "/WZ_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : 18*picobarns,
        'analyses' : ['VH',  'HTT'],
        'responsible' : 'Evan',
    },
    'WWinclusive' : {
        'datasetpath' : '/WW_TuneZ2_7TeV_pythia6_tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'pu' : 'S6',
        'x_sec' : 47*picobarns,
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'GluGluToZZTo4L' : {
        'datasetpath' : "/GluGluToZZTo4L_7TeV-gg2zz-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.00174*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'GluGluToZZTo2L2L' : {
        'datasetpath' : "/GluGluToZZTo2L2L_7TeV-gg2zz-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.00348*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'ZZTo4mu_powheg_v2' : { #v1 samples had a bug where the Z->4l peak weighting screwed the overall cross-section IAR 31.May.2012
        'datasetpath' : "/ZZTo4mu_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.03067*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'ZZTo4e_powheg_v2' : { #v1 samples had a bug where the Z->4l peak weighting screwed the overall cross-section IAR 31.May.2012
        'datasetpath' : "/ZZTo4e_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.03067*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'ZZTo4tau_powheg_v2' : { #v1 samples had a bug where the Z->4l peak weighting screwed the overall cross-section IAR 31.May.2012
        'datasetpath' : "/ZZTo4tau_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.03067*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'ZZTo2e2tau_powheg_v2' : { #v1 samples had a bug where the Z->4l peak weighting screwed the overall cross-section IAR 31.May.2012
        'datasetpath' : "/ZZTo2e2tau_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.01386*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'ZZTo2e2mu_powheg_v2' : { #v1 samples had a bug where the Z->4l peak weighting screwed the overall cross-section IAR 31.May.2012
        'datasetpath' : "/ZZTo2e2mu_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.01386*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },
    'ZZTo2mu2tau_powheg_v2' : { #v1 samples had a bug where the Z->4l peak weighting screwed the overall cross-section IAR 31.May.2012
        'datasetpath' : "/ZZTo2mu2tau_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v2/AODSIM",
        'pu' : 'S6',
        'x_sec' : 0.01386*picobarns, # from https://twiki.cern.ch/twiki/bin/view/CMS/HZZSamples7TeV
        'analyses' : ['VH',  '4L', 'HTT'],
        'responsible' : 'Ian',
    },

    ############################################################################
    #### QCD datasets                       ####################################
    ############################################################################

    'QCD_20toInf_MuPt15' : {
        'datasetpath' : '/QCD_Pt-20_MuEnrichedPt-15_TuneZ2_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'pu' : 'S6',
        'x_sec' : 2.966E8*picobarns*2.855E-4,
        'analyses' : ['HTT', 'VH', 'Tau', 'Mu'],
        'responsible' : 'Josh',
    },

    ############################################################################
    #### Signal datasets                    ####################################
    ############################################################################

    'VH_100' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-100_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : (1.186 + 0.6313 + 0.1638)*picobarns*8.36e-2,
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_110' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-110_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : (0.8754 + 0.4721 + 0.1257)*picobarns*8.02e-2,
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_115' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-115_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : (0.7546 + 0.4107 + 0.1106)*picobarns*7.65e-2,
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_120' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-120_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'pu' : 'S6',
        'x_sec' : (0.6561 + 0.3598 + 0.09756)*picobarns*7.1e-2,
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_125' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-125_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.6729 + 0.3158 + 0.08634)*picobarns*6.37e-2,
        'pu' : 'S6',
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_130' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-130_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.5008 + 0.2778 + 0.07658)*picobarns*5.48e-2,
        'pu' : 'S6',
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_135' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-135_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.4390 + 0.2453 + 0.06810)*picobarns*4.52e-2,
        'pu' : 'S6',
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_140' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-140_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.3857 + 0.2172 + 0.06072)*picobarns*3.54e-2,
        'pu' : 'S6',
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_145' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-145_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.3406 + 0.1930 + 0.05435)*picobarns*2.61e-2,
        'pu' : 'S6',
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_150' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-150_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.3001 + 0.1713 + 0.04869)*picobarns*1.78e-2,
        'pu' : 'S6',
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },
    'VH_160' : {
        'datasetpath' :"/WH_ZH_TTH_HToTauTau_M-160_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.2291 + 0.1334 + 0.03942)*picobarns*3.96E-03,
        'pu' : 'S6',
        'analyses' : ['VH', 'HTT'],
        'responsible' : 'Evan',
    },

    ############################################################################
    #### VH -> WW dataset                   ####################################
    ############################################################################

    'VH_120_HWW' : {
        'datasetpath' :"/WH_ZH_TTH_HToWW_M-120_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.6561 + 0.3598 + 0.09756)*picobarns*1.43E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'VH_130_HWW' : {
        'datasetpath' :"/WH_ZH_TTH_HToWW_M-130_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.5008 + 0.2778 + 0.07658)*picobarns*3.05E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'VH_135_HWW' : {
        'datasetpath' :"/WH_ZH_TTH_HToWW_M-135_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM	",
        'x_sec' : (0.4390 + 0.2453 + 0.06810)*picobarns*4.03E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'VH_140_HWW' : {
        'datasetpath' :"/WH_ZH_TTH_HToWW_M-140_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.3857 + 0.2172 + 0.06072)*picobarns*5.03E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'VH_150_HWW' : {
        'datasetpath' :"/WH_ZH_TTH_HToWW_M-150_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM	",
        'x_sec' : (0.3001 + 0.1713 + 0.04869)*picobarns*6.98E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'VH_160_HWW' : {
        'datasetpath' :"/WH_ZH_TTH_HToWW_M-160_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM	",
        'x_sec' : (0.2291 + 0.1334 + 0.03942)*picobarns*9.08E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'WH_110_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-110_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.8754*picobarns*cube(br_w_leptons)*4.82E-02,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_115_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-115_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.7546*picobarns*cube(br_w_leptons)*8.67E-02,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_120_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-120_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.6561*picobarns*cube(br_w_leptons)*1.43E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_125_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-125_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.5729*picobarns*cube(br_w_leptons)*2.16E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_130_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-130_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.5008*picobarns*cube(br_w_leptons)*3.05E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_135_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-135_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.4390*picobarns*cube(br_w_leptons)*4.03E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_140_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-140_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.3857*picobarns*cube(br_w_leptons)*5.03E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_145_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-145_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.3406*picobarns*cube(br_w_leptons)*6.02E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_150_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-150_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.3001*picobarns*cube(br_w_leptons)*6.98E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_155_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-155_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.2646*picobarns*cube(br_w_leptons)*7.95E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },
    'WH_160_HWW3l' : {
        'datasetpath' : "/WH_HToWW_3l_M-160_7TeV-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 0.2291*picobarns*cube(br_w_leptons)*9.08E-01,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    ############################################################################
    #### ZH private production              ####################################
    ############################################################################

    #'ZH_110_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-110-7TeV-pythia6/friis-RECO_crab_reco_TT_110-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.4721*picobarns*8.02e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_115_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-115-7TeV-pythia6/friis-RECO_crab_reco_TT_115-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.4107*picobarns*7.65e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_120_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-120-7TeV-pythia6/friis-RECO_crab_reco_TT_120-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.3598*picobarns*7.1e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_125_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-125-7TeV-pythia6/friis-RECO_crab_reco_TT_125-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.3158*picobarns*6.37e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_130_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-130-7TeV-pythia6/friis-RECO_crab_reco_TT_130-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.2778*picobarns*5.48e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_135_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-135-7TeV-pythia6/friis-RECO_crab_reco_TT_135-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.2453*picobarns*4.52e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_140_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-140-7TeV-pythia6/friis-RECO_crab_reco_TT_140-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.2172*picobarns*3.54e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_145_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-145-7TeV-pythia6/friis-RECO_crab_reco_TT_145-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.1930*picobarns*2.61e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_150_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-150-7TeV-pythia6/friis-RECO_crab_reco_TT_150-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.1713*picobarns*1.78e-2,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_155_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-155-7TeV-pythia6/friis-RECO_crab_reco_TT_155-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_160_HTT' : {
        #'datasetpath' : "/ZH_ZToLL_HToTauTau_M-160-7TeV-pythia6/friis-RECO_crab_reco_TT_160-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 0.1334*picobarns*3.96e-3,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},


    #'ZH_110_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-110-7TeV-pythia6/friis-RECO_crab_reco_WW_110-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_115_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-115-7TeV-pythia6/friis-RECO_crab_reco_WW_115-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_120_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-120-7TeV-pythia6/friis-RECO_crab_reco_WW_120-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_125_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-125-7TeV-pythia6/friis-RECO_crab_reco_WW_125-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_130_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-130-7TeV-pythia6/friis-RECO_crab_reco_WW_130-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_135_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-135-7TeV-pythia6/friis-RECO_crab_reco_WW_135-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_140_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-140-7TeV-pythia6/friis-RECO_crab_reco_WW_140-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_145_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-145-7TeV-pythia6/friis-RECO_crab_reco_WW_145-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_150_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-150-7TeV-pythia6/friis-RECO_crab_reco_WW_150-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_155_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-155-7TeV-pythia6/friis-RECO_crab_reco_WW_155-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    #'ZH_160_HWW' : {
        #'datasetpath' : "/ZH_ZToLL_HToWW_WWTo2LNu-160-7TeV-pythia6/friis-RECO_crab_reco_WW_160-01715716b3165466edf30580d661ec8b/USER",
        #'x_sec' : 999,
        #'pu' : 'S6',
        #'analyses' : ['VH'],
        #'dbs' : "cms_dbs_ph_analysis_01",
    #},

    ############################################################################
    #### H2Tau samples                      ####################################
    ############################################################################

    ######################### VBF H-->tau+tau ##################################
    'VBF_H2Tau_M-100' : {
        'datasetpath' : '/VBF_HToTauTau_M-100_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-105' : {
        'datasetpath' : '/VBF_HToTauTau_M-105_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-110' : {
        'datasetpath' : '/VBF_HToTauTau_M-110_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-115' : {
        'datasetpath' : '/VBF_HToTauTau_M-115_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-120' : {
        'datasetpath' : '/VBF_HToTauTau_M-120_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-125' : {
        'datasetpath' : '/VBF_HToTauTau_M-125_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-130' : {
        'datasetpath' : '/VBF_HToTauTau_M-130_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-135' : {
        'datasetpath' : '/VBF_HToTauTau_M-135_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-140' : {
        'datasetpath' : '/VBF_HToTauTau_M-140_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-145' : {
        'datasetpath' : '/VBF_HToTauTau_M-145_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'VBF_H2Tau_M-150' : {
        'datasetpath' : '/VBF_HToTauTau_M-150_7TeV-powheg-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },

    ######################### GGH SM H-->tau+tau ##################################

    'GGH_H2Tau_M-100' : {
        'datasetpath' : '/GluGluToHToTauTau_M-100_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-105' : {
        'datasetpath' : '/GluGluToHToTauTau_M-105_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-110' : {
        'datasetpath' : '/GluGluToHToTauTau_M-110_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-115' : {
        'datasetpath' : '/GluGluToHToTauTau_M-115_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-120' : {
        'datasetpath' : '/GluGluToHToTauTau_M-120_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-125' : {
        'datasetpath' : '/GluGluToHToTauTau_M-125_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-130' : {
        'datasetpath' : '/GluGluToHToTauTau_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-135' : {
        'datasetpath' : '/GluGluToHToTauTau_M-135_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-140' : {
        'datasetpath' : '/GluGluToHToTauTau_M-140_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-145' : {
        'datasetpath' : '/GluGluToHToTauTau_M-145_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGH_H2Tau_M-150' : {
        'datasetpath' : '/GluGluToHToTauTau_M-150_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },

    ######################### GGH MSSM H-->tau+tau ##################################

    'GGHMSSM_H2Tau_M-90' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-90_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-100' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-100_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-120' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-120_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-130' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-130_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-140' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-140_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-160' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-160_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-180' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-180_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-200' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-200_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-250' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-250_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-300' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-300_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-350' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-350_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-400' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-400_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-450' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-450_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-500' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-500_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-600' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-600_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-700' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-700_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-800' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-800_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-900' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-900_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'GGHMSSM_H2Tau_M-1000' : {
        'datasetpath' : '/SUSYGluGluToHToTauTau_M-1000_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },

	######################### BBH MSSM H-->tau+tau ##################################

    'BBH_H2Tau_M-90' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-90_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-100' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-100_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-120' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-120_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-130' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-130_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-140' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-140_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-160' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-160_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-180' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-180_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-200' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-200_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-250' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-250_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-300' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-300_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-350' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-350_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-400' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-400_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-450' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-450_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-500' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-500_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-600' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-600_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-700' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-700_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-800' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-800_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-900' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-900_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },
    'BBH_H2Tau_M-1000' : {
        'datasetpath' : '/SUSYBBHToTauTau_M-1000_7TeV-pythia6-tauola/Fall11-PU_S6_START42_V14B-v1/AODSIM',
        'x_sec' : -999*picobarns,
        'pu' : 'S6',
        'analyses' : ['HTT'],
        'responsible' : 'Josh',
    },

    ############################################################################
    #### Obscure VH backgrounds             ####################################
    ############################################################################

    'TTWToLplus' : {
        'datasetpath' :"/TTWTo2Lplus2Nu_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.006841)*picobarns,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'TTZToLplus' : {
        'datasetpath' :"/TTZTo2Lplus2Nu_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.002024)*picobarns,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'WWWTo2Lplus' : {
        'datasetpath' :"/WWWTo2Lplus2Nu_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        # These have some problems with the xsec for this sample
        #'x_sec' : (0.008957)*picobarns,
        # Just basically turn it off.
        'x_sec' : (0.017)*picobarns*cube(br_w_leptons),
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'TTWToLminus' : {
        'datasetpath' :"/TTWTo2Lminus2Nu_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.002705)*picobarns,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'TTZToLminus' : {
        'datasetpath' :"/TTZTo2Lminus2Nu_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : (0.001946)*picobarns,
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    'WWWTo2Lminus' : {
        'datasetpath' :"/WWWTo2Lminus2Nu_7TeV-madgraph/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        #'x_sec' : (0.004109)*picobarns,
        # These have some problems with the xsec for this sample
        # Just basically turn it off.
        'x_sec' : (0.017)*picobarns*cube(br_w_leptons),
        'pu' : 'S6',
        'analyses' : ['VH'],
        'responsible' : 'Evan',
    },

    ############################################################################
    #### UF ZZ private production (M_ll>4)  ####################################
    ############################################################################

    'uf_zz2e2m' : {
        'datasetpath' : "/Summer11/zz2e2m_powheg_GENSIMRECO_v2/USER",
        'x_sec' : (0.1525)*picobarns,
		'pu' : 'S6', #todo: check this
        'analyses' : ['4L'],
        'dbs' : "cms_dbs_ph_analysis_02",
        'responsible' : 'Ian',
    },
    'uf_zz2e2t' : {
        'datasetpath' : "/Summer11/zz2e2tau_powheg_GENSIMRECO_v2/USER",
        'x_sec' : (0.1523)*picobarns,
		'pu' : 'S6', #todo: check this
        'analyses' : ['4L'],
        'dbs' : "cms_dbs_ph_analysis_02",
        'responsible' : 'Ian',
    },
    'uf_zz2m2t' : {
        'datasetpath' : "/Summer11/zz2mu2tau_powheg_GENSIMRECO_v2/USER",
        'x_sec' : (0.1517)*picobarns,
		'pu' : 'S6', #todo: check this
        'analyses' : ['4L'],
        'dbs' : "cms_dbs_ph_analysis_02",
        'responsible' : 'Ian',
    },
    'uf_zz4e' : {
        'datasetpath' : "/Summer11/zz4e_powheg_GENSIMRECO_v2/USER",
        'x_sec' : (0.0664)*picobarns,
		'pu' : 'S6', #todo: check this
        'analyses' : ['4L'],
        'dbs' : "cms_dbs_ph_analysis_02",
        'responsible' : 'Ian',
    },
    'uf_zz4m' : {
        'datasetpath' : "/Summer11/zz4mu_powheg_GENSIMRECO_v2/USER",
        'x_sec' : (0.0661)*picobarns,
		'pu' : 'S6', #todo: check this
        'analyses' : ['4L'],
        'dbs' : "cms_dbs_ph_analysis_02",
        'responsible' : 'Ian',
    },
    'uf_zz4t' : {
        'datasetpath' : "/Summer11/zz4tau_powheg_GENSIMRECO_v2/USER",
        'x_sec' : (0.0659)*picobarns,
		'pu' : 'S6', #todo: check this
        'analyses' : ['4L'],
        'dbs' : "cms_dbs_ph_analysis_02",
        'responsible' : 'Ian',

    },
    ############################################################################
    #### H->ZZ->4L signal samples           ####################################
    ############################################################################
    #todo: add lots of HZZ samples IAR 31.May.2012
    'ggH_ZZ_4l_115' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-115_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_120' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-120_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_130' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_140' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-140_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_150' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-150_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_160' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-160_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_170' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-170_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_180' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-180_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_190' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-190_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_200' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-200_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_210' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-210_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_220' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-220_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_230' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-230_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_250' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-250_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_275' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-275_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_300' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-30_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_325' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-325_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_350' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-350_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_375' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-375_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_400' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-400_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_425' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-425_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_450' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-450_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_475' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-475_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_500' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-500_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_525' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-525_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_550' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-550_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_575' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-575_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'ggH_ZZ_4l_600' : {
        'datasetpath' : "/GluGluToHToZZTo4L_M-600_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_115' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-115_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_120' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-120_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_130' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_140' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-140_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_150' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-150_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_160' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-160_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_170' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-170_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_180' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-180_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_190' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-190_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_200' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-200_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_210' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-210_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_220' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-220_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_230' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-230_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_250' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-250_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_275' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-275_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_300' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-30_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_325' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-325_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_350' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-350_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_375' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-375_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_400' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-400_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_425' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-425_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_450' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-450_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_475' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-475_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_500' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-500_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_525' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-525_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_550' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-550_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_575' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-575_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },
    'VBF_ZZ_4l_600' : {
        'datasetpath' : "/VBF_ToHToZZTo4L_M-600_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM",
        'x_sec' : 1,
        'pu' : 'S6',
        'analyses' : ['4L'],
        'responsible' : 'Ian',
    },

}

# Add HToBB
for mass in range(100, 150, 5):
  ver=4
  if mass>135 :
    ver=1
  datadefs['WH_WToLNu_HToBB_M-%i' % mass]= {
    'datasetpath' :'/WH_WToLNu_HToBB_M-%i_7TeV-powheg-herwigpp/Fall11-PU_S6_START42_V14B-v%i/AODSIM' % (mass, ver),
    'pu' : 'S6',
    'x_sec' : -999,
    'analyses' : ['VH', 'HBB'],
    'responsible' : 'Tapas',
    }
  

# Add all the datasets
# Following https://twiki.cern.ch/twiki/bin/viewauth/CMS/Collisions2011Analysis
def build_data_set(pd, analyses, who):
  subsample_dict = {
      'data_%s_Run2011B_PromptReco_v1' % pd : {
          'datasetpath' : "/%s/Run2011B-PromptReco-v1/AOD" % pd,
          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt",
          'firstRun' : 175832,
          'lastRun' : 180296,
          'analyses' : analyses,
          'responsible' : who,
      },
      'data_%s_Run2011A_PromptReco_v6_1409' % pd : {
          'datasetpath' : "/%s/Run2011A-PromptReco-v6/AOD" % pd,
          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt",
          'firstRun' : 172620,
          'lastRun' : 175770,
          'analyses' : analyses,
          'responsible' : who,
      },
      'data_%s_Run2011A_05Aug2011_v1' % pd : {
          'datasetpath' : "/%s/Run2011A-05Aug2011-v1/AOD" % pd,
          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt",
          'firstRun' : 170053,
          'lastRun' : 172619,
          'analyses' : analyses,
          'responsible' : who,
      },
      'data_%s_Run2011A_PromptReco_v4' % pd : {
          'datasetpath' : "/%s/Run2011A-PromptReco-v4/AOD" % pd,
          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt",
          'firstRun' : 165071,
          'lastRun' : 168437,
          'analyses' : analyses,
          'responsible' : who,
      },
      'data_%s_Run2011A_May10ReReco_v1' % pd : {
          'datasetpath' : "/%s/Run2011A-May10ReReco-v1/AOD" % pd,
          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt",
          'firstRun' : 160404,
          'lastRun' : 163869,
          'analyses' : analyses,
          'responsible' : who,
      },
      'data_%s_Run2011A_16Jan2012_v1' % pd : {
          'datasetpath' : "/%s/Run2011A-16Jan2012-v1/AOD" % pd,
          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt",
          'firstRun' : 160431,
          'lastRun' : 180252,
          'analyses' : analyses,
          'responsible' : who,
      },
      'data_%s_Run2011A_16Jan2012_v1' % pd : {
          'datasetpath' : "/%s/Run2011A-16Jan2012-v1/AOD" % pd,
          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_v2.txt",
          'firstRun' : 160431,
          'lastRun' : 180252,
          'analyses' : analyses,
          'responsible' : who,
      },
# Used for 2011 ZZ results, but now the 16JanReReco is default
#      'data_%s_Run2011A_Jul05ReReco_v1' % pd : {
#          'datasetpath' : "/%s/Run2011A-05Jul2011ReReco-ECAL-v1/AOD" % pd,
#          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Jul05JSON.txt",
#          'firstRun' : 160404,
#          'lastRun' : 167913,
#          'analyses' : analyses,
#          'responsible' : who,
#      },
#      'data_%s_Run2011A_Oct03ReReco_v1' % pd : {
#          'datasetpath' : "/%s/Run2011A-03Oct2011-v1/AOD" % pd,
#          'lumi_mask' : "FinalStateAnalysis/RecoTools/data/masks/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt",
#          'firstRun' : 172620,
#          'lastRun' : 175770,
#          'analyses' : analyses,
#          'responsible' : who,
#      },
  }
  sample_dict = {
    'data_%s' % pd : subsample_dict.keys()
  }
  return subsample_dict, sample_dict

# Build all the PDs we use
data_DoubleMu, list_DoubleMu = build_data_set('DoubleMu', ['VH', 'Mu'], 'Ian')
datadefs.update(data_DoubleMu)
data_name_map.update(list_DoubleMu)

data_MuEG, list_MuEG = build_data_set('MuEG', ['VH', 'HTT', 'Mu'], 'Evan')
datadefs.update(data_MuEG)
data_name_map.update(list_MuEG)

data_DoubleE, list_DoubleE = build_data_set('DoubleElectron', ['VH',], 'Ian')
datadefs.update(data_DoubleE)
data_name_map.update(list_DoubleE)

data_SingleMu, list_SingleMu = build_data_set('SingleMu', ['Tau', 'Mu'], 'Maria')
datadefs.update(data_SingleMu)
data_name_map.update(list_SingleMu)

data_SingleElectron, list_SingleElectron = build_data_set('SingleElectron', ['Tau', 'E', 'Wjets'], 'Maria')
datadefs.update(data_SingleElectron)
data_name_map.update(list_SingleElectron)

data_TauPlusX, list_TauPlusX = build_data_set('TauPlusX', ['HTT', ], 'Josh')
datadefs.update(data_TauPlusX)
data_name_map.update(list_TauPlusX)

    # Add all the embedded samples
embedded_samples = [
    "/DoubleMu/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau116_ptmu1_13had1_17_v2-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011B_PR_v1_embedded_trans1_tau115_ptelec1_17had1_17_v2-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_PR_v4_embedded_trans1_tau116_ptmu1_13had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_PR_v4_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_May10thRR_v1_embedded_trans1_tau116_ptmu1_13had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_May10thRR_v1_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_Aug05thRR_v1_embedded_trans1_tau116_ptmu1_13had1_17_v2-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_Aug05thRR_v1_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_03Oct2011_v1_embedded_trans1_tau116_ptmu1_13had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/USER",
    "/DoubleMu/StoreResults-DoubleMu_2011A_03Oct2011_v1_embedded_trans1_tau115_ptelec1_17had1_17_v1-f456bdbb960236e5c696adfe9b04eaae/USER",
]
_embed_matcher = re.compile('.*(?P<name>2011.*)_v\d-f4.*')
for embedded_sample in embedded_samples:
    match = _embed_matcher.match(embedded_sample)
    assert(match)
    name = match.group('name')
    datadefs[name] = {
        'datasetpath' : embedded_sample,
        'analyses' : ['HTT'],
        'responsible' : 'Evan',
        'xsec' : -999,
        'pu' : 'data',
    }
