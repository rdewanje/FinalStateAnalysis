'''

Make inclusive e-mu (Z + ttbar) control plots

'''

import os
import glob
from FinalStateAnalysis.PlotTools.Plotter import Plotter

jobid = os.environ['jobid']

output_dir = os.path.join('results', jobid, 'plots', 'em')

samples = [
    'Zjets_M50',
    'WZ*',
    'WW*',
    'ZZ*',
    'TT*',
    'WplusJets*',
    "data_MuEG*",
]

files = []
lumifiles = []

for x in samples:
    files.extend(glob.glob('results/%s/AnalyzeEM/%s.root' % (jobid, x)))
    lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))

plotter = Plotter(files, lumifiles, output_dir)

# Override ordering
plotter.mc_samples = [
    'TTplusJets_madgraph',
    'WplusJets_madgraph',
    'Zjets_M50',
    'WZJetsTo3LNu*',
    'WW*',
    'ZZJetsTo4L*',
]

sqrts = 7 if '7TeV' in jobid else 8

plotter.plot_mc_vs_data('em', 'emMass', rebin=10, leftside=False,
                        xaxis='m_{e#mu} (GeV)')
plotter.add_cms_blurb(sqrts)
plotter.save('mass')

plotter.plot_mc_vs_data('em', 'mPt')
plotter.save('mPt')
plotter.plot_mc_vs_data('em', 'ePt')
plotter.save('ePt')

plotter.plot_mc_vs_data('em', 'mAbsEta')
plotter.save('mAbsEta')
plotter.plot_mc_vs_data('em', 'eAbsEta')
plotter.save('eAbsEta')

plotter.plot_mc_vs_data('em', 'nvtx')
plotter.save('nvtx')

plotter.plot_mc_vs_data('em', 'bjetCSVVeto')
plotter.save('bjetCSVVeto')

plotter.plot_mc_vs_data('em', 'bjetVeto')
plotter.save('bjetVeto')
