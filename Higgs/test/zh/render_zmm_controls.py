#!/usr/bin/env python
'''

Render the Zmumu control plots

Author: Evan K. Friis, UW

'''

from RecoLuminosity.LumiDB import argparse
import json
import logging
import os
import re
import sys

# Steal the args so ROOT doesn't mess them up!
args = sys.argv[:]
sys.argv = [sys.argv[0]]

import rootpy.io as io
from rootpy.plotting import views
from rootpy.plotting import Canvas

import FinalStateAnalysis.MetaData.data_views import data_views

if __name__ == "__main__":
    log = logging.getLogger("render_fake_rate_controls")
    view_builder = logging.getLogger("get_views")
    view_builder.setLevel(logging.DEBUG)
    logging.basicConfig(level=logging.DEBUG)

    parser = argparse.ArgumentParser()
    parser.add_argument('meta', help='File with meta information')
    parser.add_argument('output', help='Output directory')
    parser.add_argument('files', metavar='file', nargs='+',
                        help = 'Histogram files')
    args = parser.parse_args(args[1:])

    meta_info = None
    log.info("Opening meta file: %s", args.meta)
    with open(args.meta) as meta_file:
        meta_info = json.load(meta_file)

    log.info("Building views")
    data_views = data_views.data_views(
        args.files,
        'data_DoubleMu',
    )
    data_view = views.FunctorView(data_views['data_DoubleMu']['view'],
                                  lambda x: x.Rebin(1))

    mc_samples_to_stack = ['Zjets', 'QCDMu', 'Wjets', 'ttjets']
    # Build the MC stack
    mc_view = views.StackView(*[
        views.FunctorView(data_views[x]['view'], lambda x: x.Rebin(1))
        for x in mc_samples_to_stack], sort=True)

    log.debug("Zjets unweighted entries: %f",
              data_views['Zjets']['unweighted_view'].Get(
                  'plots/entry_count').Integral()
             )


    # We just need to figure out the directory structure from any old file
    layout_filename = data_views.values()[0]['subsamples'].values()[0]['filename']
    log.info("Getting file layout from %s", layout_filename)

    canvas = Canvas()

    plot_list = open(os.path.join(args.output, 'plot_list.txt'), 'w')

    with io.open(layout_filename, 'r') as layout_file:
        for path, subdirs, histos in layout_file.walk(class_pattern='TH1F'):
            for histo in histos:
                data_histo = data_view.Get(os.path.join(path, histo))
                mc_histo = mc_view.Get(os.path.join(path, histo))
                mc_histo.Draw()
                data_histo.Draw('same')
                mc_histo.SetMaximum(1.2*max(
                    mc_histo.GetMaximum(), data_histo.GetMaximum()))

                plot_filename = os.path.join(
                    args.output,
                    path.replace('/', '_') + '_' + histo + '.pdf'
                )
                plot_list.write(plot_filename + '\n')
                canvas.SaveAs(plot_filename)
