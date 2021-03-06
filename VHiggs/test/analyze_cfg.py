import FWCore.ParameterSet.Config as cms

import FinalStateAnalysis.Utilities.TauVarParsing as TauVarParsing
options = TauVarParsing.TauVarParsing(
    skipEvents=0, # For debugging
    puScenario='S4',
)

options.outputFile="vhiggs.root"
options.parseArguments()

process = cms.Process("VHiggsAna2")

# Input in FWLITE mode
process.fwliteInput = cms.PSet(fileNames = cms.vstring(options.inputFiles))
# Input in FWK mode
process.source = cms.Source(
    "PoolSource", fileNames = cms.untracked.vstring(options.inputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents),
    #lumisToProcess = options.buildPoolSourceLumiMask()
)

if options.eventsToProcess:
    process.source.eventsToProcess = cms.untracked.VEventRange(
        options.eventsToProcess)

plots_filename = options.outputFile.replace('.root', '.plots.root')
process.fwliteOutput = cms.PSet(fileName = cms.string(options.outputFile))
process.TFileService = cms.Service(
    "TFileService", fileName = cms.string(plots_filename)
)

process.steering = cms.PSet(
    analyzers = cms.vstring(
        #'eee',
        #'eem',
        #'eet',
        'emm',
        'emt',
        #'ett',
        'mmm',
        'mmt',
        #'mtt'
    ),
    reportAfter = cms.uint32(1000),
    ignored_cuts = cms.vstring()
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents))

pu_weight = 'puWeight("2011AB", "%s")' % options.puScenario
if not options.isMC:
    pu_weight = '1.0'

# Common among all analyzers
process.common = cms.PSet(
    weights = cms.vstring(
        pu_weight
    ),
    evtSrc = cms.InputTag("patFinalStateEventProducer"),
    skimCounter = cms.InputTag("eventCount", "", "TUPLE"),
)

process.steering.ignored_cuts = cms.vstring(
    "*ID*", "*_BJetVeto", "*RelIso*", "*Veto*", "*CombinedIsolation*"
)

from FinalStateAnalysis.VHiggs.selectionEEE import selections as selectionsEEE
from FinalStateAnalysis.VHiggs.selectionEEM import selections as selectionsEEM
from FinalStateAnalysis.VHiggs.selectionEET import selections as selectionsEET
from FinalStateAnalysis.VHiggs.selectionEMM import selections as selectionsEMM
from FinalStateAnalysis.VHiggs.selectionEMT import selections as selectionsEMT
from FinalStateAnalysis.VHiggs.selectionETT import selections as selectionsETT
from FinalStateAnalysis.VHiggs.selectionMMM import selections as selectionsMMM
from FinalStateAnalysis.VHiggs.selectionMMT import selections as selectionsMMT
from FinalStateAnalysis.VHiggs.selectionMTT import selections as selectionsMTT

from FinalStateAnalysis.VHiggs.selectionEEE import plots as plotsEEE
from FinalStateAnalysis.VHiggs.selectionEEM import plots as plotsEEM
from FinalStateAnalysis.VHiggs.selectionEET import plots as plotsEET
from FinalStateAnalysis.VHiggs.selectionEMM import plots as plotsEMM
from FinalStateAnalysis.VHiggs.selectionEMT import plots as plotsEMT
from FinalStateAnalysis.VHiggs.selectionETT import plots as plotsETT
from FinalStateAnalysis.VHiggs.selectionMMM import plots as plotsMMM
from FinalStateAnalysis.VHiggs.selectionMMT import plots as plotsMMT
from FinalStateAnalysis.VHiggs.selectionMTT import plots as plotsMTT

process.eee = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateElecElecElec'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            #plot = plotsEEE,
            plot = cms.PSet(histos=cms.VPSet()),
        ),
        selections = selectionsEEE
    )
)

process.eem = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateElecElecMu'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(1).pt'),
            take = cms.uint32(5),
            plot = plotsEEM,
        ),
        selections = selectionsEEM
    )
)

process.eet = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateElecElecTau'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            plot = plotsEET,
        ),
        selections = selectionsEET
    )
)

process.emm = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateElecMuMu'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            plot = plotsEMM,
        ),
        selections = selectionsEMM,
    )
)

process.emt = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateElecMuTau'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            plot = plotsEMT,
        ),
        selections = selectionsEMT,
    )
)

process.ett = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateElecTauTau'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            plot = plotsETT,
        ),
        selections = selectionsETT,
    )
)

process.mmm = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateMuMuMu'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            plot = plotsMMM,
        ),
        selections = selectionsMMM,
    )
)

process.mmt = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateMuMuTau'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            plot = plotsMMT,
        ),
        selections = selectionsMMT,
    )
)

process.mtt = cms.PSet(
    process.common,
    src = cms.InputTag('finalStateMuTauTau'),
    analysis = cms.PSet(
        ignore = process.steering.ignored_cuts,
        final = cms.PSet(
            sort = cms.string('daughter(2).pt'),
            take = cms.uint32(5),
            plot = plotsMTT,
        ),
        selections = selectionsMTT,
    )
)

################################################################################
####   Configure PU weighting   ################################################
################################################################################

# Build the filter selectors for skimming events.
for channel in process.steering.analyzers:
    module = getattr(process, channel)
    print "Removing unused histograms"
    del module.analysis.final.plot.OS
    del module.analysis.final.plot.SS
    if options.isMC:
        scenario = options.puScenario
        print "Configuring %s ntuple for %s PU re-weighting" % (
            channel, scenario)
        module.analysis.final.plot.ntuple.pu2011A = cms.string(
            "evt.puWeight('2011A', '%s')" % scenario
        )
        module.analysis.final.plot.ntuple.pu2011B = cms.string(
            "evt.puWeight('2011B', '%s')" % scenario
        )
        module.analysis.final.plot.ntuple.pu2011AB = cms.string(
            "evt.puWeight('2011AB', '%s')" % scenario
        )
    else:
        module.analysis.final.plot.ntuple.pu2011A = cms.string("1.0")
        module.analysis.final.plot.ntuple.pu2011B = cms.string("1.0")
        module.analysis.final.plot.ntuple.pu2011AB = cms.string("1.0")


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("FinalStateAnalysis.RecoTools.eventCount_cfi")
process.count = cms.Path(process.eventCount)
process.schedule = cms.Schedule(process.count)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring()
    ),
    outputCommands = cms.untracked.vstring('keep *')
)
process.end = cms.EndPath(process.out)

process.ProfilerService = cms.Service (
    "ProfilerService",
    firstEvent = cms.untracked.int32(3),
    lastEvent = cms.untracked.int32(12),
    paths = cms.untracked.vstring('emtPath')
)

# Build the filter selectors for skimming events.
for channel in process.steering.analyzers:
    # Build the filter
    filter = cms.EDFilter(
        "PATFinalStateAnalysisFilter",
        getattr(process, channel)
    )
    setattr(process, channel + "Filter", filter)
    path = cms.Path(filter)
    setattr(process, channel + "Path", path)
    process.schedule.append(path)
    process.out.SelectEvents.SelectEvents.append(channel + "Path")

process.schedule.append(process.end)
