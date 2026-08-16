import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import (
    pythia8CommonSettingsBlock,
)
from Configuration.Generator.MCTunesRun3ECM13p6TeV.PythiaCP5Settings_cfi import (
    pythia8CP5SettingsBlock,
)


# Unfiltered Run 3 minimum-bias fragment. CMSSW_14_0_21_patch1 contains the
# Run 3 CP5 tune but does not ship a named 13.6 TeV MinBias fragment.
generator = cms.EDFilter(
    "Pythia8ConcurrentGeneratorFilter",
    filterEfficiency=cms.untracked.double(1.0),
    maxEventsToPrint=cms.untracked.int32(0),
    pythiaHepMCVerbosity=cms.untracked.bool(False),
    pythiaPylistVerbosity=cms.untracked.int32(0),
    comEnergy=cms.double(13600.0),
    PythiaParameters=cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters=cms.vstring("SoftQCD:inelastic = on"),
        parameterSets=cms.vstring(
            "pythia8CommonSettings",
            "pythia8CP5Settings",
            "processParameters",
        ),
    ),
)

ProductionFilterSequence = cms.Sequence(generator)
