# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step3_RECO --processName RERECO --eventcontent RECOSIM --datatier GEN-SIM-RECO --conditions 140X_mcRun3_2024_realistic_v26 --step RAW2DIGI,L1Reco,RECO,RECOSIM --geometry DB:Extended --era Run3_2024 --filein file:digiraw.root --fileout file:reco.root --mc --number 10 --no_exec --outputCommand drop *,keep *_muons_*_*,keep *_generalTracks_*_*,keep *_offlineBeamSpot_*_*,keep *_mix_MergedTrackTruth_*,keep *_genParticles_*_*,keep *_rpcRecHits_*_*,keep *_gemRecHits_*_*,keep *_gemSegments_*_*,keep *_dt4DSegments_*_*,keep *_cscSegments_*_*,keep *_simSiPixelDigis_*_*,keep *_simSiStripDigis_*_*,keep *_siPixelClusters_*_*,keep *_siStripClusters_*_*
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.setType('outputFile', VarParsing.varType.string)
options.setDefault('inputFiles', ['file:digiraw.root'])
options.setDefault('outputFile', 'reco.root')
options.parseArguments()

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024

process = cms.Process('RERECO',Run3_2024)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3_RECO nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun3_2024_realistic_v26', '')
process.RECOSIMoutput.outputCommands.append('drop *')
process.RECOSIMoutput.outputCommands.append('keep *_muons_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_generalTracks_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_offlineBeamSpot_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_mix_MergedTrackTruth_*')
process.RECOSIMoutput.outputCommands.append('keep *_genParticles_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_rpcRecHits_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_gemRecHits_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_gemSegments_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_dt4DSegments_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_cscSegments_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_simSiPixelDigis_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_simSiStripDigis_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_siPixelClusters_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_siStripClusters_*_*')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.endjob_step,process.RECOSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)



# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
