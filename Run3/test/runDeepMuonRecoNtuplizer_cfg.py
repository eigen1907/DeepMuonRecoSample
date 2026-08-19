import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.AlCa.GlobalTag import GlobalTag


options = VarParsing("analysis")

options.setType("outputFile", VarParsing.varType.string)
options.setDefault("inputFiles", ["file:reco.root"])
options.setDefault("outputFile", "ntuple.root")

options.register(
    "isMC",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on MC if True, collision data if False",
)

options.parseArguments()


process = cms.Process("DeepMuonRecoSample", Run3_2024)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


if options.isMC:
    global_tag = "140X_mcRun3_2024_realistic_v26"
else:
    global_tag = "140X_dataRun3_v17"

process.GlobalTag = GlobalTag(
    process.GlobalTag,
    global_tag,
    "",
)


process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
)


process.load(
    "DeepMuonRecoSample.Ntuplizer.deepMuonRecoNtuplizer_cfi"
)

process.deepMuonRecoNtuplizer.isMC = cms.bool(options.isMC)


if options.isMC:
    process.load(
        "SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi"
    )
    process.load(
        "SimTracker.TrackerHitAssociation.tpClusterProducer_cfi"
    )

    process.tpClusterProducer.trackingParticleSrc = cms.InputTag(
        "mix",
        "MergedTrackTruth",
    )
    process.tpClusterProducer.pixelSimLinkSrc = cms.InputTag(
        "simSiPixelDigis"
    )
    process.tpClusterProducer.stripSimLinkSrc = cms.InputTag(
        "simSiStripDigis"
    )
    process.tpClusterProducer.pixelClusterSrc = cms.InputTag(
        "siPixelClusters"
    )
    process.tpClusterProducer.stripClusterSrc = cms.InputTag(
        "siStripClusters"
    )

    process.deepMuonRecoNtuplizer.trackingParticles = cms.InputTag(
        "mix",
        "MergedTrackTruth",
    )


process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string(options.outputFile),
)


if options.isMC:
    process.p = cms.Path(
        process.tpClusterProducer
        * process.quickTrackAssociatorByHits
        * process.deepMuonRecoNtuplizer
    )
else:
    process.p = cms.Path(
        process.deepMuonRecoNtuplizer
    )