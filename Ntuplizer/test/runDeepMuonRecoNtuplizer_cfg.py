import FWCore.ParameterSet.Config as cms

process = cms.Process("DeepMuonRecoSample")

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun4_realistic_v4', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

process.tpClusterProducer.pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel")
process.tpClusterProducer.phase2OTSimLinkSrc = cms.InputTag("simSiPixelDigis", "Tracker")
process.tpClusterProducer.pixelClusterSrc = cms.InputTag("siPixelClusters")
process.tpClusterProducer.phase2OTClusterSrc = cms.InputTag("siPhase2Clusters")

process.load("DeepMuonRecoSample.Ntuplizer.deepMuonRecoNtuplizer_cfi")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.p = cms.Path(
    process.tpClusterProducer *
    process.quickTrackAssociatorByHits *
    process.deepMuonRecoNtuplizer
)
