import FWCore.ParameterSet.Config as cms


deepMuonRecoNtuplizer = cms.EDAnalyzer(
    "DeepMuonRecoNtuplizer",
    isMC=cms.bool(True),
    muons=cms.InputTag("muons"),
    tracks=cms.InputTag("generalTracks"),
    trackingParticles=cms.InputTag("mix", "MergedTrackTruth"),
    associator=cms.InputTag("quickTrackAssociatorByHits"),
    rpcRecHits=cms.InputTag("rpcRecHits"),
    gemRecHits=cms.InputTag("gemRecHits"),
    gemSegments=cms.InputTag("gemSegments"),
    dtSegments=cms.InputTag("dt4DSegments"),
    cscSegments=cms.InputTag("cscSegments"),
)
