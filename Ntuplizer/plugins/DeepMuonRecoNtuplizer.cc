// -*- C++ -*-
// Package:    DeepMuonReco/Ntuplizer
// Class:      DeepMuonRecoNtuplizer

#include <memory>
#include <vector>
#include <map>
#include <iostream>
#include <cmath>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "TTree.h"

class DeepMuonRecoNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DeepMuonRecoNtuplizer(const edm::ParameterSet&);
  ~DeepMuonRecoNtuplizer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void clearVectors();

  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<TrackingParticleCollection> tpToken_;
  edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> associatorToken_;
  
  edm::EDGetTokenT<RPCRecHitCollection> rpcHitToken_;
  edm::EDGetTokenT<GEMRecHitCollection> gemHitToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> dtSegmentToken_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentToken_;

  edm::ESGetToken<RPCGeometry, MuonGeometryRecord> rpcGeomToken_;
  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomToken_;
  edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;

  TTree* tree_;

  std::vector<double> trackPt, trackEta, trackPhi; 
  std::vector<double> trackPx, trackPy, trackPz;
  std::vector<double> trackVx, trackVy, trackVz;
  std::vector<double> trackQOverP, trackLambda;
  std::vector<double> trackDxy, trackDsz;
  
  std::vector<double> trackQOverPErr, trackLambdaErr;
  std::vector<double> trackPhiErr, trackDxyErr, trackDszErr;
  
  std::vector<double> trackChi2, trackNdof;
  std::vector<int> trackCharge;
  std::vector<int> trackNAlgo;

  std::vector<int> trackIsGoodTrack;
  std::vector<int> trackIsRecoMuon; 
  std::vector<int> trackIsTrkMuon;
  std::vector<int> trackIsGlbMuon;  
  std::vector<int> trackIsPFMuon;

  std::vector<int> trackIsMatchedMuon;
  std::vector<int> trackMatchTpIdx;
  std::vector<float> trackMatchQuality;

  std::vector<double> tpPt, tpEta, tpPhi;
  std::vector<int> tpPdgId, tpCharge, tpStatus;

  std::vector<int> rpcHitRawId;
  std::vector<double> rpcHitPosX, rpcHitPosY, rpcHitPosZ;
  std::vector<double> rpcHitPosXErr, rpcHitPosYErr;
  std::vector<int> rpcHitCls, rpcHitBX;

  std::vector<int> gemHitRawId;
  std::vector<double> gemHitPosX, gemHitPosY, gemHitPosZ;
  std::vector<double> gemHitPosXErr, gemHitPosYErr;
  std::vector<int> gemHitCls, gemHitBX;

  std::vector<int> dtSegRawId;
  std::vector<double> dtSegPosX, dtSegPosY, dtSegPosZ;
  std::vector<double> dtSegPosXErr, dtSegPosYErr;
  std::vector<double> dtSegDirX, dtSegDirY, dtSegDirZ;
  std::vector<double> dtSegDirXErr, dtSegDirYErr;
  std::vector<double> dtSegChi2, dtSegNdof;

  std::vector<int> cscSegRawId;
  std::vector<double> cscSegPosX, cscSegPosY, cscSegPosZ;
  std::vector<double> cscSegPosXErr, cscSegPosYErr;
  std::vector<double> cscSegDirX, cscSegDirY, cscSegDirZ;
  std::vector<double> cscSegDirXErr, cscSegDirYErr;
  std::vector<double> cscSegChi2, cscSegNdof, cscSegTime;
};

DeepMuonRecoNtuplizer::DeepMuonRecoNtuplizer(const edm::ParameterSet& iConfig)
  : muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
    tpToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
    associatorToken_(consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getParameter<edm::InputTag>("associator"))),
    
    rpcHitToken_(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("rpcRecHits"))),
    gemHitToken_(consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"))),
    dtSegmentToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("dtSegments"))),
    cscSegmentToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegments"))),
    
    rpcGeomToken_(esConsumes()),
    gemGeomToken_(esConsumes()),
    dtGeomToken_(esConsumes()),
    cscGeomToken_(esConsumes()) {
  usesResource("TFileService");
}

void DeepMuonRecoNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "DeepMuonReco Tree");

  tree_->Branch("track_pt", &trackPt);   
  tree_->Branch("track_eta", &trackEta);  
  tree_->Branch("track_phi", &trackPhi);
  tree_->Branch("track_px", &trackPx);
  tree_->Branch("track_py", &trackPy);
  tree_->Branch("track_pz", &trackPz);
  tree_->Branch("track_vx", &trackVx);
  tree_->Branch("track_vy", &trackVy);
  tree_->Branch("track_vz", &trackVz);
  tree_->Branch("track_qoverp", &trackQOverP);
  tree_->Branch("track_lamda", &trackLambda);
  tree_->Branch("track_dxy", &trackDxy);
  tree_->Branch("track_dsz", &trackDsz);
  tree_->Branch("track_qoverp_err", &trackQOverPErr);
  tree_->Branch("track_lambda_err", &trackLambdaErr);
  tree_->Branch("track_phi_err", &trackPhiErr);
  tree_->Branch("track_dxy_err", &trackDxyErr);
  tree_->Branch("track_dsz_err", &trackDszErr);
  tree_->Branch("track_charge", &trackCharge);
  tree_->Branch("track_chi2", &trackChi2);
  tree_->Branch("track_ndof", &trackNdof);
  tree_->Branch("track_n_algo", &trackNAlgo);
  tree_->Branch("track_is_good_track", &trackIsGoodTrack);
  tree_->Branch("track_is_reco_muon", &trackIsRecoMuon);
  tree_->Branch("track_is_trk_muon", &trackIsTrkMuon);
  tree_->Branch("track_is_glb_muon", &trackIsGlbMuon);
  tree_->Branch("track_is_pf_muon", &trackIsPFMuon);
  tree_->Branch("track_is_matched_muon", &trackIsMatchedMuon);
  tree_->Branch("track_match_quality", &trackMatchQuality);
  tree_->Branch("track_match_tp_idx", &trackMatchTpIdx);

  tree_->Branch("tp_pt", &tpPt);
  tree_->Branch("tp_eta", &tpEta);
  tree_->Branch("tp_phi", &tpPhi);
  tree_->Branch("tp_pdg_id", &tpPdgId);
  tree_->Branch("tp_charge", &tpCharge);
  tree_->Branch("tp_status", &tpStatus);

  tree_->Branch("rpc_hit_rawid", &rpcHitRawId);
  tree_->Branch("rpc_hit_pos_x", &rpcHitPosX);
  tree_->Branch("rpc_hit_pos_y", &rpcHitPosY);
  tree_->Branch("rpc_hit_pos_z", &rpcHitPosZ);
  tree_->Branch("rpc_hit_pos_x_err", &rpcHitPosXErr);
  tree_->Branch("rpc_hit_pos_y_err", &rpcHitPosYErr);
  tree_->Branch("rpc_hit_cls", &rpcHitCls);
  tree_->Branch("rpc_hit_bx", &rpcHitBX);

  tree_->Branch("gem_hit_rawid", &gemHitRawId);
  tree_->Branch("gem_hit_pos_x", &gemHitPosX);
  tree_->Branch("gem_hit_pos_y", &gemHitPosY);
  tree_->Branch("gem_hit_pos_z", &gemHitPosZ);
  tree_->Branch("gem_hit_pos_x_err", &gemHitPosXErr);
  tree_->Branch("gem_hit_pos_y_err", &gemHitPosYErr);
  tree_->Branch("gem_hit_cls", &gemHitCls);
  tree_->Branch("gem_hit_bx", &gemHitBX);

  tree_->Branch("dt_seg_rawid", &dtSegRawId);
  tree_->Branch("dt_seg_pos_x", &dtSegPosX);
  tree_->Branch("dt_seg_pos_y", &dtSegPosY);
  tree_->Branch("dt_seg_pos_z", &dtSegPosZ);
  tree_->Branch("dt_seg_pos_x_err", &dtSegPosXErr);
  tree_->Branch("dt_seg_pos_y_err", &dtSegPosYErr);
  tree_->Branch("dt_seg_dir_x", &dtSegDirX);
  tree_->Branch("dt_seg_dir_y", &dtSegDirY);
  tree_->Branch("dt_seg_dir_z", &dtSegDirZ);
  tree_->Branch("dt_seg_dir_x_err", &dtSegDirXErr);
  tree_->Branch("dt_seg_dir_y_err", &dtSegDirYErr);
  tree_->Branch("dt_seg_chi2", &dtSegChi2);
  tree_->Branch("dt_seg_ndof", &dtSegNdof);

  tree_->Branch("csc_seg_rawid", &cscSegRawId);
  tree_->Branch("csc_seg_pos_x", &cscSegPosX);
  tree_->Branch("csc_seg_pos_y", &cscSegPosY);
  tree_->Branch("csc_seg_pos_z", &cscSegPosZ);
  tree_->Branch("csc_seg_pos_x_err", &cscSegPosXErr);
  tree_->Branch("csc_seg_pos_y_err", &cscSegPosYErr);
  tree_->Branch("csc_seg_dir_x", &cscSegDirX);
  tree_->Branch("csc_seg_dir_y", &cscSegDirY);
  tree_->Branch("csc_seg_dir_z", &cscSegDirZ);
  tree_->Branch("csc_seg_dir_x_err", &cscSegDirXErr);
  tree_->Branch("csc_seg_dir_y_err", &cscSegDirYErr);
  tree_->Branch("csc_seg_chi2", &cscSegChi2);
  tree_->Branch("csc_seg_ndof", &cscSegNdof);
  tree_->Branch("csc_seg_time", &cscSegTime);
}

void DeepMuonRecoNtuplizer::clearVectors() {
  trackPt.clear(); trackEta.clear(); trackPhi.clear();
  trackPx.clear(); trackPy.clear(); trackPz.clear();
  trackVx.clear(); trackVy.clear(); trackVz.clear();
  trackQOverP.clear(); trackLambda.clear();
  trackDxy.clear(); trackDsz.clear();
  trackQOverPErr.clear(); trackLambdaErr.clear(); trackPhiErr.clear();
  trackDxyErr.clear(); trackDszErr.clear();
  trackChi2.clear(); trackNdof.clear();
  trackCharge.clear(); trackNAlgo.clear();
  trackIsGoodTrack.clear(); trackIsMatchedMuon.clear();
  trackIsRecoMuon.clear(); trackIsTrkMuon.clear(); 
  trackIsGlbMuon.clear(); trackIsPFMuon.clear();
  trackMatchTpIdx.clear(); trackMatchQuality.clear();

  tpPt.clear(); tpEta.clear(); tpPhi.clear();
  tpPdgId.clear(); tpCharge.clear(); tpStatus.clear();

  rpcHitRawId.clear(); rpcHitPosX.clear(); rpcHitPosY.clear(); rpcHitPosZ.clear();
  rpcHitPosXErr.clear(); rpcHitPosYErr.clear(); rpcHitCls.clear(); rpcHitBX.clear();

  gemHitRawId.clear(); gemHitPosX.clear(); gemHitPosY.clear(); gemHitPosZ.clear();
  gemHitPosXErr.clear(); gemHitPosYErr.clear(); gemHitCls.clear(); gemHitBX.clear();

  dtSegRawId.clear(); dtSegPosX.clear(); dtSegPosY.clear(); dtSegPosZ.clear();
  dtSegPosXErr.clear(); dtSegPosYErr.clear(); 
  dtSegDirX.clear(); dtSegDirY.clear(); dtSegDirZ.clear();
  dtSegDirXErr.clear(); dtSegDirYErr.clear(); dtSegChi2.clear(); dtSegNdof.clear();

  cscSegRawId.clear(); cscSegPosX.clear(); cscSegPosY.clear(); cscSegPosZ.clear();
  cscSegPosXErr.clear(); cscSegPosYErr.clear(); 
  cscSegDirX.clear(); cscSegDirY.clear(); cscSegDirZ.clear();
  cscSegDirXErr.clear(); cscSegDirYErr.clear(); cscSegChi2.clear(); cscSegNdof.clear(); cscSegTime.clear();
}

void DeepMuonRecoNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  clearVectors();
  
  edm::Handle<RPCRecHitCollection> rpcHandle;
  iEvent.getByToken(rpcHitToken_, rpcHandle);

  edm::Handle<GEMRecHitCollection> gemHandle;
  iEvent.getByToken(gemHitToken_, gemHandle);

  edm::Handle<DTRecSegment4DCollection> dtHandle;
  iEvent.getByToken(dtSegmentToken_, dtHandle);

  edm::Handle<CSCSegmentCollection> cscHandle;
  iEvent.getByToken(cscSegmentToken_, cscHandle);

  const auto& rpcGeom = iSetup.getData(rpcGeomToken_);
  const auto& gemGeom = iSetup.getData(gemGeomToken_);
  const auto& dtGeom  = iSetup.getData(dtGeomToken_);
  const auto& cscGeom = iSetup.getData(cscGeomToken_);

  if (rpcHandle.isValid()) {
    for (const auto &hit : *rpcHandle) {
      RPCDetId rpcId(hit.geographicalId());
      const auto rpcDet = rpcGeom.idToDet(rpcId);

      LocalPoint localPos = hit.localPosition();
      LocalError localPosErr = hit.localPositionError();
      GlobalPoint globalPos = rpcDet->surface().toGlobal(localPos);

      rpcHitRawId.push_back(hit.geographicalId());
      rpcHitPosX.push_back(globalPos.x());
      rpcHitPosY.push_back(globalPos.y());
      rpcHitPosZ.push_back(globalPos.z());
      rpcHitPosXErr.push_back(std::sqrt(localPosErr.xx()));
      rpcHitPosYErr.push_back(std::sqrt(localPosErr.yy()));
      rpcHitCls.push_back(hit.clusterSize());
      rpcHitBX.push_back(hit.BunchX());
    }
  }

  if (gemHandle.isValid()) {
    for (const auto &hit : *gemHandle) {
      GEMDetId gemId(hit.geographicalId());
      const auto gemDet = gemGeom.idToDet(gemId);

      LocalPoint localPos = hit.localPosition();
      LocalError localPosErr = hit.localPositionError();
      GlobalPoint globalPos = gemDet->surface().toGlobal(localPos);

      gemHitRawId.push_back(hit.geographicalId());
      gemHitPosX.push_back(globalPos.x());
      gemHitPosY.push_back(globalPos.y());
      gemHitPosZ.push_back(globalPos.z());
      gemHitPosXErr.push_back(std::sqrt(localPosErr.xx()));
      gemHitPosYErr.push_back(std::sqrt(localPosErr.yy()));
      gemHitCls.push_back(hit.clusterSize());
      gemHitBX.push_back(hit.BunchX());
    }
  }

  if (dtHandle.isValid()) {
    for (const auto &seg : *dtHandle) {
      DTChamberId dtId(seg.geographicalId());
      const auto dtDet = dtGeom.idToDet(dtId);

      LocalPoint localPos = seg.localPosition();
      LocalError localPosErr = seg.localPositionError();
      GlobalPoint globalPos = dtDet->surface().toGlobal(localPos);

      LocalVector localDir = seg.localDirection();
      LocalError localDirErr = seg.localDirectionError();
      GlobalVector globalDir = dtDet->surface().toGlobal(localDir);

      dtSegRawId.push_back(seg.geographicalId());
      dtSegPosX.push_back(globalPos.x());
      dtSegPosY.push_back(globalPos.y());
      dtSegPosZ.push_back(globalPos.z());
      dtSegPosXErr.push_back(std::sqrt(localPosErr.xx()));
      dtSegPosYErr.push_back(std::sqrt(localPosErr.yy()));
      dtSegDirX.push_back(globalDir.x());
      dtSegDirY.push_back(globalDir.y());
      dtSegDirZ.push_back(globalDir.z());
      dtSegDirXErr.push_back(std::sqrt(localDirErr.xx()));
      dtSegDirYErr.push_back(std::sqrt(localDirErr.yy()));
      dtSegChi2.push_back(seg.chi2());
      dtSegNdof.push_back(seg.degreesOfFreedom());
    }
  }

  if (cscHandle.isValid()) {
    for (const auto &seg : *cscHandle) {
      CSCDetId cscId(seg.geographicalId());
      const auto cscDet = cscGeom.idToDet(cscId);

      LocalPoint localPos = seg.localPosition();
      LocalError localPosErr = seg.localPositionError();
      GlobalPoint globalPos = cscDet->surface().toGlobal(localPos);

      LocalVector localDir = seg.localDirection();
      LocalError localDirErr = seg.localDirectionError();
      GlobalVector globalDir = cscDet->surface().toGlobal(localDir);

      cscSegRawId.push_back(seg.geographicalId());
      cscSegPosX.push_back(globalPos.x());
      cscSegPosY.push_back(globalPos.y());
      cscSegPosZ.push_back(globalPos.z());
      cscSegPosXErr.push_back(std::sqrt(localPosErr.xx()));
      cscSegPosYErr.push_back(std::sqrt(localPosErr.yy()));
      cscSegDirX.push_back(globalDir.x());
      cscSegDirY.push_back(globalDir.y());
      cscSegDirZ.push_back(globalDir.z());
      cscSegDirXErr.push_back(std::sqrt(localDirErr.xx()));
      cscSegDirYErr.push_back(std::sqrt(localDirErr.yy()));
      cscSegChi2.push_back(seg.chi2());
      cscSegNdof.push_back(seg.degreesOfFreedom());
      cscSegTime.push_back(seg.time());
    }
  }

  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackToken_, tracks);

  edm::Handle<TrackingParticleCollection> tps;
  iEvent.getByToken(tpToken_, tps);

  edm::Handle<reco::TrackToTrackingParticleAssociator> associatorHandle;
  iEvent.getByToken(associatorToken_, associatorHandle);
  const auto& associator = *associatorHandle;

  edm::RefToBaseVector<reco::Track> trackRefs;
  for (size_t i = 0; i < tracks->size(); ++i) {
    trackRefs.push_back(edm::RefToBase<reco::Track>(reco::TrackRef(tracks, i)));
  }
  
  edm::RefVector<TrackingParticleCollection> tpRefs;
  for (size_t i = 0; i < tps->size(); ++i) {
    tpRefs.push_back(TrackingParticleRef(tps, i));
  }

  reco::RecoToSimCollection recoToSims = associator.associateRecoToSim(trackRefs, tpRefs);
  
  std::map<unsigned int, const reco::Muon*> trackToMuonMap;
  for (const auto& mu : *muons) {
    if (mu.innerTrack().isNonnull()) {
      trackToMuonMap[mu.innerTrack().key()] = &mu;
    }
  }

  for (size_t i = 0; i < tracks->size(); ++i) {
    reco::TrackRef trkRef(tracks, i);

    trackPt.push_back(trkRef->pt());
    trackEta.push_back(trkRef->eta());
    trackPhi.push_back(trkRef->phi());
    trackPx.push_back(trkRef->px());
    trackPy.push_back(trkRef->py());
    trackPz.push_back(trkRef->pz());
    
    trackVx.push_back(trkRef->vx());
    trackVy.push_back(trkRef->vy());
    trackVz.push_back(trkRef->vz());
    
    trackQOverP.push_back(trkRef->qoverp());
    trackLambda.push_back(trkRef->lambda());
    trackDxy.push_back(trkRef->dxy());
    trackDsz.push_back(trkRef->dsz());
    
    trackQOverPErr.push_back(trkRef->qoverpError());
    trackLambdaErr.push_back(trkRef->lambdaError());
    trackPhiErr.push_back(trkRef->phiError());
    trackDxyErr.push_back(trkRef->dxyError());
    trackDszErr.push_back(trkRef->dszError());

    trackChi2.push_back(trkRef->chi2());
    trackNdof.push_back(trkRef->ndof());
    trackCharge.push_back(trkRef->charge());
    trackNAlgo.push_back(trkRef->algo()); 

    int isGoodTrack = 0;
    int isRecoMuon = 0;
    int isTrkMuon = 0;
    int isGlbMuon = 0;
    int isPFMuon = 0;

    if ((trkRef->pt() > 0.5) && (trkRef->p() > 2.5) && (std::abs(trkRef->eta()) < 3.0)) {
      isGoodTrack = 1;
    }

    if (trackToMuonMap.count(trkRef.key())) {
      const reco::Muon* mu = trackToMuonMap[trkRef.key()];
      isRecoMuon = 1;
      isTrkMuon = mu->isTrackerMuon();
      isGlbMuon = mu->isGlobalMuon();
      isPFMuon = mu->isPFMuon();
    }

    trackIsRecoMuon.push_back(isRecoMuon);
    trackIsTrkMuon.push_back(isTrkMuon);
    trackIsGlbMuon.push_back(isGlbMuon);
    trackIsPFMuon.push_back(isPFMuon);

    int isMatchedMuon = 0;
    int matchedTpIdx = -1;
    float matchQual = 0.0;
    edm::RefToBase<reco::Track> trkRefBase(trkRef);
    if (recoToSims.find(trkRefBase) != recoToSims.end()) {
      const auto& tpQualPairs = recoToSims[trkRefBase];
      if (!tpQualPairs.empty()) {
        const auto& tpQualPair = tpQualPairs.front();
        const auto& bestTpRef = tpQualPair.first;
        matchQual = tpQualPair.second;
        
        if (std::abs(bestTpRef->pdgId()) == 13 && bestTpRef->status() == 1) {
          isMatchedMuon = 1;
        }
        
        for(size_t k = 0; k < tps->size(); ++k) {
          if(bestTpRef == TrackingParticleRef(tps, k)) {
            matchedTpIdx = k;
            break;
          }
        }
      }
    }
    trackIsGoodTrack.push_back(isGoodTrack);
    trackIsMatchedMuon.push_back(isMatchedMuon);
    trackMatchTpIdx.push_back(matchedTpIdx);
    trackMatchQuality.push_back(matchQual);
  }

  for (size_t i = 0; i < tps->size(); ++i) {
    TrackingParticleRef tpRef(tps, i);
    tpPt.push_back(tpRef->pt());
    tpEta.push_back(tpRef->eta());
    tpPhi.push_back(tpRef->phi());
    tpPdgId.push_back(tpRef->pdgId());
    tpCharge.push_back(tpRef->charge());
    tpStatus.push_back(tpRef->status());
  }

  tree_->Fill();
}

void DeepMuonRecoNtuplizer::endJob() {}

void DeepMuonRecoNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("muons", edm::InputTag("muons"));
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("trackingParticles", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("associator", edm::InputTag("quickTrackAssociatorByHits"));
  
  desc.add<edm::InputTag>("rpcRecHits", edm::InputTag("rpcRecHits"));
  desc.add<edm::InputTag>("gemRecHits", edm::InputTag("gemRecHits"));
  desc.add<edm::InputTag>("dtSegments", edm::InputTag("dt4DSegments"));
  desc.add<edm::InputTag>("cscSegments", edm::InputTag("cscSegments"));
  
  descriptions.add("deepMuonRecoNtuplizer", desc);
}

DEFINE_FWK_MODULE(DeepMuonRecoNtuplizer);