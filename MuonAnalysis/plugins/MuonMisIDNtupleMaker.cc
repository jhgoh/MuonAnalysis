#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

struct SVFitResult
{
  SVFitResult():isValid(false) {}

  bool isValid;
  reco::Candidate::LorentzVector p4, leg1, leg2;
  int q1, q2;
  reco::Particle::Point vertex;
  reco::Vertex::CovarianceMatrix cov;
  double chi2, ndof, lxy, vz;
};

typedef math::XYZTLorentzVector LV;

class MuonMisIDNtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  MuonMisIDNtupleMaker(const edm::ParameterSet& pset);
  virtual ~MuonMisIDNtupleMaker() {};
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  SVFitResult fitSV(const reco::Particle::Point& pvPos, const reco::Vertex::CovarianceMatrix& pvCov,
                    const reco::TransientTrack& transTrack1,
                    const reco::TransientTrack& transTrack2) const;
  int muonIdBit(const reco::Muon& mu, const reco::Vertex& vertex) const;
  int genCategory(const reco::GenParticle& p) const;
  template <typename T>
  std::pair<int, int> matchTwo(const LV& lv1, const LV& lv2, T& coll) const;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

  // Constants for the vertex fit
  const double pionMass = 0.1396;
  const double kaonMass = 0.4937;
  const double protonMass = 0.9383;

  const bool applyGenFilter_, useBeamSpot_;
  int pdgId_;

  double cut_minVtxRawMass_, cut_maxVtxRawMass_;
  double cut_minVtxMass_, cut_maxVtxMass_, cut_minVtxLxy_, cut_maxVtxLxy_;
  const double cut_minTrkPt_, cut_maxTrkEta_;
  const int cut_minTrkNHit_;
  const double cut_maxTrkChi2_, cut_minTrkSigXY_, cut_minTrkSigZ_;
  const double cut_maxVtxDCA_, cut_maxVtxChi2_, cut_minVtxSignif_;

  double mass1_, mass2_;
  int pdgId1_, pdgId2_;

  // Trees and histograms
  typedef std::vector<int> vint;
  typedef std::vector<float> vfloat;

  TTree* tree_;
  unsigned char b_run, b_lumi;
  unsigned long long int b_event;
  //double b_genWeight, b_puWeight;
  unsigned char b_nPV, b_nSV, b_nGen;

  float b_vtx_mass, b_vtx_pt, b_vtx_lxy, b_vtx_vz;
  vfloat* b_trk_pt, * b_trk_eta, * b_trk_phi;
  vint* b_trk_pdgId;

  vfloat* b_mu_pt, * b_mu_dR;
  vint* b_mu_q, * b_mu_id;

  float b_gen_dR;

  TH1D* hN_;
  TH1D* hM_;
};

MuonMisIDNtupleMaker::MuonMisIDNtupleMaker(const edm::ParameterSet& pset):
  applyGenFilter_(pset.getUntrackedParameter<bool>("applyGenFilter", false)),
  useBeamSpot_(pset.getUntrackedParameter<bool>("useBeamSpot", false)),
  cut_minVtxLxy_(pset.getUntrackedParameter<double>("minVtxLxy")),
  cut_maxVtxLxy_(pset.getUntrackedParameter<double>("maxVtxLxy")),
  cut_minTrkPt_(pset.getUntrackedParameter<double>("minTrkPt")),
  cut_maxTrkEta_(pset.getUntrackedParameter<double>("maxTrkEta")),
  cut_minTrkNHit_(pset.getUntrackedParameter<int>("minTrkNHit")),
  cut_maxTrkChi2_(pset.getUntrackedParameter<double>("maxTrkChi2")),
  cut_minTrkSigXY_(pset.getUntrackedParameter<double>("minTrkSigXY")),
  cut_minTrkSigZ_(pset.getUntrackedParameter<double>("minTrkSigZ")),
  cut_maxVtxDCA_(pset.getUntrackedParameter<double>("maxVtxDCA")),
  cut_maxVtxChi2_(pset.getUntrackedParameter<double>("maxVtxChi2")),
  cut_minVtxSignif_(pset.getUntrackedParameter<double>("minVtxSignif"))
{
  genParticleToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticles"));
  vertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"));
  beamSpotToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  trackToken_ = consumes<reco::TrackCollection>(pset.getParameter<edm::InputTag>("tracks"));
  pfCandToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("pfCandidates"));
  muonToken_ = consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));

  const string vtxType = pset.getUntrackedParameter<string>("vtxType");
  if ( vtxType == "kshort" ) {
    pdgId_ = 310;
    pdgId1_ = pdgId2_ = 211;
    mass1_ = pionMass; mass2_ = pionMass;
    cut_minVtxRawMass_ = 0.35; cut_maxVtxRawMass_ = 0.65;
    cut_minVtxMass_ = 0.40; cut_maxVtxMass_ = 0.60;
  }
  else if ( vtxType == "phi" ) {
    pdgId_ = 333;
    pdgId1_ = pdgId2_ = 321;
    mass1_ = kaonMass; mass2_ = kaonMass;
    cut_minVtxRawMass_ = 0.96; cut_maxVtxRawMass_ = 1.08;
    cut_minVtxMass_ = 0.98; cut_maxVtxMass_ = 1.06;
  }
  else if ( vtxType == "lambda" ) {
    pdgId_ = 3122;
    pdgId1_ = 2212; pdgId2_ = 211;
    mass1_ = protonMass; mass2_ = pionMass;
    cut_minVtxRawMass_ = 1.04; cut_maxVtxRawMass_ = 1.24;
    cut_minVtxMass_ = 1.06; cut_maxVtxMass_ = 1.22;
  }
  else if ( vtxType == "d0" ) {
    pdgId_ = 421;
    pdgId1_ = 321; pdgId2_ = 211;
    mass1_ = kaonMass; mass2_ = pionMass;
    cut_minVtxRawMass_ = 1.6; cut_maxVtxRawMass_ = 2.1;
    cut_minVtxMass_ = 1.7; cut_maxVtxMass_ = 2.0;
  }
  else {
    pdgId_ = pset.getUntrackedParameter<int>("pdgId");
    mass1_ = pset.getUntrackedParameter<double>("mass1");
    mass2_ = pset.getUntrackedParameter<double>("mass2");
    pdgId1_ = pset.getUntrackedParameter<int>("pdgId1");
    pdgId2_ = pset.getUntrackedParameter<int>("pdgId2");
    cut_minVtxRawMass_ = pset.getUntrackedParameter<double>("minRawMass");
    cut_maxVtxRawMass_ = pset.getUntrackedParameter<double>("maxRawMass");
    cut_minVtxMass_ = pset.getUntrackedParameter<double>("minMass");
    cut_maxVtxMass_ = pset.getUntrackedParameter<double>("maxMass");
  }

  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");

  tree_->Branch("run", &b_run, "run/b"); // 8 bit unsigned integer
  tree_->Branch("lumi", &b_lumi, "lumi/b"); // 8 bit unsigned integer
  tree_->Branch("event", &b_event, "event/l"); // 64bit unsigned integer

  //tree_->Branch("genWeight", &b_genWeight, "genWeight/F");
  //tree_->Branch("puWeight", &b_puWeight, "puWeight/F");
  tree_->Branch("nPV", &b_nPV, "nPV/b");
  tree_->Branch("nSV", &b_nSV, "nSV/b");
  tree_->Branch("nGen", &b_nGen, "nGen/b");

  tree_->Branch("vtx_mass", &b_vtx_mass, "vtx_mass/F");
  tree_->Branch("vtx_pt"  , &b_vtx_pt  , "vtx_pt/F"  );
  tree_->Branch("vtx_lxy" , &b_vtx_lxy , "vtx_lxy/F" );
  tree_->Branch("vtx_vz"  , &b_vtx_vz  , "vtx_vz/F"  );
  b_trk_pdgId = new vint();
  b_trk_pt = new vfloat();
  b_trk_eta = new vfloat();
  b_trk_phi = new vfloat();
  tree_->Branch("trk_pdgId", "std::vector<int>", &b_trk_pdgId);
  tree_->Branch("trk_pt", "std::vector<float>", &b_trk_pt);
  tree_->Branch("trk_eta", "std::vector<float>", &b_trk_eta);
  tree_->Branch("trk_phi", "std::vector<float>", &b_trk_phi);
  b_mu_q = new vint();
  b_mu_id = new vint();
  b_mu_dR = new vfloat();
  b_mu_pt = new vfloat();
  tree_->Branch("mu_q", "std::vector<int>", &b_mu_q);
  tree_->Branch("mu_id", "std::vector<int>", &b_mu_id);
  tree_->Branch("mu_dR", "std::vector<float>", &b_mu_dR);
  tree_->Branch("mu_pt", "std::vector<float>", &b_mu_pt);

  tree_->Branch("gen_dR", &b_gen_dR, "gen_dR/F");

  hN_ = fs->make<TH1D>("hN", "hN", 100, 0, 100);
  hM_ = fs->make<TH1D>("hM", "hM", 100, cut_minVtxMass_, cut_maxVtxMass_);

  hN_->SetMinimum(0);
  hM_->SetMinimum(0);
}

void MuonMisIDNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  b_run = event.id().run();
  b_lumi = event.id().luminosityBlock();
  b_event = event.id().event();

  //b_genWeight = b_puWeight = 0;

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  edm::Handle<reco::TrackCollection> trackHandle;
  event.getByToken(trackToken_, trackHandle);

  edm::Handle<pat::PackedCandidateCollection> pfCandHandle;
  event.getByToken(pfCandToken_, pfCandHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);
  b_nPV = vertexHandle->size();
  const reco::Vertex pv = vertexHandle->at(0);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  event.getByToken(beamSpotToken_, beamSpotHandle);
  const reco::BeamSpot& beamSpot = *beamSpotHandle;

  const reco::Particle::Point pvPos = useBeamSpot_ ? beamSpot.position() : pv.position();
  const reco::Vertex::CovarianceMatrix pvCov = useBeamSpot_ ? beamSpot.covariance3D() : pv.covariance();

  edm::Handle<reco::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  b_nGen = 0;
  //std::vector<reco::GenParticle> genMuHads;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  std::vector<const reco::GenParticle*> resonances;
  if ( !event.isRealData() ) {
    //bool hasResonance = false;
    event.getByToken(genParticleToken_, genParticleHandle);
    for ( auto& p : *genParticleHandle ) {
      if ( !p.isLastCopy() ) continue;
      if ( abs(p.pdgId()) == pdgId_ ) {
        bool isDuplicated = false;
        for ( int i=0, n=p.numberOfDaughters(); i<n; ++i ) {
          if ( pdgId_ == std::abs(p.daughter(i)->pdgId()) ) { isDuplicated = true; break; }
        }
        if ( isDuplicated ) continue;
        resonances.push_back(&p);
      }

      //if ( p.status() != 1 ) continue;
      //const int aid = abs(p.pdgId());
      //if ( p.charge() != 0 and (aid != 13 or aid > 100) ) genMuHads.push_back(p);
    }
    if ( applyGenFilter_ and resonances.empty() ) return;
    b_nGen = resonances.size();
  }

  // Collect transient tracks
  std::vector<reco::TransientTrack> transTracks;
  if ( trackHandle.isValid() ) {
    for ( auto track = trackHandle->begin(); track != trackHandle->end(); ++track ) {
      if ( track->pt() < 0.35 or std::abs(track->eta()) > cut_maxTrkEta_ ) continue;
      // Apply basic track quality cuts
      if ( !track->quality(reco::TrackBase::loose) or
          track->normalizedChi2() >= cut_maxTrkChi2_ or track->numberOfValidHits() < cut_minTrkNHit_ ) continue;
      const double ipSigXY = std::abs(track->dxy(pvPos)/track->dxyError());
      const double ipSigZ = std::abs(track->dz(pvPos)/track->dzError());
      if ( ipSigXY < cut_minTrkSigXY_ or ipSigZ < cut_minTrkSigZ_  ) continue;
      auto transTrack = trackBuilder->build(&*track);
      transTracks.push_back(transTrack);
    }
  }
  else if ( pfCandHandle.isValid() ) {
    for ( auto cand = pfCandHandle->begin(); cand != pfCandHandle->end(); ++cand ) {
      if ( cand->pt() < 0.35 or std::abs(cand->eta()) > cut_maxTrkEta_ ) continue;
      auto track = cand->pseudoTrack();
      // Apply basic track quality cuts
      if ( !track.quality(reco::TrackBase::loose) or
            track.normalizedChi2() >= cut_maxTrkChi2_ or track.numberOfValidHits() < cut_minTrkNHit_ ) continue;
      const double ipSigXY = std::abs(track.dxy(pvPos)/track.dxyError());
      const double ipSigZ = std::abs(track.dz(pvPos)/track.dzError());
      if ( ipSigXY < cut_minTrkSigXY_ or ipSigZ < cut_minTrkSigZ_  ) continue;
      auto transTrack = trackBuilder->build(track);
      transTracks.push_back(transTrack);
    }
  }

  // Collect vertices after the fitting
  std::vector<SVFitResult> svs;
  const bool isSameFlav = (pdgId1_ == pdgId2_);
  for ( auto itr1 = transTracks.begin(); itr1 != transTracks.end(); ++itr1 ) {
    const reco::Track& track1 = itr1->track();
    if ( isSameFlav and track1.charge() < 0 ) continue;
    const double e1 = sqrt(mass1_*mass1_ + track1.momentum().mag2());

    for ( auto itr2 = transTracks.begin(); itr2 != transTracks.end(); ++itr2 ) {
      if ( itr1 == itr2 ) continue;
      const reco::Track& track2 = itr2->track();
      if ( track1.charge() == track2.charge() ) continue;
      if ( std::abs(deltaPhi(track1.phi(), track2.phi())) > 3.14 ) continue;
      if ( isSameFlav and track2.charge() > 0 ) continue;
      if ( track1.pt() < cut_minTrkPt_ and track2.pt() < cut_minTrkPt_ ) continue; // at least one track should pass minimum pt cut
      const double e2 = sqrt(mass2_*mass2_ + track2.momentum().mag2());

      const double px = track1.px() + track2.px();
      const double py = track1.py() + track2.py();
      const double pz = track1.pz() + track2.pz();
      const double p2 = px*px + py*py + pz*pz;
      const double e = e1+e2;

      const double rawMass = sqrt(e*e - p2);
      if ( rawMass < cut_minVtxRawMass_ or rawMass > cut_maxVtxRawMass_ ) continue;

      auto res = fitSV(pvPos, pvCov, *itr1, *itr2);
      if ( !res.isValid ) continue;

      svs.push_back(res);
    }
  }
  b_nSV = svs.size();
  hN_->Fill(svs.size());

  // Loop over the SV fit results to fill tree
  for ( const auto& sv : svs ) {
    b_vtx_mass = sv.p4.mass();
    b_vtx_pt = sv.p4.pt();
    b_vtx_lxy = sv.lxy;
    b_vtx_vz = sv.vz;

    // Fill track variables
    *b_trk_pdgId = {sv.q1*pdgId1_, sv.q2*pdgId2_};
    *b_trk_pt = {float(sv.leg1.pt()), float(sv.leg2.pt())};
    *b_trk_eta = {float(sv.leg1.eta()), float(sv.leg2.eta())};
    *b_trk_phi = {float(sv.leg1.phi()), float(sv.leg2.phi())};

    // Match muons to the SV legs
    *b_mu_q = {0, 0};
    *b_mu_id = {0, 0};
    *b_mu_dR = {-1, -1};
    *b_mu_pt = {0, 0};
    auto muonIdxPair = matchTwo(sv.leg1, sv.leg2, *muonHandle);
    const int muonIdx1 = muonIdxPair.first, muonIdx2 = muonIdxPair.second;
    if ( muonIdx1 >= 0 ) {
      const auto& mu = muonHandle->at(muonIdx1);
      b_mu_q->at(0) = mu.charge();
      b_mu_pt->at(0) = mu.pt();
      b_mu_id->at(0) = muonIdBit(mu, pv);
      b_mu_dR->at(0) = deltaR(mu.p4(), sv.leg1);
    }
    if ( muonIdx2 >= 0 ) {
      const auto& mu = muonHandle->at(muonIdx2);
      b_mu_q->at(1) = mu.charge();
      b_mu_pt->at(1) = mu.pt();
      b_mu_id->at(1) = muonIdBit(mu, pv);
      b_mu_dR->at(1) = deltaR(mu.p4(), sv.leg1);
    }

    // Match gen muons or hadrons to the SV legs
    double minDR = 0.3;
    const reco::GenParticle* matchedGen = 0;
    for ( auto& x : resonances ) {
      const double dR = deltaR(x->p4(), sv.p4);
      if ( dR < minDR ) { matchedGen = x; minDR = dR; }
    }
    b_gen_dR = matchedGen ? minDR : -1;
/*
    //b_gen1 = b_gen2 = LV();
    //b_genPdgId1 = b_genPdgId2 = b_genType1 = b_genType2 = 0;
    //b_genDR1 = b_genDR2 = -1;
    auto genIdxPair = matchTwo(sv.leg1, sv.leg2, genMuHads);
    const int genIdx1 = genIdxPair.first, genIdx2 = genIdxPair.second;
    if ( genIdx1 >= 0 ) {
      const auto& gp = genMuHads.at(genIdx1);
      b_gen1 = gp.p4();
      b_genPdgId1 = gp.pdgId();
      b_genType1 = genCategory(gp);
      b_genDR1 = deltaR(gp.p4(), sv.leg1);
    }
    if ( genIdx2 >= 0 ) {
      const auto& gp = genMuHads.at(genIdx2);
      b_gen2 = gp.p4();
      b_genPdgId2 = gp.pdgId();
      b_genType2 = genCategory(gp);
      b_genDR2 = deltaR(gp.p4(), sv.leg2);
    }
*/

    hM_->Fill(b_vtx_mass);
    tree_->Fill();
  }
}

SVFitResult MuonMisIDNtupleMaker::fitSV(const reco::Particle::Point& pvPos, const reco::Vertex::CovarianceMatrix& pvCov,
                                        const reco::TransientTrack& transTrack1,
                                        const reco::TransientTrack& transTrack2) const
{
  SVFitResult result;

  try {
    if ( !transTrack1.impactPointTSCP().isValid() or !transTrack2.impactPointTSCP().isValid() ) return result;
    auto ipState1 = transTrack1.impactPointTSCP().theState();
    auto ipState2 = transTrack2.impactPointTSCP().theState();
    //if ( std::abs(ipState1.position().z()-pv.z()) > 1 or
    //     std::abs(ipState2.position().z()-pv.z()) > 1 ) return SVFitResult();

    ClosestApproachInRPhi cApp;
    cApp.calculate(ipState1, ipState2);
    if ( !cApp.status() ) return result;

    const float dca = std::abs(cApp.distance());
    if ( dca < 0. or dca > cut_maxVtxDCA_ ) return result;

    GlobalPoint cxPt = cApp.crossingPoint();
    if ( std::hypot(cxPt.x(), cxPt.y()) > 120. or std::abs(cxPt.z()) > 300. ) return result;

    TrajectoryStateClosestToPoint caState1 = transTrack1.trajectoryStateClosestToPoint(cxPt);
    TrajectoryStateClosestToPoint caState2 = transTrack2.trajectoryStateClosestToPoint(cxPt);
    if ( !caState1.isValid() or !caState2.isValid() ) return result;
    if ( caState1.momentum().dot(caState2.momentum()) < 0 ) return result;

    std::vector<reco::TransientTrack> transTracks = {transTrack1, transTrack2};

    KalmanVertexFitter fitter(true);
    TransientVertex tsv = fitter.vertex(transTracks);
    if ( !tsv.isValid() or tsv.totalChiSquared() < 0. ) return result;

    reco::Vertex sv = tsv;
    if ( sv.normalizedChi2() > cut_maxVtxChi2_ ) return result;

    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
    typedef ROOT::Math::SVector<double, 3> SVector3;

    GlobalPoint vtxPos(sv.x(), sv.y(), sv.z());
    SMatrixSym3D totalCov = pvCov + sv.covariance();
    SVector3 distanceVectorXY(sv.x() - pvPos.x(), sv.y() - pvPos.y(), 0.);

    const double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
    const double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVectorXY)) / rVtxMag;
    if( rVtxMag < cut_minVtxLxy_ or rVtxMag > cut_maxVtxLxy_ or rVtxMag / sigmaRvtxMag < cut_minVtxSignif_ ) return result;

    //SVector3 distanceVector3D(sv.x() - pvx, sv.y() - pvy, sv.z() - pvz);
    //const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);

    // Cuts finished, now we create the candidates and push them back into the collections.
    int q1 = 0, q2 = 0;
    GlobalVector mom1, mom2;
    if ( !tsv.hasRefittedTracks() ) {
      q1 =  transTrack1.trajectoryStateClosestToPoint(vtxPos).charge();
      q2 =  transTrack2.trajectoryStateClosestToPoint(vtxPos).charge();
      mom1 = transTrack1.trajectoryStateClosestToPoint(vtxPos).momentum();
      mom2 = transTrack2.trajectoryStateClosestToPoint(vtxPos).momentum();
    }
    else {
      auto refTracks = tsv.refittedTracks();
      if ( refTracks.size() < 2 ) return result;

      q1 =  refTracks.at(0).trajectoryStateClosestToPoint(vtxPos).charge();
      q2 =  refTracks.at(1).trajectoryStateClosestToPoint(vtxPos).charge();
      mom1 = refTracks.at(0).trajectoryStateClosestToPoint(vtxPos).momentum();
      mom2 = refTracks.at(1).trajectoryStateClosestToPoint(vtxPos).momentum();
    }
    if ( mom1.mag() <= 0 or mom2.mag() <= 0 ) return result;
    const GlobalVector mom = mom1+mom2;

    const double candE1 = hypot(mom1.mag(), mass1_);
    const double candE2 = hypot(mom2.mag(), mass2_);
    const double vtxChi2 = sv.chi2();
    const double vtxNdof = sv.ndof();

    reco::Particle::Point vtx(sv.x(), sv.y(), sv.z());
    const reco::Vertex::CovarianceMatrix vtxCov(sv.covariance());

    const LV candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
    if ( cut_minVtxMass_ > candLVec.mass() or cut_maxVtxMass_ < candLVec.mass() ) return result;

    result.p4 = candLVec;
    result.vertex = vtx;
    result.leg1.SetXYZT(mom1.x(), mom1.y(), mom1.z(), candE1);
    result.leg2.SetXYZT(mom2.x(), mom2.y(), mom2.z(), candE2);
    result.q1 = q1;
    result.q2 = q2;
    result.chi2 = vtxChi2;
    result.ndof = vtxNdof;
    result.cov = vtxCov;
    result.lxy = rVtxMag;
    result.vz = std::abs(pvPos.z()-vtx.z());
    result.isValid = true;
  } catch ( std::exception& e ) { return SVFitResult(); }

  return result;
}

int MuonMisIDNtupleMaker::muonIdBit(const reco::Muon& mu, const reco::Vertex& vtx) const
{
  int result = 0;

  if ( muon::isLooseMuon(mu)       ) result |= 1<<0;
  if ( muon::isMediumMuon(mu)      ) result |= 1<<1;
  if ( muon::isTightMuon(mu, vtx)  ) result |= 1<<2;
  if ( muon::isSoftMuon(mu, vtx)   ) result |= 1<<3;

  return result;
}

int MuonMisIDNtupleMaker::genCategory(const reco::GenParticle& p) const
{
  return 0;
}

template<typename T>
std::pair<int, int> MuonMisIDNtupleMaker::matchTwo(const LV& lv1, const LV& lv2, T& coll) const
{
  std::map<double, int> idxs1, idxs2;
  for ( int i=0, n=coll.size(); i<n; ++i ) {
    const auto& p = coll.at(i);
    const double dR1 = deltaR(lv1, p.p4());
    const double dR2 = deltaR(lv2, p.p4());

    if ( dR1 < 0.3 ) idxs1[dR1] = i;
    if ( dR2 < 0.3 ) idxs2[dR2] = i;
  }
  int idx1 = idxs1.empty() ? -1 : idxs1.begin()->second;
  int idx2 = idxs2.empty() ? -1 : idxs2.begin()->second;
  // Special care for duplication
  if ( idx1 == idx2 and idx1 != -1 ) {
    const double dR1 = idxs1.begin()->first;
    const double dR2 = idxs2.begin()->first;
    if ( dR1 >= dR2 ) {
      if ( idxs1.size() > 1 ) idx1 = std::next(idxs1.begin())->second;
      else idx1 = -1;
    }
    else {
      if ( idxs2.size() > 1 ) idx2 = std::next(idxs2.begin())->second;
      else idx2 = -1;
    }
  }

  return make_pair(idx1, idx2);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonMisIDNtupleMaker);

