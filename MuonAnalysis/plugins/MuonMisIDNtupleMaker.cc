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
  double chi2, ndof, lxy;
};

typedef math::XYZTLorentzVector LV;

class MuonMisIDNtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  MuonMisIDNtupleMaker(const edm::ParameterSet& pset);
  virtual ~MuonMisIDNtupleMaker() {};
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  SVFitResult fitSV(const reco::Vertex& pv,
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

  // Constants for the vertex fit
  const double pionMass = 0.1396;
  const double kaonMass = 0.4937;
  const double protonMass = 0.9383;

  const double trkMinPt_, trkMaxEta_;
  double vtxMinRawMass_, vtxMaxRawMass_;
  double vtxMinMass_, vtxMaxMass_, vtxMinLxy_, vtxMaxLxy_;
  const double trkChi2_, trkDCA_, vtxChi2_, vtxSignif_;

  double mass1_, mass2_;
  int pdgId1_, pdgId2_;

  // Trees and histograms
  TTree* tree_;
  int b_run, b_lumi, b_event;
  double b_genWeight, b_puWeight;
  int b_nPV;

  double b_mass, b_pt, b_lxy;
  LV b_track1, b_track2;

  int b_muQ1, b_muQ2, b_pdgId1, b_pdgId2;
  LV b_mu1, b_mu2;
  int b_muId1, b_muId2;
  double b_muDR1, b_muDR2;

  int b_genPdgId1, b_genPdgId2, b_genType1, b_genType2;
  LV b_gen1, b_gen2;
  double b_genDR1, b_genDR2;

  TH1D* hN_;
  TH1D* hM_, * hMAll_;
};

MuonMisIDNtupleMaker::MuonMisIDNtupleMaker(const edm::ParameterSet& pset):
  trkMinPt_(pset.getParameter<double>("trkMinPt")),
  trkMaxEta_(pset.getParameter<double>("trkMaxEta")),
  vtxMinLxy_(pset.getParameter<double>("vtxMinLxy")),
  vtxMaxLxy_(pset.getParameter<double>("vtxMaxLxy")),
  trkChi2_(pset.getParameter<double>("trkChi2")),
  trkDCA_(pset.getParameter<double>("trkDCA")),
  vtxChi2_(pset.getParameter<double>("vtxChi2")),
  vtxSignif_(pset.getParameter<double>("vtxSignif"))
{
  genParticleToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticles"));
  vertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"));
  trackToken_ = consumes<reco::TrackCollection>(pset.getParameter<edm::InputTag>("tracks"));
  pfCandToken_ = consumes<pat::PackedCandidateCollection>(pset.getParameter<edm::InputTag>("pfCandidates"));
  muonToken_ = consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));

  const string vtxType = pset.getParameter<string>("vtxType");
  if ( vtxType == "kshort" ) {
    pdgId1_ = pdgId2_ = 211;
    mass1_ = pionMass; mass2_ = pionMass;
    vtxMinRawMass_ = 0.40; vtxMaxRawMass_ = 0.60;
    vtxMinMass_ = 0.43; vtxMaxMass_ = 0.57;
  }
  else if ( vtxType == "phi" ) {
    pdgId1_ = pdgId2_ = 321;
    mass1_ = kaonMass; mass2_ = kaonMass;
    vtxMinRawMass_ = 0.95; vtxMaxRawMass_ = 1.08;
    vtxMinMass_ = 0.98; vtxMaxMass_ = 1.06;
  }
  else if ( vtxType == "lambda" ) {
    pdgId1_ = 2212; pdgId2_ = 211;
    mass1_ = protonMass; mass2_ = pionMass;
    vtxMinRawMass_ = 1.091; vtxMaxRawMass_ = 1.139;
    vtxMinMass_ = 1.101; vtxMaxMass_ = 1.129;
  }
  else {
    mass1_ = pset.getParameter<double>("mass1");
    mass2_ = pset.getParameter<double>("mass2");
    pdgId1_ = pset.getParameter<int>("pdgId1");
    pdgId2_ = pset.getParameter<int>("pdgId2");
    vtxMinRawMass_ = pset.getParameter<double>("minRawMass");
    vtxMaxRawMass_ = pset.getParameter<double>("maxRawMass");
    vtxMinMass_ = pset.getParameter<double>("minMass");
    vtxMaxMass_ = pset.getParameter<double>("maxMass");
  }

  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");

  tree_->Branch("run", &b_run, "run/I");
  tree_->Branch("lumi", &b_lumi, "lumi/I");
  tree_->Branch("event", &b_event, "event/I");

  tree_->Branch("genWeight", &b_genWeight, "genWeight/D");
  tree_->Branch("puWeight", &b_puWeight, "puWeight/D");
  tree_->Branch("nPV", &b_nPV, "nPV/I");

  tree_->Branch("mass", &b_mass, "mass/D");
  tree_->Branch("pt"  , &b_pt  , "pt/D"  );
  tree_->Branch("lxy" , &b_lxy , "lxy/D" );
  tree_->Branch("pdgId1", &b_pdgId1, "pdgId1/I");
  tree_->Branch("pdgId2", &b_pdgId2, "pdgId2/I");
  tree_->Branch("track1", "math::XYZTLorentzVector", &b_track1);
  tree_->Branch("track2", "math::XYZTLorentzVector", &b_track2);

  tree_->Branch("muQ1", &b_muQ1, "muQ1/I");
  tree_->Branch("muQ2", &b_muQ2, "muQ2/I");
  tree_->Branch("muId1", &b_muId1, "muId1/I");
  tree_->Branch("muId2", &b_muId2, "muId2/I");
  tree_->Branch("muDR1", &b_muDR1, "muDR1/D");
  tree_->Branch("muDR2", &b_muDR2, "muDR2/D");
  tree_->Branch("mu1", "math::XYZTLorentzVector", &b_mu1);
  tree_->Branch("mu2", "math::XYZTLorentzVector", &b_mu2);

  tree_->Branch("genPdgId1", &b_genPdgId1, "genPdgId1/I");
  tree_->Branch("genPdgId2", &b_genPdgId2, "genPdgId2/I");
  tree_->Branch("genType1", &b_genType1, "genType1/I");
  tree_->Branch("genType2", &b_genType2, "genType2/I");
  tree_->Branch("genDR1", &b_genDR1, "genDR1/D");
  tree_->Branch("genDR2", &b_genDR2, "genDR2/D");
  tree_->Branch("gen1", "math::XYZTLorentzVector", &b_gen1);
  tree_->Branch("gen2", "math::XYZTLorentzVector", &b_gen2);

  hN_ = fs->make<TH1D>("hN", "hN", 100, 0, 100);
  hMAll_ = fs->make<TH1D>("hMAll", "hMAll", 100, vtxMinMass_, vtxMaxMass_);
  hM_ = fs->make<TH1D>("hM", "hM", 100, vtxMinMass_, vtxMaxMass_);
}

void MuonMisIDNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  b_run = event.id().run();
  b_lumi = event.id().luminosityBlock();
  b_event = event.id().event();

  b_genWeight = b_puWeight = -999;
  b_nPV = -999;
  b_mass = b_pt = b_lxy = -999;
  b_pdgId1 = b_pdgId2 = -999;
  b_muQ1 = b_muQ2 = b_muId1 = b_muId2 = b_muDR1 = b_muDR2 = -999;
  b_genPdgId1 = b_genPdgId2 = b_genType1 = b_genType2 = b_genDR1 = b_genDR2 = -999;

  b_track1 = b_track2 = b_mu1 = b_mu2 = b_gen1 = b_gen2 = LV();

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  edm::Handle<reco::TrackCollection> trackHandle;
  event.getByToken(trackToken_, trackHandle);

  edm::Handle<pat::PackedCandidateCollection> pfCandHandle;
  event.getByToken(pfCandToken_, pfCandHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);
  reco::Vertex pv = vertexHandle->at(0);

  edm::Handle<reco::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  std::vector<reco::GenParticle> genMuHads;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  if ( !event.isRealData() ) {
    event.getByToken(genParticleToken_, genParticleHandle);
    for ( auto& p : *genParticleHandle ) {
      if ( p.status() != 1 ) continue;
      const int aid = abs(p.pdgId());
      if ( p.charge() != 0 and (aid != 13 or aid > 100) ) genMuHads.push_back(p);
    }
  }

  // Collect transient tracks
  std::vector<reco::TransientTrack> transTracks;
  if ( trackHandle.isValid() ) {
    for ( auto track = trackHandle->begin(); track != trackHandle->end(); ++track ) {
      if ( track->pt() < trkMinPt_ or std::abs(track->eta()) > trkMaxEta_ ) continue;
      // Apply basic track quality cuts
      if ( !track->quality(reco::TrackBase::loose) or
          track->normalizedChi2() >= 5 or track->numberOfValidHits() < 6 ) continue;
      auto transTrack = trackBuilder->build(&*track);
      transTracks.push_back(transTrack);
    }
  }
  else if ( pfCandHandle.isValid() ) {
    for ( auto cand = pfCandHandle->begin(); cand != pfCandHandle->end(); ++cand ) {
      if ( cand->pt() < trkMinPt_ or std::abs(cand->eta()) > trkMaxEta_ ) continue;

      auto track = cand->pseudoTrack();
      auto transTrack = trackBuilder->build(track);
      transTracks.push_back(transTrack);
    }
  }

  // Collect vertices after the fitting
  std::vector<SVFitResult> svs;
  for ( auto itr1 = transTracks.begin(); itr1 != transTracks.end(); ++itr1 ) {
    const reco::Track& track1 = itr1->track();
    if ( pdgId1_ == pdgId2_ and track1.charge() < 0 ) continue;
    const double e1 = sqrt(mass1_*mass1_ + track1.momentum().mag2());

    for ( auto itr2 = transTracks.begin(); itr2 != transTracks.end(); ++itr2 ) {
      if ( itr1 == itr2 ) continue;
      const reco::Track& track2 = itr2->track();
      if ( pdgId1_ == pdgId2_ and track2.charge() > 0 ) continue;
      const double e2 = sqrt(mass2_*mass2_ + track2.momentum().mag2());

      const double px = track1.px() + track2.px();
      const double py = track1.px() + track2.px();
      const double pz = track1.px() + track2.px();
      const double p2 = px*px + py*py + pz*pz;
      const double e = e1+e2;

      const double rawMass = sqrt(e*e - p2);
      if ( rawMass < vtxMinRawMass_ or rawMass > vtxMaxRawMass_ ) continue;

      auto res = fitSV(pv, *itr1, *itr2);
      if ( !res.isValid ) continue;

      svs.push_back(res);
    }
  }
  hN_->Fill(svs.size());

  // Loop over the SV fit results and find best one in the event
  SVFitResult sv;
  for ( auto& isv : svs ) {
    hMAll_->Fill(isv.p4.mass());

    if ( !sv.isValid or isv.p4.pt() > sv.p4.pt() ) {
      sv = isv;
    }
  }
  if ( !sv.isValid ) return;

  // Fill variables for the best SV
  b_mass = sv.p4.mass();
  b_pt = sv.p4.pt();
  b_lxy = sv.lxy;
  b_pdgId1 = sv.q1*pdgId1_;
  b_pdgId2 = sv.q2*pdgId2_;
  b_track1 = sv.leg1;
  b_track2 = sv.leg2;
  hM_->Fill(b_mass);

  // Match muons to the SV legs
  auto muonIdxPair = matchTwo(sv.leg1, sv.leg2, *muonHandle);
  const int muonIdx1 = muonIdxPair.first, muonIdx2 = muonIdxPair.second;
  if ( muonIdx1 >= 0 ) {
    const auto& mu = muonHandle->at(muonIdx1);
    b_muQ1 = mu.charge();
    b_mu1 = mu.p4();
    b_muId1 = muonIdBit(mu, pv);
    b_muDR1 = deltaR(mu.p4(), sv.leg1);
  }
  if ( muonIdx2 >= 0 ) {
    const auto& mu = muonHandle->at(muonIdx2);
    b_muQ2 = mu.charge();
    b_mu2 = mu.p4();
    b_muId2 = muonIdBit(mu, pv);
    b_muDR2 = deltaR(mu.p4(), sv.leg2);
  }

  // Match gen muons or hadrons to the SV legs
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

  tree_->Fill();
}

SVFitResult MuonMisIDNtupleMaker::fitSV(const reco::Vertex& pv,
                                        const reco::TransientTrack& transTrack1,
                                        const reco::TransientTrack& transTrack2) const
{
  SVFitResult result;

  try {
    auto ipState1 = transTrack1.impactPointTSCP().theState();
    auto ipState2 = transTrack2.impactPointTSCP().theState();

    ClosestApproachInRPhi cApp;
    cApp.calculate(ipState1, ipState2);
    if ( !cApp.status() ) return result;

    const float dca = std::abs(cApp.distance());
    if ( dca < 0. or dca > trkDCA_ ) return result;

    GlobalPoint cxPt = cApp.crossingPoint();
    if ( std::hypot(cxPt.x(), cxPt.y()) > 120. or std::abs(cxPt.z()) > 300. ) return result;

    TrajectoryStateClosestToPoint caState1 = transTrack1.trajectoryStateClosestToPoint(cxPt);
    TrajectoryStateClosestToPoint caState2 = transTrack2.trajectoryStateClosestToPoint(cxPt);
    if ( !caState1.isValid() or !caState2.isValid() ) return result;

    std::vector<reco::TransientTrack> transTracks = {transTrack1, transTrack2};

    KalmanVertexFitter fitter(true);
    TransientVertex tsv = fitter.vertex(transTracks);
    if ( !tsv.isValid() or tsv.totalChiSquared() < 0. ) return result;

    reco::Vertex sv = tsv;
    if ( sv.normalizedChi2() > vtxChi2_ ) return result;

    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
    typedef ROOT::Math::SVector<double, 3> SVector3;

    GlobalPoint vtxPos(sv.x(), sv.y(), sv.z());
    SMatrixSym3D totalCov = pv.covariance() + sv.covariance();
    SVector3 distanceVectorXY(sv.x() - pv.x(), sv.y() - pv.y(), 0.);

    const double rVtxMag = ROOT::Math::Mag(distanceVectorXY);
    const double sigmaRvtxMag = sqrt(ROOT::Math::Similarity(totalCov, distanceVectorXY)) / rVtxMag;
    if( rVtxMag < vtxMinLxy_ or rVtxMag > vtxMaxLxy_ or rVtxMag / sigmaRvtxMag < vtxSignif_ ) return result;

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
    if ( vtxMinMass_ > candLVec.mass() or vtxMaxMass_ < candLVec.mass() ) return result;

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
    result.isValid = true;
  } catch ( std::exception& e ) { return result; }

  return result;
}

int MuonMisIDNtupleMaker::muonIdBit(const reco::Muon& mu, const reco::Vertex& vtx) const
{
  int result = 0;

  if ( muon::isLooseMuon(mu)       ) result |= 1<<0;
  if ( muon::isMediumMuon(mu)      ) result |= 1<<1;
  if ( muon::isTightMuon(mu, vtx)  ) result |= 1<<2;

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

