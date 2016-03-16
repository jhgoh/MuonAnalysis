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
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

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

  //const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
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

  int b_muQ1, b_muQ2, b_pdgId1, b_pdgId2;
  math::XYZTLorentzVector b_mu1, b_mu2, b_track1, b_track2;
  int b_muId1, b_muId2;
  double b_muDR1, b_muDR2;

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
  vertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"));
  trackToken_ = consumes<reco::TrackCollection>(pset.getParameter<edm::InputTag>("tracks"));
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
  b_muQ1 = b_muQ2 = b_muId1 = b_muId2 = -999;
  b_muDR1 = b_muDR2 = -999;

  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);

  edm::Handle<reco::TrackCollection> trackHandle;
  event.getByToken(trackToken_, trackHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  event.getByToken(vertexToken_, vertexHandle);
  reco::Vertex pv = vertexHandle->at(0);

  edm::Handle<reco::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  // Collect transient tracks
  std::vector<reco::TransientTrack> transTracks;
  for ( auto track = trackHandle->begin(); track != trackHandle->end(); ++track ) {
    if ( track->pt() < trkMinPt_ or std::abs(track->eta()) > trkMaxEta_ ) continue;
    // Apply basic track quality cuts
    if ( !track->quality(reco::TrackBase::loose) or
         track->normalizedChi2() >= 5 or track->numberOfValidHits() < 6 ) continue;
    auto transTrack = trackBuilder->build(&*track);
    transTracks.push_back(transTrack);
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
  hM_->Fill(b_mass);

  // Match muons to the SV legs
  std::map<double, reco::MuonRef> muRefs1, muRefs2;
  for ( int i=0, n=muonHandle->size(); i<n; ++i ) {
    const auto& mu = muonHandle->at(i);
    const double dR1 = deltaR(sv.leg1, mu.p4());
    const double dR2 = deltaR(sv.leg2, mu.p4());

    if ( dR1 < 0.3 ) muRefs1[dR1] = reco::MuonRef(muonHandle, i);
    if ( dR2 < 0.3 ) muRefs2[dR2] = reco::MuonRef(muonHandle, i);
  }
  reco::MuonRef muRef1, muRef2;
  if      ( !muRefs1.empty() ) muRef1 = muRefs1.begin()->second;
  else if ( !muRefs2.empty() ) muRef2 = muRefs2.begin()->second;
  // Special care for duplication
  if ( muRef1.isNonnull() and muRef1 == muRef2 ) {
    const double dR1 = muRefs1.begin()->first;
    const double dR2 = muRefs2.begin()->first;
    if ( dR1 > dR2 ) {
      if ( muRefs1.size() > 1 ) muRef1 = std::next(muRefs1.begin())->second;
      else muRef1 = reco::MuonRef();
    }
    else {
      if ( muRefs2.size() > 1 ) muRef2 = std::next(muRefs2.begin())->second;
      else muRef2 = reco::MuonRef();
    }
  }

  if ( muRef1.isNonnull() ) {
    b_muQ1 = muRef1->charge();
    b_mu1 = muRef1->p4();
    b_muId1 = muonIdBit(*muRef1, pv);
    b_muDR1 = deltaR(muRef1->p4(), sv.leg1);
  }
  if ( muRef2.isNonnull() ) {
    b_muQ2 = muRef2->charge();
    b_mu2 = muRef2->p4();
    b_muId2 = muonIdBit(*muRef2, pv);
    b_muDR2 = deltaR(muRef2->p4(), sv.leg2);
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

    const math::XYZTLorentzVector candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2);
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
  if ( muon::isMediumMuon(mu) ) result |= 1<<1;
  if ( muon::isTightMuon(mu, vtx)  ) result |= 1<<2;

  return result;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonMisIDNtupleMaker);

