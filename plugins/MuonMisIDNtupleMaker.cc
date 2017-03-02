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
  reco::Candidate::LorentzVector p4, leg1, leg2, leg3;
  int q1, q2, q3;
  reco::Particle::Point vertex;
  reco::Vertex::CovarianceMatrix cov;
  double chi2, ndof, lxy, vz;
};

typedef math::XYZTLorentzVector LV;
typedef std::vector<int> vint;
typedef std::vector<float> vfloat;
typedef std::vector<bool> vbool;
typedef edm::ValueMap<float> VMapF;

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
  SVFitResult fitSV(const reco::Particle::Point& pvPos, const reco::Vertex::CovarianceMatrix& pvCov,
                    const reco::TransientTrack& transTrack1,
                    const reco::TransientTrack& transTrack2,
                    const reco::TransientTrack& transTrack3) const;
  void fillMuonId(reco::MuonRef muRef, const reco::Vertex& vertex, vbool* results[],
                  std::vector<edm::Handle<VMapF> >& vMapHandles, std::vector<vbool*>& results2, int idx) const;
  template <typename T1, typename T2>
  vint matchByDR(const T1& etas, const T1& phis, T2& coll) const;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfCandToken_;
  edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

  std::vector<edm::EDGetTokenT<VMapF> > vMapTokens_;

  // Constants for the vertex fit
  const double muonMass = 0.1057;
  const double pionMass = 0.1396;
  const double kaonMass = 0.4937;
  const double protonMass = 0.9383;

  const bool applyGenFilter_, useBeamSpot_;
  int pdgId_;

  double cut_minVtxRawMass_, cut_maxVtxRawMass_;
  double cut_minVtxRawMass12_, cut_maxVtxRawMass12_;
  double cut_minVtxMass12_, cut_maxVtxMass12_;
  double cut_minVtxMass_, cut_maxVtxMass_, cut_minVtxLxy_, cut_maxVtxLxy_, cut_minVtxLxyz_, cut_maxVtxLxyz_;
  const double cut_minTrkPt_, cut_maxTrkEta_;
  const int cut_minTrkNHit_;
  const double cut_maxTrkChi2_, cut_minTrkSigXY_, cut_minTrkSigZ_;
  const double cut_maxVtxDCA_, cut_maxVtxChi2_, cut_minVtxSignif_;
  const double cut_minVtxSignif3D_;

  double mass1_, mass2_, mass3_;
  int pdgId1_, pdgId2_, pdgId3_;

  // Trees and histograms
  TTree* tree_;
  unsigned char b_run, b_lumi;
  unsigned long long int b_event;
  //double b_genWeight, b_puWeight;
  unsigned char b_nPV, b_nSV, b_nGen;

  float b_vtx_mass, b_vtx_pt, b_vtx_lxy, b_vtx_vz, b_vtx_chi2, b_vtx_ndof;
  float b_vtx_mass12, b_vtx_mass23, b_vtx_mass13; // 3 track case for Dalitz
  vfloat* b_trk_pt, * b_trk_eta, * b_trk_phi;
  vint* b_trk_pdgId;

  vfloat* b_mu_pt, * b_mu_dR;
  vint* b_mu_q;
  static const int nId_ = 26;
  vbool* b_mu_id[nId_];
  std::vector<vbool*> b_mu_vids_;

  float b_gen_dR;

  TH1D* hN_;
  TH1D* hM_;
};

MuonMisIDNtupleMaker::MuonMisIDNtupleMaker(const edm::ParameterSet& pset):
  applyGenFilter_(pset.getUntrackedParameter<bool>("applyGenFilter", false)),
  useBeamSpot_(pset.getUntrackedParameter<bool>("useBeamSpot", false)),
  cut_minVtxLxy_(pset.getUntrackedParameter<double>("minVtxLxy")),
  cut_maxVtxLxy_(pset.getUntrackedParameter<double>("maxVtxLxy")),
  cut_minVtxLxyz_(pset.getUntrackedParameter<double>("minVtxLxyz")),
  cut_maxVtxLxyz_(pset.getUntrackedParameter<double>("maxVtxLxyz")),
  cut_minTrkPt_(pset.getUntrackedParameter<double>("minTrkPt")),
  cut_maxTrkEta_(pset.getUntrackedParameter<double>("maxTrkEta")),
  cut_minTrkNHit_(pset.getUntrackedParameter<int>("minTrkNHit")),
  cut_maxTrkChi2_(pset.getUntrackedParameter<double>("maxTrkChi2")),
  cut_minTrkSigXY_(pset.getUntrackedParameter<double>("minTrkSigXY")),
  cut_minTrkSigZ_(pset.getUntrackedParameter<double>("minTrkSigZ")),
  cut_maxVtxDCA_(pset.getUntrackedParameter<double>("maxVtxDCA")),
  cut_maxVtxChi2_(pset.getUntrackedParameter<double>("maxVtxChi2")),
  cut_minVtxSignif_(pset.getUntrackedParameter<double>("minVtxSignif")),
  cut_minVtxSignif3D_(pset.getUntrackedParameter<double>("minVtxSignif3D"))
{
  genParticleToken_ = consumes<reco::GenParticleCollection>(pset.getParameter<edm::InputTag>("genParticles"));
  vertexToken_ = consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"));
  beamSpotToken_ = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  trackToken_ = consumes<reco::TrackCollection>(pset.getParameter<edm::InputTag>("tracks"));
  pfCandToken_ = consumes<pat::PackedCandidateCollection>(edm::InputTag("pfCandidates"));
  muonToken_ = consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons"));

  const string vtxType = pset.getUntrackedParameter<string>("vtxType");
  pdgId3_ = 0;
  if ( vtxType == "kshort" ) {
    pdgId_ = 310;
    pdgId1_ = pdgId2_ = 211;
    mass1_ = pionMass; mass2_ = pionMass;
    cut_minVtxRawMass_ = 0.43; cut_maxVtxRawMass_ = 0.58;
    cut_minVtxMass_ = 0.45; cut_maxVtxMass_ = 0.56;
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
    cut_minVtxRawMass_ = 1.09; cut_maxVtxRawMass_ = 1.14;
    cut_minVtxMass_ = 1.1; cut_maxVtxMass_ = 1.13;
  }
  else if ( vtxType == "jpsi" ) {
    pdgId_ = 443;
    pdgId1_ = pdgId2_ = 13;
    mass1_ = mass2_ = muonMass;
    cut_minVtxRawMass_ = 3.09-0.2; cut_maxVtxRawMass_ = 3.09+0.2;
    cut_minVtxMass_ = 3.09-0.1; cut_maxVtxMass_ = 3.09+0.1;
  }
  else if ( vtxType == "D0" ) {
    pdgId_ = 421;
    pdgId1_ = 321; pdgId2_ = 211;
    mass1_ = kaonMass; mass2_ = pionMass;
    cut_minVtxRawMass_ = 1.6; cut_maxVtxRawMass_ = 2.1;
    cut_minVtxMass_ = 1.7; cut_maxVtxMass_ = 2.0;
  }
  else if ( vtxType == "D+" ) {
    pdgId_ = 421;
    pdgId1_ = 321; pdgId2_ = 211; pdgId3_ = 211;
    mass1_ = kaonMass; mass2_ = pionMass; mass3_ = pionMass;
    cut_minVtxRawMass_ = 1.6; cut_maxVtxRawMass_ = 2.1;
    cut_minVtxRawMass12_ = 0; cut_maxVtxRawMass12_ = 999; // No cut on kpi since the grouping is not real
    cut_minVtxMass12_ = 0; cut_maxVtxMass12_ = 999; // No cut on kpi since the grouping is not real
    cut_minVtxMass_ = 1.7; cut_maxVtxMass_ = 2.0;
  }
  else if ( vtxType == "B+" ) {
    pdgId_ = 521;
    pdgId1_ = pdgId2_ = 13; pdgId3_ = 321;
    mass1_ = mass2_ = muonMass; mass3_ = kaonMass;
    cut_minVtxRawMass_ = 5.0; cut_maxVtxRawMass_ = 5.5;
    cut_minVtxRawMass12_ = 3.09-0.2; cut_maxVtxRawMass12_ = 3.09+0.2; // Cut on jpsi
    cut_minVtxMass12_ = 3.09-0.1; cut_maxVtxMass12_ = 3.09+0.1; // Cut on jpsi
    cut_minVtxMass_ = 5.1; cut_maxVtxMass_ = 5.4;
  }
  else {
    pdgId_ = pset.getUntrackedParameter<int>("pdgId");
    mass1_ = pset.getUntrackedParameter<double>("mass1");
    mass2_ = pset.getUntrackedParameter<double>("mass2");
    mass3_ = pset.getUntrackedParameter<double>("mass3", 0);
    pdgId1_ = pset.getUntrackedParameter<int>("pdgId1");
    pdgId2_ = pset.getUntrackedParameter<int>("pdgId2");
    pdgId3_ = pset.getUntrackedParameter<int>("pdgId2", 0);
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
  tree_->Branch("vtx_ndof", &b_vtx_ndof, "vtx_ndof/F");
  tree_->Branch("vtx_chi2", &b_vtx_chi2, "vtx_chi2/F");
  if ( pdgId3_ != 0 ) {
    tree_->Branch("vtx_mass12", &b_vtx_mass12, "vtx_mass12/F");
    tree_->Branch("vtx_mass23", &b_vtx_mass23, "vtx_mass23/F");
    tree_->Branch("vtx_mass13", &b_vtx_mass13, "vtx_mass13/F");
  }

  b_trk_pdgId = new vint();
  b_trk_pt = new vfloat();
  b_trk_eta = new vfloat();
  b_trk_phi = new vfloat();
  tree_->Branch("trk_pdgId", "std::vector<int>", &b_trk_pdgId);
  tree_->Branch("trk_pt", "std::vector<float>", &b_trk_pt);
  tree_->Branch("trk_eta", "std::vector<float>", &b_trk_eta);
  tree_->Branch("trk_phi", "std::vector<float>", &b_trk_phi);
  b_mu_q = new vint();
  b_mu_dR = new vfloat();
  b_mu_pt = new vfloat();
  tree_->Branch("mu_q", "std::vector<int>", &b_mu_q);
  tree_->Branch("mu_dR", "std::vector<float>", &b_mu_dR);
  tree_->Branch("mu_pt", "std::vector<float>", &b_mu_pt);

  for ( int i=0; i<nId_; ++i ) b_mu_id[i] = new vbool();

  int nId = 0;
  tree_->Branch("mu_isTight" , "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isMedium", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isLoose" , "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isSoft"  , "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isHighPt", "std::vector<bool>", &b_mu_id[nId++]);

  tree_->Branch("mu_isGLB", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTRK", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isSTA", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isRPC", "std::vector<bool>", &b_mu_id[nId++]);

  tree_->Branch("mu_isGLBPT", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTMLastLoose", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTMLastTight", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTM2DLoose", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTM2DTight", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isOneLoose", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isOneTight", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isLastLowPtLoose", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isLastLowPtTight", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isGMTkChi2Compat", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isGMStaChi2Compat", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isGMTkKinkTight", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTMLastAngLoose", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTMLastAngTight", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTMOneAngLoose", "std::vector<bool>", &b_mu_id[nId++]);
  tree_->Branch("mu_isTMOneAngTight", "std::vector<bool>", &b_mu_id[nId++]);

  tree_->Branch("mu_isRPCMuLoose", "std::vector<bool>", &b_mu_id[nId++]);

  assert(nId == nId_);

  for ( auto idMapSet : pset.getParameter<std::vector<edm::ParameterSet> >("idMaps") ) {
    const auto name = idMapSet.getUntrackedParameter<std::string>("name");
    vMapTokens_.push_back(consumes<VMapF>(idMapSet.getParameter<edm::InputTag>("src")));
    b_mu_vids_.push_back(new vbool());
    tree_->Branch(("mu_"+name).c_str(), "std::vector<bool>", b_mu_vids_.back());
  }

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

  std::vector<edm::Handle<VMapF> > vMapHandles;
  for ( auto& token : vMapTokens_ ) {
    edm::Handle<VMapF> handle;
    event.getByToken(token, handle);
    vMapHandles.push_back(handle);
  }

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
  std::vector<bool> isMuonSeededFlags;
  if ( trackHandle.isValid() ) {
    for ( auto track = trackHandle->begin(); track != trackHandle->end(); ++track ) {
      if ( track->pt() < 0.5 or std::abs(track->eta()) > cut_maxTrkEta_ ) continue;
      // Apply basic track quality cuts
      if ( !track->quality(reco::TrackBase::loose) or
          track->normalizedChi2() >= cut_maxTrkChi2_ or track->numberOfValidHits() < cut_minTrkNHit_ ) continue;
      const double ipSigXY = std::abs(track->dxy(pvPos)/track->dxyError());
      const double ipSigZ = std::abs(track->dz(pvPos)/track->dzError());
      if ( ipSigXY < cut_minTrkSigXY_ or ipSigZ < cut_minTrkSigZ_  ) continue;
      auto transTrack = trackBuilder->build(&*track);
      transTracks.push_back(transTrack);
      const bool isMuonSeeded = track->originalAlgo() != reco::TrackBase::muonSeededStepOutIn; // Avoid bias from muon seeded track
      isMuonSeededFlags.push_back(isMuonSeeded);
    }
  }
  else if ( pfCandHandle.isValid() ) {
    for ( auto cand = pfCandHandle->begin(); cand != pfCandHandle->end(); ++cand ) {
      if ( cand->pt() < 0.5 or std::abs(cand->eta()) > cut_maxTrkEta_ ) continue;
      auto track = cand->pseudoTrack();
      // Apply basic track quality cuts
      if ( !track.quality(reco::TrackBase::loose) or
            track.normalizedChi2() >= cut_maxTrkChi2_ or track.numberOfValidHits() < cut_minTrkNHit_ ) continue;
      const double ipSigXY = std::abs(track.dxy(pvPos)/track.dxyError());
      const double ipSigZ = std::abs(track.dz(pvPos)/track.dzError());
      if ( ipSigXY < cut_minTrkSigXY_ or ipSigZ < cut_minTrkSigZ_  ) continue;
      auto transTrack = trackBuilder->build(track);
      transTracks.push_back(transTrack);
      const bool isMuonSeeded = track.originalAlgo() != reco::TrackBase::muonSeededStepOutIn; // Avoid bias from muon seeded track
      isMuonSeededFlags.push_back(isMuonSeeded);
    }
  }

  // Collect vertices after the fitting
  std::vector<SVFitResult> svs;
  const bool isSameFlav = (pdgId1_ == pdgId2_);
  for ( auto itr1 = transTracks.begin(); itr1 != transTracks.end(); ++itr1 ) {
    if ( std::abs(pdgId1_) != 13 and isMuonSeededFlags.at(itr1-transTracks.begin()) ) continue;
    const reco::Track& track1 = itr1->track();
    if ( isSameFlav and track1.charge() < 0 ) continue;
    const double e1 = sqrt(mass1_*mass1_ + track1.momentum().mag2());

    for ( auto itr2 = transTracks.begin(); itr2 != transTracks.end(); ++itr2 ) {
      if ( std::abs(pdgId2_) != 13 and isMuonSeededFlags.at(itr2-transTracks.begin()) ) continue;
      if ( itr1 == itr2 ) continue;
      const reco::Track& track2 = itr2->track();
      if ( track1.charge() == track2.charge() ) continue;
      if ( isSameFlav ) {
        if ( track2.charge() > 0 ) continue;
      }
      else {
        // heavier particle takes larger portion of momentum
        if ( mass1_ > mass2_ ) {
          if ( track1.pt() < track2.pt() ) continue;
        }
        else {
          if ( track2.pt() < track1.pt() ) continue;
        }
      }
      const double e2 = sqrt(mass2_*mass2_ + track2.momentum().mag2());

      const auto mom = track1.momentum() + track2.momentum();
      const double e = e1+e2;
      const double rawMass = sqrt(e*e - mom.mag2());

      if ( pdgId3_ == 0 ) { // V->TT case
        if ( track1.pt() < cut_minTrkPt_ and track2.pt() < cut_minTrkPt_ ) continue; // at least one track should pass minimum pt cut
        if ( rawMass < cut_minVtxRawMass_ or rawMass > cut_maxVtxRawMass_ ) continue;

        auto res = fitSV(pvPos, pvCov, *itr1, *itr2);
        if ( !res.isValid ) continue;

        svs.push_back(res);
      }
      else { // V->TTT three body decay case
        // Works only for no-tertiory vertex cases. ex) D+->Kpipi or B+->JpsiK+
        // Assume track1 and track2 are opposite signed : D+ -> (K,pi),pi and B+ -> (mu,mu),K

        // Do the kinematic preselection on the sub-vertex
        if ( rawMass < cut_minVtxRawMass12_ or rawMass > cut_maxVtxRawMass12_ ) continue;

        for ( auto itr3 = transTracks.begin(); itr3 != transTracks.end(); ++itr3 ) {
          if ( std::abs(pdgId3_) != 13 and isMuonSeededFlags.at(itr3-transTracks.begin()) ) continue;
          if ( itr1 == itr3 or itr2 == itr3 ) continue;
          const reco::Track& track3 = itr3->track();
          if ( track1.pt() < cut_minTrkPt_ and track2.pt() < cut_minTrkPt_ and track3.pt() < cut_minTrkPt_ ) continue;

          const auto momF = mom + track3.momentum();
          const double e3 = sqrt(mass3_*mass3_ + track3.momentum().mag2());
          const double eF = e+e3;

          const double rawMassF = sqrt(eF*eF - momF.mag2());
          if ( rawMassF < cut_minVtxRawMass_ or rawMassF > cut_maxVtxRawMass_ ) continue;

          auto res = fitSV(pvPos, pvCov, *itr1, *itr2, *itr3);
          if ( !res.isValid ) continue;

          svs.push_back(res);
        }
      }
    }
  }
  b_nSV = svs.size();
  hN_->Fill(svs.size());

  // Loop over the SV fit results to fill tree
  for ( const auto& sv : svs ) {
    b_vtx_mass = sv.p4.mass();
    b_vtx_pt = sv.p4.pt();
    b_vtx_ndof = sv.ndof;
    b_vtx_chi2 = sv.chi2;
    b_vtx_lxy = sv.lxy;
    b_vtx_vz = sv.vz;

    if ( pdgId3_ != 0 ) {
      b_vtx_mass12 = (sv.leg1+sv.leg2).mass();
      b_vtx_mass23 = (sv.leg2+sv.leg3).mass();
      b_vtx_mass13 = (sv.leg1+sv.leg3).mass();
    }

    // Fill track variables
    *b_trk_pdgId = {sv.q1*pdgId1_, sv.q2*pdgId2_};
    *b_trk_pt = {float(sv.leg1.pt()), float(sv.leg2.pt())};
    *b_trk_eta = {float(sv.leg1.eta()), float(sv.leg2.eta())};
    *b_trk_phi = {float(sv.leg1.phi()), float(sv.leg2.phi())};
    if ( pdgId3_ != 0 ) {
      b_trk_pdgId->push_back(sv.q3*pdgId3_);
      b_trk_pt->push_back(sv.leg3.pt());
      b_trk_eta->push_back(sv.leg3.eta());
      b_trk_phi->push_back(sv.leg3.phi());
    }

    // Match muons to the SV legs
    auto muMatch = matchByDR(*b_trk_eta, *b_trk_phi, *muonHandle);
    if ( pdgId3_ == 0 ) {
      *b_mu_q = vint({0, 0});
      *b_mu_dR = vfloat({-1, -1});
      *b_mu_pt = vfloat({0, 0});
      for ( int i=0; i<nId_; ++i ) *b_mu_id[i] = vbool({false, false});
      for ( int i=0, n=vMapHandles.size(); i<n; ++i ) *b_mu_vids_[i] = vbool({false, false});
    }
    else {
      *b_mu_q = vint({0, 0, 0});
      *b_mu_dR = vfloat({-1, -1, -1});
      *b_mu_pt = vfloat({0, 0, 0});
      for ( int i=0; i<nId_; ++i ) *b_mu_id[i] = vbool({false, false, false});
      for ( int i=0, n=vMapHandles.size(); i<n; ++i ) *b_mu_vids_[i] = vbool({false, false, false});
    }
    if ( muMatch[0] >= 0 ) {
      reco::MuonRef muRef(muonHandle, muMatch[0]);
      b_mu_q->at(0) = muRef->charge();
      b_mu_pt->at(0) = muRef->pt();
      b_mu_dR->at(0) = deltaR(muRef->p4(), sv.leg1);
      fillMuonId(muRef, pv, b_mu_id, vMapHandles, b_mu_vids_, 0);
    }
    if ( muMatch[1] >= 0 ) {
      reco::MuonRef muRef(muonHandle, muMatch[1]);
      b_mu_q->at(1) = muRef->charge();
      b_mu_pt->at(1) = muRef->pt();
      b_mu_dR->at(1) = deltaR(muRef->p4(), sv.leg1);
      fillMuonId(muRef, pv, b_mu_id, vMapHandles, b_mu_vids_, 1);
    }
    if ( pdgId3_ != 0 and muMatch[2] >= 0 ) {
      reco::MuonRef muRef(muonHandle, muMatch[2]);
      b_mu_q->at(2) = muRef->charge();
      b_mu_pt->at(2) = muRef->pt();
      b_mu_dR->at(2) = deltaR(muRef->p4(), sv.leg1);
      fillMuonId(muRef, pv, b_mu_id, vMapHandles, b_mu_vids_, 2);
    }

    // Match gen muons or hadrons to the SV legs
    double minDR = 0.3;
    const reco::GenParticle* matchedGen = 0;
    for ( auto& x : resonances ) {
      const double dR = deltaR(x->p4(), sv.p4);
      if ( dR < minDR ) { matchedGen = x; minDR = dR; }
    }
    b_gen_dR = matchedGen ? minDR : -1;

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
    if ( std::hypot(cxPt.x(), cxPt.y()) > 120. or std::abs(cxPt.z()) > 300. ) return result; // CA in the tracker volume, 1.2m, 3m
 
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

    SVector3 distanceVector3D(sv.x() - pvPos.x(), sv.y() - pvPos.y(), sv.z() - pvPos.z());
    const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
    const double sigmaRvtxMag3D = sqrt(ROOT::Math::Similarity(totalCov, distanceVector3D)) / rVtxMag3D;
    if( rVtxMag3D < cut_minVtxLxyz_ or rVtxMag3D > cut_maxVtxLxyz_ or rVtxMag3D / sigmaRvtxMag3D < cut_minVtxSignif3D_ ) return result;

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
    if ( mom1.mag() <= 1 or mom2.mag() <= 1 ) return result;
    const GlobalVector mom = mom1+mom2;

    const double candE1 = hypot(mom1.mag(), mass1_);
    const double candE2 = hypot(mom2.mag(), mass2_);
    const double vtxChi2 = sv.chi2();
    const double vtxNdof = sv.ndof();

    reco::Particle::Point vtx(sv.x(), sv.y(), sv.z());
    if ( sv.x()*mom.x()+sv.y()*mom.y()+sv.z()*mom.z() <= 0 ) return result;

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

SVFitResult MuonMisIDNtupleMaker::fitSV(const reco::Particle::Point& pvPos, const reco::Vertex::CovarianceMatrix& pvCov,
                                        const reco::TransientTrack& transTrack1,
                                        const reco::TransientTrack& transTrack2,
                                        const reco::TransientTrack& transTrack3) const
{
  SVFitResult result;

  try {
    if ( !transTrack1.impactPointTSCP().isValid() or
         !transTrack2.impactPointTSCP().isValid() or
         !transTrack3.impactPointTSCP().isValid() ) return result;
    auto ipState1 = transTrack1.impactPointTSCP().theState();
    auto ipState2 = transTrack2.impactPointTSCP().theState();
    auto ipState3 = transTrack3.impactPointTSCP().theState();
    //if ( std::abs(ipState1.position().z()-pv.z()) > 1 or
    //     std::abs(ipState2.position().z()-pv.z()) > 1 ) return SVFitResult();

    ClosestApproachInRPhi cApp12, cApp23, cApp13;
    cApp12.calculate(ipState1, ipState2);
    cApp23.calculate(ipState2, ipState3);
    cApp13.calculate(ipState1, ipState3);
    if ( !cApp12.status() or !cApp23.status() or !cApp13.status() ) return result;

    const float dca12 = std::abs(cApp12.distance());
    const float dca23 = std::abs(cApp23.distance());
    const float dca13 = std::abs(cApp13.distance());
    if ( dca12 < 0. or dca12 > cut_maxVtxDCA_ or
         dca23 < 0. or dca23 > cut_maxVtxDCA_ or
         dca13 < 0. or dca13 > cut_maxVtxDCA_ ) return result;

    const GlobalPoint cxPt12 = cApp12.crossingPoint();
    const GlobalPoint cxPt23 = cApp23.crossingPoint();
    const GlobalPoint cxPt13 = cApp13.crossingPoint();
    if ( std::hypot(cxPt12.x(), cxPt12.y()) > 120. or std::abs(cxPt12.z()) > 300. or
         std::hypot(cxPt23.x(), cxPt23.y()) > 120. or std::abs(cxPt23.z()) > 300. or
         std::hypot(cxPt13.x(), cxPt13.y()) > 120. or std::abs(cxPt13.z()) > 300. ) return result;
    // FIXME: approximate crossing point as an average of 3 crossing points
    const GlobalPoint cxPt((cxPt12.x()+cxPt23.x()+cxPt13.x())/3.,
                           (cxPt12.y()+cxPt23.y()+cxPt13.y())/3.,
                           (cxPt12.z()+cxPt23.z()+cxPt13.z())/3.);

    TrajectoryStateClosestToPoint caState1 = transTrack1.trajectoryStateClosestToPoint(cxPt);
    TrajectoryStateClosestToPoint caState2 = transTrack2.trajectoryStateClosestToPoint(cxPt);
    TrajectoryStateClosestToPoint caState3 = transTrack3.trajectoryStateClosestToPoint(cxPt);
    if ( !caState1.isValid() or !caState2.isValid() or !caState3.isValid() ) return result;
    if ( caState1.momentum().dot(caState2.momentum()) < 0 or
         caState1.momentum().dot(caState3.momentum()) < 0 or
         caState2.momentum().dot(caState3.momentum()) < 0 ) return result;

    std::vector<reco::TransientTrack> transTracks = {transTrack1, transTrack2, transTrack3};

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

    SVector3 distanceVector3D(sv.x() - pvPos.x(), sv.y() - pvPos.y(), sv.z() - pvPos.z());
    const double rVtxMag3D = ROOT::Math::Mag(distanceVector3D);
    const double sigmaRvtxMag3D = sqrt(ROOT::Math::Similarity(totalCov, distanceVector3D)) / rVtxMag3D;
    if( rVtxMag3D < cut_minVtxLxyz_ or rVtxMag3D > cut_maxVtxLxyz_ or rVtxMag3D / sigmaRvtxMag3D < cut_minVtxSignif3D_ ) return result;

    // Cuts finished, now we create the candidates and push them back into the collections.
    int q1 = 0, q2 = 0, q3 = 0;
    GlobalVector mom1, mom2, mom3;
    if ( !tsv.hasRefittedTracks() ) {
      q1 =  transTrack1.trajectoryStateClosestToPoint(vtxPos).charge();
      q2 =  transTrack2.trajectoryStateClosestToPoint(vtxPos).charge();
      q3 =  transTrack3.trajectoryStateClosestToPoint(vtxPos).charge();
      mom1 = transTrack1.trajectoryStateClosestToPoint(vtxPos).momentum();
      mom2 = transTrack2.trajectoryStateClosestToPoint(vtxPos).momentum();
      mom3 = transTrack3.trajectoryStateClosestToPoint(vtxPos).momentum();
    }
    else {
      auto refTracks = tsv.refittedTracks();
      if ( refTracks.size() < 3 ) return result;

      q1 =  refTracks.at(0).trajectoryStateClosestToPoint(vtxPos).charge();
      q2 =  refTracks.at(1).trajectoryStateClosestToPoint(vtxPos).charge();
      q3 =  refTracks.at(2).trajectoryStateClosestToPoint(vtxPos).charge();
      mom1 = refTracks.at(0).trajectoryStateClosestToPoint(vtxPos).momentum();
      mom2 = refTracks.at(1).trajectoryStateClosestToPoint(vtxPos).momentum();
      mom3 = refTracks.at(2).trajectoryStateClosestToPoint(vtxPos).momentum();
    }
    if ( mom1.mag() <= 1 or mom2.mag() <= 1 or mom3.mag() <= 1 ) return result;
    const GlobalVector mom12 = mom1+mom2;
    const GlobalVector mom = mom12+mom3;

    const double candE1 = hypot(mom1.mag(), mass1_);
    const double candE2 = hypot(mom2.mag(), mass2_);
    const double candE3 = hypot(mom3.mag(), mass3_);
    const double vtxChi2 = sv.chi2();
    const double vtxNdof = sv.ndof();

    reco::Particle::Point vtx(sv.x(), sv.y(), sv.z());
    if ( sv.x()*mom.x()+sv.y()*mom.y()+sv.z()*mom.z() <= 0 ) return result;
    const reco::Vertex::CovarianceMatrix vtxCov(sv.covariance());

    const LV cand12LVec(mom12.x(), mom12.y(), mom12.z(), candE1+candE2);
    if ( cut_minVtxMass12_ > cand12LVec.mass() or cut_maxVtxMass12_ < cand12LVec.mass() ) return result;
    const LV candLVec(mom.x(), mom.y(), mom.z(), candE1+candE2+candE3);
    if ( cut_minVtxMass_ > candLVec.mass() or cut_maxVtxMass_ < candLVec.mass() ) return result;

    result.p4 = candLVec;
    result.vertex = vtx;
    result.leg1.SetXYZT(mom1.x(), mom1.y(), mom1.z(), candE1);
    result.leg2.SetXYZT(mom2.x(), mom2.y(), mom2.z(), candE2);
    result.leg3.SetXYZT(mom3.x(), mom3.y(), mom3.z(), candE3);
    result.q1 = q1;
    result.q2 = q2;
    result.q3 = q3;
    result.chi2 = vtxChi2;
    result.ndof = vtxNdof;
    result.cov = vtxCov;
    result.lxy = rVtxMag;
    result.vz = std::abs(pvPos.z()-vtx.z());
    result.isValid = true;
  } catch ( std::exception& e ) { return SVFitResult(); }

  return result;
}

void MuonMisIDNtupleMaker::fillMuonId(reco::MuonRef muRef, const reco::Vertex& vtx, vbool* result[],
                                      std::vector<edm::Handle<VMapF> >& vMapHandles, std::vector<vbool*>& results2, int idx) const
{
  const auto& mu = *muRef;

  int nId = 0;
  result[nId++]->at(idx) = muon::isTightMuon(mu, vtx);
  result[nId++]->at(idx) = muon::isMediumMuon(mu);
  result[nId++]->at(idx) = muon::isLooseMuon(mu);
  result[nId++]->at(idx) = muon::isSoftMuon(mu, vtx);
  result[nId++]->at(idx) = muon::isHighPtMuon(mu, vtx);

  result[nId++]->at(idx) = mu.isGlobalMuon();
  result[nId++]->at(idx) = mu.isTrackerMuon();
  result[nId++]->at(idx) = mu.isStandAloneMuon();
  result[nId++]->at(idx) = mu.isRPCMuon();

  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::GlobalMuonPromptTight);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMLastStationLoose);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMLastStationTight);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TM2DCompatibilityLoose);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TM2DCompatibilityTight);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMOneStationLoose);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMOneStationTight);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMLastStationOptimizedLowPtLoose);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMLastStationOptimizedLowPtTight);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::GMTkChiCompatibility);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::GMStaChiCompatibility);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::GMTkKinkTight);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMLastStationAngLoose);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMLastStationAngTight);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMOneStationAngLoose);
  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::TMOneStationAngTight);

  result[nId++]->at(idx) = muon::isGoodMuon(mu, muon::RPCMuLoose, reco::Muon::RPCHitAndTrackArbitration);

  assert(nId == nId_);

  for ( int i=0, n=vMapHandles.size(); i<n; ++i ) {
    auto& handle = vMapHandles.at(i);
    results2[i]->at(idx) = (*handle)[muRef] > 0.5;
  }
}

template<typename T1, typename T2>
vint MuonMisIDNtupleMaker::matchByDR(const T1& etas, const T1& phis, T2& coll) const
{
  assert(etas.size() == phis.size());
  const int n = etas.size(), m = coll.size();

  std::map<double, std::pair<int, int> > matches;
  for ( int i=0; i<n; ++i ) {
    const auto eta = etas[i], phi = phis[i];
    for ( int j=0; j<m; ++j ) {
      const auto innerTrack = coll.at(j).innerTrack();
      const auto bestTrack = coll.at(j).bestTrackRef();
      const reco::TrackBase* trk = nullptr;
      if ( innerTrack.isNonnull() ) trk = innerTrack.get();
      else if ( bestTrack.isNonnull() ) trk = bestTrack.get(); // Retry with best track if inner track is invalid
      else continue; // This should not happen, but just for safety...

      const double trkPhi = trk->phi();
      const double trkEta = trk->eta();
      const double dR2 = deltaR2(trkEta, trkPhi, eta, phi);
      matches[dR2] = make_pair(i, j);
    }
  }

  // Cleanup to make unique mapping
  vint result(n, -1);
  for ( const auto& match : matches ) {
    //const double dR2 = match.first;
    const int i = match.second.first, j = match.second.second;
    if ( result[i] == -1 ) result[i] = j;
  }

  return result;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonMisIDNtupleMaker);

