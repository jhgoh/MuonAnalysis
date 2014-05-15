#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class MuonRPCRecHitAnalyzer : public edm::EDAnalyzer
{
public:
  MuonRPCRecHitAnalyzer( const edm::ParameterSet& pset);
  ~MuonRPCRecHitAnalyzer() {};

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  edm::InputTag muonLabel_;
  //edm::InputTag rpcRecHitsLabel_;

  typedef std::vector<double> doubles;
  typedef std::vector<int> ints;
  typedef doubles* doublesP;
  typedef ints* intsP;

  TH1F* hEvent_;
  TTree* tree_;
  int run_, lumi_, event_;
  doublesP pt_, eta_, phi_;
  intsP charge_;
  intsP nHit_, nPixelHit_, nTrackerHit_, nMuonHit_, nRPCHit_;

};

MuonRPCRecHitAnalyzer::MuonRPCRecHitAnalyzer(const edm::ParameterSet& pset)
{
  muonLabel_ = pset.getParameter<edm::InputTag>("muon");

  pt_ = new doubles;
  eta_ = new doubles;
  phi_ = new doubles;
  charge_ = new ints;

  nHit_        = new ints;
  nPixelHit_   = new ints;
  nTrackerHit_ = new ints;
  nMuonHit_    = new ints;
  nRPCHit_     = new ints;

  edm::Service<TFileService> fs;
  hEvent_ = fs->make<TH1F>("hEvent", "Event count", 5, 0, 5);
  hEvent_->GetXaxis()->SetBinLabel(1, "Total");
  hEvent_->GetXaxis()->SetBinLabel(2, "Muon");
  hEvent_->GetXaxis()->SetBinLabel(3, "RPCHit");

  tree_ = fs->make<TTree>("tree", "tree");
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("event", &event_, "event/I");
  tree_->Branch("lumi" , &lumi_ , "lumi/I" );

  tree_->Branch("pt" , &pt_ );
  tree_->Branch("eta", &eta_);
  tree_->Branch("phi", &phi_);
  tree_->Branch("charge", &charge_);

  tree_->Branch("nHit"       , &nHit_       );
  tree_->Branch("nPixelHit"  , &nPixelHit_  );
  tree_->Branch("nTrackerHit", &nTrackerHit_);
  tree_->Branch("nMuonHit"   , &nMuonHit_   );
  tree_->Branch("nRPCHit"    , &nRPCHit_    );
}

void MuonRPCRecHitAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  run_ = event.id().run();
  event_ = event.id().event();
  lumi_ = event.id().luminosityBlock();
  pt_->clear();
  eta_->clear();
  phi_->clear();
  charge_->clear();
  nHit_->clear();
  nPixelHit_->clear();
  nTrackerHit_->clear();
  nMuonHit_->clear();
  nRPCHit_->clear();

  hEvent_->Fill(0);

  edm::Handle<edm::View<reco::Muon> > muonHandle;
  event.getByLabel(muonLabel_, muonHandle);

  int nMuon = 0;
  for ( unsigned int i=0, n=muonHandle->size(); i<n; ++i )
  {
    const reco::Muon& mu = muonHandle->at(i);
    if ( !mu.isGlobalMuon() ) continue;
    const reco::Track& muTrack = *mu.globalTrack();
    ++nMuon;

    pt_->push_back(mu.pt());
    eta_->push_back(mu.eta());
    phi_->push_back(mu.phi());
    charge_->push_back(mu.charge());
    nHit_->push_back(muTrack.hitPattern().numberOfValidHits());
    nPixelHit_->push_back(muTrack.hitPattern().numberOfValidPixelHits());
    nTrackerHit_->push_back(muTrack.hitPattern().numberOfValidTrackerHits());
    nMuonHit_->push_back(muTrack.hitPattern().numberOfValidMuonHits());
    nRPCHit_->push_back(muTrack.hitPattern().numberOfValidMuonRPCHits());
  }
  if ( nMuon > 0 ) hEvent_->Fill(1);

  tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonRPCRecHitAnalyzer);


