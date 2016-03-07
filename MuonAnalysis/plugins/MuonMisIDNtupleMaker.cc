#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class MuonMisIDNtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  MuonMisIDNtupleMaker(const edm::ParameterSet& pset);
  virtual ~MuonMisIDNtupleMaker() {};
  void analyze(const edm::Event& event, const edm::EventSetup&) override;

private:
  //const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muToken_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> svToken_;
  edm::EDGetTokenT<std::vector<double> > lxyToken_;

  void initTree();
  TTree* tree_;
  int run_, lumi_, event_;
  double genWeight_, puWeight_;
  int nPV_;

  double mass_, pt_, lxy_;

  int muQ1_, muQ2_, trackQ1_, trackQ2_, pdgId1_, pdgId2_;
  math::XYZTLorentzVector mu1_, mu2_, track1_, track2_;

  TH1F* hM_;
};

MuonMisIDNtupleMaker::MuonMisIDNtupleMaker(const edm::ParameterSet& pset)
{
  muToken_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("mu"));
  auto svLabel = pset.getParameter<edm::InputTag>("sv");
  svToken_ = consumes<reco::VertexCompositeCandidateCollection>(svLabel);
  lxyToken_ = consumes<std::vector<double> >(edm::InputTag(svLabel.label(), "lxy"));

  usesResource("TFileService");

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "tree");

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/I");

  tree_->Branch("genWeight", &genWeight_, "genWeight/D");
  tree_->Branch("puWeight", &puWeight_, "puWeight/D");
  tree_->Branch("nPV", &nPV_, "nPV/I");

  tree_->Branch("mass", &mass_, "mass/D");
  tree_->Branch("pt"  , &pt_  , "pt/D"  );
  tree_->Branch("lxy" , &lxy_ , "lxy/D" );

  tree_->Branch("muQ1", &muQ1_, "muQ1/I");
  tree_->Branch("muQ2", &muQ2_, "muQ2/I");
  tree_->Branch("trackQ1", &trackQ1_, "trackQ1/I");
  tree_->Branch("trackQ2", &trackQ2_, "trackQ2/I");
  tree_->Branch("pdgId1", &pdgId1_, "pdgId1/I");
  tree_->Branch("pdgId2", &pdgId2_, "pdgId2/I");

  hM_ = fs->make<TH1F>("hM", "hM", 100, 0.4, 0.6);
}

void MuonMisIDNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup&)
{
  run_ = event.id().run();
  lumi_ = event.id().luminosityBlock();
  event_ = event.id().event();

  edm::Handle<edm::View<reco::Muon> > muHandle;
  event.getByToken(muToken_, muHandle);

  edm::Handle<reco::VertexCompositeCandidateCollection> svHandle;
  event.getByToken(svToken_, svHandle);

  edm::Handle<std::vector<double> > lxyHandle;
  event.getByToken(lxyToken_, lxyHandle);

  for ( int i=0, n=svHandle->size(); i<n; ++i ) {
    //const double lxy = lxyHandle->at(i);
    auto& sv = svHandle->at(i);
    hM_->Fill(sv.mass());
  }

  tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonMisIDNtupleMaker);

