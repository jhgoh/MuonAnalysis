#ifndef BESTAnalysis_TTbarLeptonJet_EventTupleProducerData_H
#define BESTAnalysis_TTbarLeptonJet_EventTupleProducerData_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <vector>
#include <string>

template<typename Lepton>
class EventTupleProducerData : public edm::EDAnalyzer
{
public:
  EventTupleProducerData(const edm::ParameterSet& pset);
  ~EventTupleProducerData();

  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup);

private:
  bool hasMother(const reco::Candidate* p, const int pdgId);

private:
  bool doMCMatch_;

  // Input objects
  edm::InputTag genLabel_;
  edm::InputTag leptonLabel_;
  edm::InputTag jetLabel_;
  edm::InputTag metLabel_;
  edm::InputTag vertexLabel_;
  std::string bTagType_;

  // Cuts
  StringCutObjectSelector<Lepton, true>* isGoodLepton_;
  StringCutObjectSelector<pat::Jet, true>* isGoodJet_;

  // Output tree
  TTree* tree_;
  int run_, lumi_, event_;
  math::XYZTLorentzVector lepton_, met_;
  int charge_;
  int nPV_;
  std::vector<math::XYZTLorentzVector>* jets_;
  std::vector<double>* bTags_;
  std::vector<int>* jetMCBits_;

};

template<typename Lepton>
EventTupleProducerData<Lepton>::EventTupleProducerData(const edm::ParameterSet& pset)
{
  doMCMatch_ = pset.getParameter<bool>("doMCMatch");

  // Input labels
  genLabel_ = pset.getParameter<edm::InputTag>("gen");
  jetLabel_ = pset.getParameter<edm::InputTag>("jet");
  metLabel_ = pset.getParameter<edm::InputTag>("met");
  leptonLabel_ = pset.getParameter<edm::InputTag>("lepton");
  vertexLabel_ = pset.getParameter<edm::InputTag>("vertex");

  std::string leptonCut = pset.getParameter<std::string>("leptonCut");
  isGoodLepton_ = new StringCutObjectSelector<Lepton, true>(leptonCut);
  std::string jetCut = pset.getParameter<std::string>("jetCut");
  isGoodJet_ = new StringCutObjectSelector<pat::Jet, true>(jetCut);

  bTagType_ = pset.getParameter<std::string>("bTagType");

  // Output histograms and tree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree", "Mixed event tree");
  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("lumi", &lumi_, "lumi/I");
  tree_->Branch("event", &event_, "event/I");

  tree_->Branch("lepton", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", &lepton_);
  tree_->Branch("charge", &charge_, "charge/I");
  tree_->Branch("met", "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>", &met_);
  tree_->Branch("vertex", &nPV_, "nPV_/I");
  

  jets_ = new std::vector<math::XYZTLorentzVector>();
  bTags_ = new std::vector<double>();
  jetMCBits_ = new std::vector<int>();
//  eventWeight_ = new std::vector<double>();
  tree_->Branch("jets", "std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >", &jets_);
  tree_->Branch("bTags", "std::vector<double>", &bTags_);
  tree_->Branch("jetMCBits", "std::vector<int>", &jetMCBits_);
}

template<typename Lepton>
EventTupleProducerData<Lepton>::~EventTupleProducerData()
{
  if ( jets_ ) delete jets_;
  if ( bTags_ ) delete bTags_;
  if ( jetMCBits_ ) delete jetMCBits_;
}

template<typename Lepton>
void EventTupleProducerData<Lepton>::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  jets_->clear();
  bTags_->clear();
  jetMCBits_->clear();

  edm::Handle<edm::View<reco::Vertex> > vertexHandle;
  event.getByLabel(vertexLabel_, vertexHandle);
  nPV_ = vertexHandle->size();

  edm::Handle<edm::View<Lepton> > leptonHandle;
  event.getByLabel(leptonLabel_, leptonHandle);
  int nPassingLepton = 0;
  int leptonIdx = -1;
  for ( int i=0, n=leptonHandle->size(); i<n; ++i )
  {
    if ( !(*isGoodLepton_)(leptonHandle->at(i)) ) continue;

    ++nPassingLepton;
    leptonIdx = i;
  }
  if ( nPassingLepton != 1 ) return;
  lepton_ = leptonHandle->at(leptonIdx).p4();
  charge_ = leptonHandle->at(leptonIdx).charge();

  edm::Handle<edm::View<pat::MET> > metHandle;
  event.getByLabel(metLabel_, metHandle);
  met_ = metHandle->at(0).p4();

  const reco::Candidate* genLepB, * genHadB, * genHadJ1, * genHadJ2;
  genLepB = genHadB = genHadJ1 = genHadJ2 = 0;
  if ( doMCMatch_ )
  {
    edm::Handle<reco::GenParticleCollection> genHandle;
    event.getByLabel(genLabel_, genHandle);
    if ( !genHandle.isValid() ) 
    {
      doMCMatch_ = false;
    }
    else
    {
      // Find top quark from the genParticles
      const reco::GenParticle* t1 = 0, * t2 = 0;
      const reco::GenParticle* w1 = 0, * w2 = 0;
      const reco::GenParticle* b1 = 0, * b2 = 0;
      const reco::GenParticle* l1 = 0, * l2 = 0;
      std::vector<const reco::GenParticle*> w1Decay, w2Decay;
      for ( int i=0, n=genHandle->size(); i<n; ++i )
      {
        const reco::GenParticle& p = genHandle->at(i);
        if ( p.status() != 3 ) continue;

        switch(p.pdgId())
        {
          case   6: t1 = &p; break;
          case  -6: t2 = &p; break;
          case  24: w1 = &p; break;
          case -24: w2 = &p; break;
          case   5: b1 = &p; break;
          case  -5: b2 = &p; break;
          case  15: break;
          case -15: break;
          case -11: case -13: l1 = &p; break;
          case  11: case  13: l2 = &p; break;
          case  12: case  14: case  16:
          case -12: case -14: case -16: break;
          default:
            if ( hasMother(&p, 24) ) w1Decay.push_back(&p);
            else if ( hasMother(&p, -24) ) w2Decay.push_back(&p);
        }
      }

      if ( !t1 or !t2 or !w1 or !w2 or !b1 or !b2 ) return;
      if ( (l1 and l2) or (!l1 and !l2) ) return;

      if ( l1 and w2Decay.size() > 1 )
      {
        genLepB = b1;
        genHadB = b2;

        genHadJ1 = w2Decay[0];
        genHadJ2 = w2Decay[1];
      }
      else if ( l2 and w1Decay.size() > 1 )
      {
        genLepB = b2;
        genHadB = b1;

        genHadJ1 = w1Decay[0];
        genHadJ2 = w1Decay[1];
      }
      else
      {
        cout << "FATAL: This should not happen\n";
        return;
      }
    }
  }

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  event.getByLabel(jetLabel_, jetHandle);
  for ( int i=0, n=jetHandle->size(); i<n; ++i )
  {
    const pat::Jet& jet = jetHandle->at(i);

    if ( !(*isGoodJet_)(jet) ) continue;

    jets_->push_back(jet.p4());
    bTags_->push_back(jet.bDiscriminator(bTagType_));

    int jetMCBit = 0;
    const double dRLepB = genLepB ? deltaR(jet, *genLepB) : 1e9;
    const double dRHadB = genHadB ? deltaR(jet, *genHadB) : 1e9;
    const double dRHadJ1 = genHadJ1 ? deltaR(jet, *genHadJ1) : 1e9;
    const double dRHadJ2 = genHadJ2 ? deltaR(jet, *genHadJ2) : 1e9;

    if ( dRLepB < 0.5 ) jetMCBit |= 1;
    if ( dRHadB < 0.5 ) jetMCBit |= 2;
    if ( dRHadJ1 < 0.5 or dRHadJ2 < 0.5 ) jetMCBit |= 4;

    jetMCBits_->push_back(jetMCBit);

  }
  if ( jets_->size() < 3 ) return;

  // Now put jets in current event to the event cache
  run_ = event.run();
  lumi_ = event.luminosityBlock();
  event_ = event.id().event();

  tree_->Fill();
}

template<typename Lepton>
bool EventTupleProducerData<Lepton>::hasMother(const reco::Candidate* p, const int pdgId)
{
  for ( int i=0, n=p->numberOfMothers(); i<n; ++i )
  {
    const reco::Candidate* m = p->mother(i);
    if ( !m ) return false;
    else if ( m->pdgId() == pdgId ) return true;
    else if ( hasMother(m, pdgId) ) return true;
  }

  return false;
}

typedef EventTupleProducerData<pat::Electron> EventTupleProducerDataElectron;
typedef EventTupleProducerData<pat::Muon> EventTupleProducerDataMuon;

DEFINE_FWK_MODULE(EventTupleProducerDataElectron);
DEFINE_FWK_MODULE(EventTupleProducerDataMuon);

#endif
