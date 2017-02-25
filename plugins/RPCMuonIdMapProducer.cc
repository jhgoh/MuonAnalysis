#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

//#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class RPCMuonIdMapProducer : public edm::stream::EDProducer<>
{
public:
  RPCMuonIdMapProducer(const edm::ParameterSet& pset);
  virtual ~RPCMuonIdMapProducer() {};
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  typedef std::vector<pat::Muon> PatMuonColl;
  typedef std::vector<reco::Muon> RecoMuonColl;
  edm::EDGetTokenT<PatMuonColl> patMuonToken_;
  edm::EDGetTokenT<RecoMuonColl> recoMuonToken_;

  typedef std::vector<float> vfloat;
};

RPCMuonIdMapProducer::RPCMuonIdMapProducer(const edm::ParameterSet& pset):
  patMuonToken_(consumes<PatMuonColl>(pset.getParameter<edm::InputTag>("src"))),
  recoMuonToken_(consumes<RecoMuonColl>(pset.getParameter<edm::InputTag>("src")))
{
  produces<edm::ValueMap<float> >("Loose");
  produces<edm::ValueMap<float> >("Tight");
  produces<edm::ValueMap<float> >("TwoStationLoose");
  produces<edm::ValueMap<float> >("TwoStationTight");
  produces<edm::ValueMap<float> >("LastStationLoose");
  produces<edm::ValueMap<float> >("LastStationTight");
  produces<edm::ValueMap<float> >("SecondStationLoose");
  produces<edm::ValueMap<float> >("SecondStationTight");
}

void RPCMuonIdMapProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<PatMuonColl> patMuonHandle;
  edm::Handle<RecoMuonColl> recoMuonHandle;

  const bool isPatMuon = event.getByToken(patMuonToken_, patMuonHandle);
  //const bool isRecoMuon = isPatMuon ? false : event.getByToken(recoMuonToken_, recoMuonHandle);
  if ( !isPatMuon ) event.getByToken(recoMuonToken_, recoMuonHandle);

  std::map<std::string, vfloat> idMaps = {
    {"Loose", vfloat()},
    {"Tight", vfloat()},
    {"TwoStationLoose", vfloat()},
    {"TwoStationTight", vfloat()},
    {"LastStationLoose", vfloat()},
    {"LastStationTight", vfloat()},
    {"SecondStationLoose", vfloat()},
    {"SecondStationTight", vfloat()},
  };

  const int nMuon = isPatMuon ? patMuonHandle->size() : recoMuonHandle->size();
  for ( int i=0; i<nMuon; ++i ) {
    const reco::Muon* mu = isPatMuon ? &(patMuonHandle->at(i)) : &(recoMuonHandle->at(i));
    const double aeta = std::abs(mu->eta());

    // Collect matching information
    std::set<int> matchedStLoose, matchedStTight;
    int nLayerLoose = 0, nLayerTight = 0;
    int nLayerFirstStLoose = 0, nLayerFirstStTight = 0;
    int nLayerLastStLoose = 0, nLayerLastStTight = 0;

    const auto& muMatches = mu->matches();
    for ( const auto& muMatch : muMatches ) {
      if ( muMatch.detector() != 3 ) continue;
      if ( muMatch.rpcMatches.empty() ) continue;

      const RPCDetId rpcDet(muMatch.id);
      const int region = rpcDet.region();
      const int st = rpcDet.station();
      const int ring = rpcDet.ring();
      const bool isFirstStation = (st > 1) ? false : ( (region != 0 and ring == 3) ? false : true );
      const bool isLastStation  = [&](){
        if ( st == 4 ) return true;
        if ( region == 0 ) {
          if ( std::abs(ring) == 2 and st == 3 and aeta >= 0.8 and aeta <= 0.9 ) return true; // RB3 wheel += 2
        }
        else if ( ring == 3 ) {
          if ( st == 1 and aeta >= 0.9 and aeta <= 1.0  ) return true; // RE1/3
          if ( st == 2 and aeta >= 1.0 and aeta <= 1.15 ) return true; // RE2/3
          if ( st == 3 and aeta >= 1.1 and aeta <= 1.2  ) return true; // RE3/3
        }
        return false;
      }();

      double dx = 1e9, dxErr = 1e9;
      for ( auto& rpcMatch : muMatch.rpcMatches ) {
        const double dx1 = std::abs(muMatch.x-rpcMatch.x);
        if ( dx1 < dx ) {
          dx = dx1;
          dxErr = muMatch.xErr;
        }
      }
      if ( dx < 20 or dx < 4*dxErr ) {
        matchedStLoose.insert(st);
        ++nLayerLoose;
        if ( isFirstStation ) ++nLayerFirstStLoose;
        if ( isLastStation ) ++nLayerLastStLoose;
      }
      if ( dx < 3 and dx < 4*dxErr ) {
        matchedStTight.insert(st);
        ++nLayerTight;
        if ( isFirstStation ) ++nLayerFirstStTight;
        if ( isLastStation ) ++nLayerLastStTight;
      }
    }

    // Calculate RPC ID based on the matching information
    idMaps["Loose"].push_back(nLayerLoose>=2);
    idMaps["Tight"].push_back(nLayerTight>=2);
    idMaps["TwoStationLoose"].push_back(matchedStLoose.size()>=2);
    idMaps["TwoStationTight"].push_back(matchedStTight.size()>=2);
    idMaps["LastStationLoose"].push_back(nLayerLastStLoose>0);
    idMaps["LastStationTight"].push_back(nLayerLastStTight>0);
    idMaps["SecondStationLoose"].push_back(nLayerLoose>=2 and (nLayerLoose-nLayerFirstStLoose)>0);
    idMaps["SecondStationTight"].push_back(nLayerTight>=2 and (nLayerTight-nLayerFirstStTight)>0);
  
  }
  
  std::auto_ptr<edm::ValueMap<float> > out;
  for ( auto key = idMaps.begin(); key != idMaps.end(); ++key ) {
    const auto& name = key->first;
    const auto& idValues = key->second;

    out.reset(new edm::ValueMap<float>());
    edm::ValueMap<float>::Filler filler(*out);
    if ( isPatMuon ) filler.insert(patMuonHandle, idValues.begin(), idValues.end());
    else filler.insert(recoMuonHandle, idValues.begin(), idValues.end());
    filler.fill();
    event.put(out, name);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCMuonIdMapProducer);

