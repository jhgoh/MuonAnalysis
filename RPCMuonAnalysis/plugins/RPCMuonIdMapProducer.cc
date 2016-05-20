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
  typedef std::vector<pat::Muon> MuonColl;
  edm::EDGetTokenT<MuonColl> muToken_;

  typedef std::vector<float> vfloat;
};

RPCMuonIdMapProducer::RPCMuonIdMapProducer(const edm::ParameterSet& pset):
  muToken_(consumes<MuonColl>(pset.getParameter<edm::InputTag>("src")))
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
  edm::Handle<MuonColl> muHandle;
  event.getByToken(muToken_, muHandle);

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

  for ( int i=0, n=muHandle->size(); i<n; ++i ) {
    edm::Ref<MuonColl> muRef(muHandle, i);
    const double aeta = std::abs(muRef->eta());

    // Collect matching information
    std::set<int> matchedStLoose, matchedStTight;
    int nLayerLoose = 0, nLayerTight = 0;
    int nLayerFirstStLoose = 0, nLayerFirstStTight = 0;
    int nLayerLastStLoose = 0, nLayerLastStTight = 0;

    const auto& muMatches = muRef->matches();
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
        if ( region == 1 ) {
          if ( std::abs(ring) == 2 and st == 3 and aeta >= 0.8 and aeta <= 0.9 ) return true; // RB3 wheel += 2
        }
        else if ( ring == 3 ) {
          if ( st == 1 and aeta >= 0.9 and aeta <= 1.0  ) return true; // RE1/3
          if ( st == 2 and aeta >= 1.0 and aeta <= 1.15 ) return true; // RE2/3
          if ( st == 3 and aeta >= 1.1 and aeta <= 1.2  ) return true; // RE3/3
        }
        return false;
      }();

      const double dx = muMatch.dist();
      const double dxErr = muMatch.distErr();
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
    filler.insert(muHandle, idValues.begin(), idValues.end());
    filler.fill();
    event.put(out, name);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCMuonIdMapProducer);

