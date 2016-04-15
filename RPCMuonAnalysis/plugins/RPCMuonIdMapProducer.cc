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
    {"LastStationLoose", vfloat()},
    {"LastStationTight", vfloat()},
    {"SecondStationLoose", vfloat()},
    {"SecondStationTight", vfloat()},
  };

  for ( int i=0, n=muHandle->size(); i<n; ++i ) {
    edm::Ref<MuonColl> muRef(muHandle, i);

    // Collect matching information
    std::set<int> matchedLoose, matchedTight;

    const auto& muMatches = muRef->matches();
    for ( const auto& muMatch : muMatches ) {
      if ( muMatch.detector() != 3 ) continue;
      if ( muMatch.rpcMatches.empty() ) continue;

      const int st = muMatch.station();

      const double dx = muMatch.dist();
      const double dxErr = muMatch.distErr();
      if ( dx < 20 or dx < 4*dxErr ) matchedLoose.insert(st);
      if ( dx < 3 and dx < 4*dxErr ) matchedTight.insert(st);
    }

    // Calculate RPC ID based on the matching information
    const int nMatchLoose = matchedLoose.size();
    const int nMatchTight = matchedTight.size();
    int nMatchLooseSSt = nMatchLoose, nMatchTightSSt = nMatchTight;
    const int minSt = RPCDetId::minStationId;
    const int maxSt = RPCDetId::maxStationId;
    if ( matchedLoose.count( minSt) != 0 ) --nMatchLooseSSt;
    if ( matchedLoose.count(-minSt) != 0 ) --nMatchLooseSSt;
    if ( matchedTight.count( minSt) != 0 ) --nMatchTightSSt;
    if ( matchedTight.count(-minSt) != 0 ) --nMatchTightSSt;

    idMaps["Loose"].push_back(nMatchLoose>=2);
    idMaps["Tight"].push_back(nMatchTight>=2);
    idMaps["LastStationLoose"].push_back(matchedLoose.count(maxSt) != 0);
    idMaps["LastStationTight"].push_back(matchedTight.count(maxSt) != 0);
    idMaps["SecondStationLoose"].push_back(nMatchLooseSSt>=2);
    idMaps["SecondStationTight"].push_back(nMatchTightSSt>=2);
  
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

