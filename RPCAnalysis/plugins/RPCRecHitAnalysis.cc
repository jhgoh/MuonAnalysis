#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class RPCRecHitAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  RPCRecHitAnalysis(const edm::ParameterSet& pset);
  virtual ~RPCRecHitAnalysis() {};
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  edm::EDGetTokenT<RPCRecHitCollection> recHitToken_;

  TH1D* hBx_;
  TH1D* hNHit_, * hCLs_, * hAvgCLs_;
  TH1D* hBx0_NHit_, * hBx0_CLs_, * hBx0_AvgCLs_;
  TH1D* hBxX_NHit_, * hBxX_CLs_, * hBxX_AvgCLs_;
};

RPCRecHitAnalysis::RPCRecHitAnalysis(const edm::ParameterSet& pset):
  recHitToken_(consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("recHit")))
{
  usesResource("TFileService");

  edm::Service<TFileService> fs;
  hBx_     = fs->make<TH1D>("hBx", "Bunch crossing", 21, -10.5, 10.5);
  hNHit_   = fs->make<TH1D>("hNHit", "Number of hits", 100, 0, 100);
  hCLs_    = fs->make<TH1D>("hCLs", "Cluster size", 100, 0, 100);
  hAvgCLs_ = fs->make<TH1D>("hAvgCLs", "Average cluster size", 100, 0, 5);
  auto dirBx0 = fs->mkdir("Bx0");
  hBx0_NHit_   = dirBx0.make<TH1D>("hNHit", "Number of hits", 100, 0, 100);
  hBx0_CLs_    = dirBx0.make<TH1D>("hCLs", "Cluster size", 100, 0, 100);
  hBx0_AvgCLs_ = dirBx0.make<TH1D>("hAvgCLs", "Average cluster size", 100, 0, 5);
  auto dirBxX = fs->mkdir("BxX");
  hBxX_NHit_   = dirBxX.make<TH1D>("hNHit", "Number of hits", 100, 0, 100);
  hBxX_CLs_    = dirBxX.make<TH1D>("hCLs", "Cluster size", 100, 0, 100);
  hBxX_AvgCLs_ = dirBxX.make<TH1D>("hAvgCLs", "Average cluster size", 100, 0, 5);
}

void RPCRecHitAnalysis::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<RPCRecHitCollection> recHitHandle;
  event.getByToken(recHitToken_, recHitHandle);

  const int nHitAll = recHitHandle->size();
  int nHitBx0 = 0;
  double sumCLs = 0, sumCLsBx0 = 0, sumCLsBxX = 0;
  for ( const auto& recHit : *recHitHandle ) {
    const int cls = recHit.clusterSize();
    const int bx = recHit.BunchX();

    sumCLs += cls;
    hCLs_->Fill(cls);
    hBx_->Fill(bx);

    if ( bx == 0 ) {
      ++nHitBx0;
      sumCLsBx0 += cls;
      hBx0_CLs_->Fill(cls);
    }
    else {
      sumCLsBxX += cls;
      hBxX_CLs_->Fill(cls);
    }
  }
  hNHit_->Fill(nHitAll);
  hAvgCLs_->Fill(sumCLs/nHitAll);

  hBx0_NHit_->Fill(nHitBx0);
  hBx0_AvgCLs_->Fill(sumCLsBx0/nHitBx0);

  const int nHitBxX = nHitAll-nHitBx0;
  hBxX_NHit_->Fill(nHitBxX);
  hBxX_AvgCLs_->Fill(nHitBxX == 0 ? -1 : sumCLsBxX/nHitBxX);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RPCRecHitAnalysis);

