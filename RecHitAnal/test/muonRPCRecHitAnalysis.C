#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>

using namespace std;

typedef std::vector<double> doubles;
typedef std::vector<int> ints;

void analyzeRPCTree()
{
  TFile* f = TFile::Open("RPCAnalysis.root");
  TTree* tree = (TTree*)f->Get("RPC/tree");
  doubles* pt = new doubles;
  doubles* eta = new doubles;
  doubles* phi = new doubles;
  ints* nHit = new ints;
  ints* nTrackerHit = new ints;
  ints* nPixelHit = new ints;
  ints* nMuonHit = new ints;
  ints* nRPCHit = new ints;

  tree->SetBranchAddress("pt", &pt);
  tree->SetBranchAddress("eta", &eta);
  tree->SetBranchAddress("phi", &phi);
  tree->SetBranchAddress("nHit", &nHit);
  tree->SetBranchAddress("nTrackerHit", &nTrackerHit);
  tree->SetBranchAddress("nPixelHit", &nPixelHit);
  tree->SetBranchAddress("nMuonHit", &nMuonHit);
  tree->SetBranchAddress("nRPCHit", &nRPCHit);

  TH1F* hEta = new TH1F("hEta", "hEta;Pseudorapidity #eta;Events", 100, -2.5, 2.5);
  TH2F* hEtaVsNHit = new TH2F("hEtaVsNHit", "hEtaVsNHit;Pseudorapidity #eta;Number of hits", 100, -2.5, 2.5, 100, 0, 100);
  TH2F* hEtaVsNRPCHit = new TH2F("hEtaVsNRPCHit", "hEtaVsNRPCHit;Pseudorapidity #eta;Number of RPC hits", 100, -2.5, 2.5, 100, 0, 100);

  for ( int entry=0, nEntry=tree->GetEntries(); entry<nEntry; ++entry )
  {
    tree->GetEntry(entry);
    const int nMuon = pt->size();
    for ( int i=0; i<nMuon; ++i )
    {
      hEta->Fill(eta->at(i));
      hEtaVsNHit->Fill(eta->at(i), nHit->at(i));
      hEtaVsNRPCHit->Fill(eta->at(i), nRPCHit->at(i));
    }
  }
    
  TCanvas* cEta = new TCanvas("cEta", "cEta", 500, 500);
  hEta->Draw();

  TCanvas* cEtaVsNHit = new TCanvas("cEtaVsNHit", "cEtaVsNHit", 500, 500);
  hEtaVsNHit->Draw("COLZ");

  TCanvas* cEtaVsNRPCHit = new TCanvas("cEtaVsNRPCHit", "cEtaVsNRPCHit", 500, 500);
  hEtaVsNRPCHit->Draw("COLZ");
}
