#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "TChain.h"
#include "TProfile.h"

using namespace std;

typedef std::vector<double> doubles;
typedef std::vector<int> ints;

void muonRPCRecHitAnalysis()
{
  gROOT->ProcessLine(".x /users/jhgoh/work/rootlogon.C");

  TChain* chain = new TChain("RPC/tree");
  chain->Add("Run2012A.root");
  chain->Add("Run2012B.root");
  //TFile* f = TFile::Open("Run2012A.root");
  //TTree* chain = (TTree*)f->Get("RPC/chain");
  doubles* pt = new doubles;
  doubles* eta = new doubles;
  doubles* phi = new doubles;
  ints* nHit = new ints;
  ints* nTrackerHit = new ints;
  ints* nPixelHit = new ints;
  ints* nMuonHit = new ints;
  ints* nRPCHit = new ints;

  chain->SetBranchAddress("pt", &pt);
  chain->SetBranchAddress("eta", &eta);
  chain->SetBranchAddress("phi", &phi);
  chain->SetBranchAddress("nHit", &nHit);
  chain->SetBranchAddress("nTrackerHit", &nTrackerHit);
  chain->SetBranchAddress("nPixelHit", &nPixelHit);
  chain->SetBranchAddress("nMuonHit", &nMuonHit);
  chain->SetBranchAddress("nRPCHit", &nRPCHit);

  TH1F* hEta = new TH1F("hEta", "hEta;Pseudorapidity #eta;Events", 100, -2.5, 2.5);
  TH1F* hPhi = new TH1F("hPhi", "hPhi;Azimuthal angle #phi;Events", 100, -TMath::Pi(), TMath::Pi());

  TH2F* hEtaVsPhi = new TH2F("hEtaVsPhi", "hEtaVsPhi;Pseudorapidity #eta;Azimuthal angle #phi", 100, -2.5, 2.5, 100, -TMath::Pi(), TMath::Pi());
  TH2F* hEtaVsPhi1RPC = new TH2F("hEtaVsPhi1RPC", "hEtaVsPhi1RPC;Pseudorapidity #eta;Azimuthal angle #phi", 100, -2.5, 2.5, 100, -TMath::Pi(), TMath::Pi());
  TH2F* hEtaVsPhi2RPC = new TH2F("hEtaVsPhi2RPC", "hEtaVsPhi2RPC;Pseudorapidity #eta;Azimuthal angle #phi", 100, -2.5, 2.5, 100, -TMath::Pi(), TMath::Pi());
  TH2F* hEtaVsPhi3RPC = new TH2F("hEtaVsPhi3RPC", "hEtaVsPhi3RPC;Pseudorapidity #eta;Azimuthal angle #phi", 100, -2.5, 2.5, 100, -TMath::Pi(), TMath::Pi());

  TH2F* hEtaVsNHit = new TH2F("hEtaVsNHit", "hEtaVsNHit;Pseudorapidity #eta;Number of hits", 100, -2.5, 2.5, 100, 0, 100);
  TH2F* hEtaVsNMuonHit = new TH2F("hEtaVsNMuonHit", "hEtaVsNRPCHit;Pseudorapidity #eta;Number of Muon hits", 100, -2.5, 2.5, 70, 0, 70);
  TH2F* hEtaVsNRPCHit = new TH2F("hEtaVsNRPCHit", "hEtaVsNRPCHit;Pseudorapidity #eta;Number of RPC hits", 100, -2.5, 2.5, 10, 0, 10);
  TH2F* hEtaVsFRPCHit = new TH2F("hEtaVsFRPCHit", "hEtaVsFRPCHit;Pseudorapidity #eta;nRPCHit/nMuonHit", 100, -2.5, 2.5, 25, 0, 1);

  TProfile* pEtaVsNHit = new TProfile("pEtaVsNHit", "pEtaVsNHit;Pseudorapidity #eta;Number of hits", 100, -2.5, 2.5, 0, 100);
  TProfile* pEtaVsNMuonHit = new TProfile("pEtaVsNMuonHit", "pEtaVsNRPCHit;Pseudorapidity #eta;Number of Muon hits", 100, -2.5, 2.5, 0, 70);
  TProfile* pEtaVsNRPCHit = new TProfile("pEtaVsNRPCHit", "pEtaVsNRPCHit;Pseudorapidity #eta;Number of RPC hits", 100, -2.5, 2.5, 0, 10);
  TProfile* pEtaVsFRPCHit = new TProfile("pEtaVsFRPCHit", "pEtaVsFRPCHit;Pseudorapidity #eta;nRPCHit/nMuonHit", 100, -2.5, 2.5, 0, 1);

  TH2F* hPhiVsNRPCHit = new TH2F("hPhiVsNRPCHit", "hPhiVsNRPCHit;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 10, 0, 10);
  TH2F* hPhiVsNRPCHitA = new TH2F("hPhiVsNRPCHitA", "hPhiVsNRPCHitA;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 10, 0, 10);
  TH2F* hPhiVsNRPCHitB = new TH2F("hPhiVsNRPCHitB", "hPhiVsNRPCHitB;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 10, 0, 10);
  TH2F* hPhiVsNRPCHitC = new TH2F("hPhiVsNRPCHitC", "hPhiVsNRPCHitC;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 10, 0, 10);
  TH2F* hPhiVsNRPCHitD = new TH2F("hPhiVsNRPCHitD", "hPhiVsNRPCHitD;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 10, 0, 10);

  TProfile* pPhiVsNRPCHit = new TProfile("pPhiVsNRPCHit", "pPhiVsNRPCHit;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 10);
  TProfile* pPhiVsNRPCHitA = new TProfile("pPhiVsNRPCHitA", "pPhiVsNRPCHitA;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 10);
  TProfile* pPhiVsNRPCHitB = new TProfile("pPhiVsNRPCHitB", "pPhiVsNRPCHitB;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 10);
  TProfile* pPhiVsNRPCHitC = new TProfile("pPhiVsNRPCHitC", "pPhiVsNRPCHitC;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 10);
  TProfile* pPhiVsNRPCHitD = new TProfile("pPhiVsNRPCHitD", "pPhiVsNRPCHitD;Azimuthal angle #phi;Number of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 10);

  TH2F* hPhiVsFRPCHit  = new TH2F("hPhiVsFRPCHit" , "hPhiVsFRPCHit;Azimuthal angle #phi;Fraction of RPC hits" , 100, -TMath::Pi(), TMath::Pi(), 25, 0, 1);
  TH2F* hPhiVsFRPCHitA = new TH2F("hPhiVsFRPCHitA", "hPhiVsFRPCHitA;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 25, 0, 1);
  TH2F* hPhiVsFRPCHitB = new TH2F("hPhiVsFRPCHitB", "hPhiVsFRPCHitB;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 25, 0, 1);
  TH2F* hPhiVsFRPCHitC = new TH2F("hPhiVsFRPCHitC", "hPhiVsFRPCHitC;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 25, 0, 1);
  TH2F* hPhiVsFRPCHitD = new TH2F("hPhiVsFRPCHitD", "hPhiVsFRPCHitD;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 25, 0, 1);

  TProfile* pPhiVsFRPCHit  = new TProfile("pPhiVsFRPCHit" , "pPhiVsFRPCHit;Azimuthal angle #phi;Fraction of RPC hits" , 100, -TMath::Pi(), TMath::Pi(), 0, 1);
  TProfile* pPhiVsFRPCHitA = new TProfile("pPhiVsFRPCHitA", "pPhiVsFRPCHitA;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 1);
  TProfile* pPhiVsFRPCHitB = new TProfile("pPhiVsFRPCHitB", "pPhiVsFRPCHitB;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 1);
  TProfile* pPhiVsFRPCHitC = new TProfile("pPhiVsFRPCHitC", "pPhiVsFRPCHitC;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 1);
  TProfile* pPhiVsFRPCHitD = new TProfile("pPhiVsFRPCHitD", "pPhiVsFRPCHitD;Azimuthal angle #phi;Fraction of RPC hits", 100, -TMath::Pi(), TMath::Pi(), 0, 1);

  for ( int entry=0, nEntry=chain->GetEntries(); entry<nEntry; ++entry )
  {
    chain->GetEntry(entry);
    const int nMuon = pt->size();
    for ( int i=0; i<nMuon; ++i )
    {
      const double muonEta = eta->at(i);
      const double muonPhi = phi->at(i);
      const double muonNRPCHit = nRPCHit->at(i);
      const double muonFRPCHit = muonNRPCHit/nMuonHit->at(i);

      hEta->Fill(muonEta);
      hPhi->Fill(muonPhi);
      hEtaVsPhi->Fill(muonEta, muonPhi);
      if ( muonNRPCHit > 0 ) hEtaVsPhi1RPC->Fill(muonEta, muonPhi);
      if ( muonNRPCHit > 1 ) hEtaVsPhi2RPC->Fill(muonEta, muonPhi);
      if ( muonNRPCHit > 2 ) hEtaVsPhi3RPC->Fill(muonEta, muonPhi);
      
      hEtaVsPhi->Fill(muonEta, muonPhi);
      hEtaVsPhi->Fill(muonEta, muonPhi);

      hEtaVsNHit->Fill(muonEta, nHit->at(i));
      hEtaVsNMuonHit->Fill(muonEta, nMuonHit->at(i));
      hEtaVsNRPCHit->Fill(muonEta, muonNRPCHit);
      hEtaVsFRPCHit->Fill(muonEta, muonFRPCHit);

      hPhiVsNRPCHit->Fill(muonPhi, muonNRPCHit);
      hPhiVsFRPCHit->Fill(muonPhi, muonFRPCHit);
      pPhiVsNRPCHit->Fill(muonPhi, muonNRPCHit);
      pPhiVsFRPCHit->Fill(muonPhi, muonFRPCHit);
      if ( fabs(muonEta) < 0.8 )
      {
        hPhiVsNRPCHitA->Fill(muonPhi, muonNRPCHit);
        hPhiVsFRPCHitA->Fill(muonPhi, muonFRPCHit);
        pPhiVsNRPCHitA->Fill(muonPhi, muonNRPCHit);
        pPhiVsFRPCHitA->Fill(muonPhi, muonFRPCHit);
      }
      else if ( fabs(muonEta) < 1.2 )
      {
        hPhiVsNRPCHitB->Fill(muonPhi, muonNRPCHit);
        hPhiVsFRPCHitB->Fill(muonPhi, muonFRPCHit);
        pPhiVsNRPCHitB->Fill(muonPhi, muonNRPCHit);
        pPhiVsFRPCHitB->Fill(muonPhi, muonFRPCHit);
      }
      else if ( fabs(muonEta) < 1.6 )
      {
        hPhiVsNRPCHitC->Fill(muonPhi, muonNRPCHit);
        hPhiVsFRPCHitC->Fill(muonPhi, muonFRPCHit);
        pPhiVsNRPCHitC->Fill(muonPhi, muonNRPCHit);
        pPhiVsFRPCHitC->Fill(muonPhi, muonFRPCHit);
      }
      else if ( fabs(muonEta) < 2.5 )
      {
        hPhiVsNRPCHitD->Fill(muonPhi, muonNRPCHit);
        hPhiVsFRPCHitD->Fill(muonPhi, muonFRPCHit);
        pPhiVsNRPCHitD->Fill(muonPhi, muonNRPCHit);
        pPhiVsFRPCHitD->Fill(muonPhi, muonFRPCHit);
      }

    }
  }
    
  TCanvas* cEta = new TCanvas("cEta", "cEta", 500, 500);
  hEta->Draw();

  TCanvas* cPhi = new TCanvas("cPhi", "cPhi", 500, 500);
  hPhi->Draw();

  TCanvas* cEtaVsNHit = new TCanvas("cEtaVsNHit", "cEtaVsNHit", 500, 500);
  hEtaVsNHit->Draw("COLZ");
  pEtaVsNHit->Draw("same");

  TCanvas* cEtaVsNMuonHit = new TCanvas("cEtaVsNMuonHit", "cEtaVsNMuonHit", 500, 500);
  hEtaVsNMuonHit->Draw("COLZ");
  pEtaVsNMuonHit->Draw("same");

  TCanvas* cEtaVsNRPCHit = new TCanvas("cEtaVsNRPCHit", "cEtaVsNRPCHit", 500, 500);
  hEtaVsNRPCHit->Draw("COLZ");
  pEtaVsNRPCHit->Draw("same");

  TCanvas* cEtaVsFRPCHit = new TCanvas("cEtaVsFRPCHit", "cEtaVsFRPCHit", 500, 500);
  hEtaVsFRPCHit->Draw("COLZ");
  pEtaVsFRPCHit->Draw("same");

  TCanvas* cEtaVsPhi = new TCanvas("cEtaVsPhi", "cEtaVsPhi", 500, 500);
  hEtaVsPhi->Draw("COLZ");

  TCanvas* cEtaVsPhi1RPC = new TCanvas("cEtaVsPhi1RPC", "cEtaVsPhi1RPC", 500, 500);
  hEtaVsPhi1RPC->Draw("COLZ");

  TCanvas* cEtaVsPhi2RPC = new TCanvas("cEtaVsPhi2RPC", "cEtaVsPhi2RPC", 500, 500);
  hEtaVsPhi2RPC->Draw("COLZ");

  TCanvas* cEtaVsPhi3RPC = new TCanvas("cEtaVsPhi3RPC", "cEtaVsPhi3RPC", 500, 500);
  hEtaVsPhi3RPC->Draw("COLZ");

  TCanvas* cPhiVsNRPCHit = new TCanvas("cPhiVsNRPCHit", "cPhiVsNRPCHit", 500, 500);
  hPhiVsNRPCHit->Draw("COLZ");
  pPhiVsNRPCHit->Draw("same");

  TCanvas* cPhiVsNRPCHitA = new TCanvas("cPhiVsNRPCHitA", "cPhiVsNRPCHitA", 500, 500);
  hPhiVsNRPCHitA->Draw("COLZ");
  pPhiVsNRPCHitA->Draw("same");

  TCanvas* cPhiVsNRPCHitB = new TCanvas("cPhiVsNRPCHitB", "cPhiVsNRPCHitB", 500, 500);
  hPhiVsNRPCHitB->Draw("COLZ");
  pPhiVsNRPCHitB->Draw("same");

  TCanvas* cPhiVsNRPCHitC = new TCanvas("cPhiVsNRPCHitC", "cPhiVsNRPCHitC", 500, 500);
  hPhiVsNRPCHitC->Draw("COLZ");
  pPhiVsNRPCHitC->Draw("same");

  TCanvas* cPhiVsNRPCHitD = new TCanvas("cPhiVsNRPCHitD", "cPhiVsNRPCHitD", 500, 500);
  hPhiVsNRPCHitD->Draw("COLZ");
  pPhiVsNRPCHitD->Draw("same");

  TCanvas* cPhiVsFRPCHit = new TCanvas("cPhiVsFRPCHit", "cPhiVsFRPCHit", 500, 500);
  hPhiVsFRPCHit->Draw("COLZ");
  pPhiVsFRPCHit->Draw("same");

  TCanvas* cPhiVsFRPCHitA = new TCanvas("cPhiVsFRPCHitA", "cPhiVsFRPCHitA", 500, 500);
  hPhiVsFRPCHitA->Draw("COLZ");
  pPhiVsFRPCHitA->Draw("same");

  TCanvas* cPhiVsFRPCHitB = new TCanvas("cPhiVsFRPCHitB", "cPhiVsFRPCHitB", 500, 500);
  hPhiVsFRPCHitB->Draw("COLZ");
  pPhiVsFRPCHitB->Draw("same");

  TCanvas* cPhiVsFRPCHitC = new TCanvas("cPhiVsFRPCHitC", "cPhiVsFRPCHitC", 500, 500);
  hPhiVsFRPCHitC->Draw("COLZ");
  pPhiVsFRPCHitC->Draw("same");

  TCanvas* cPhiVsFRPCHitD = new TCanvas("cPhiVsFRPCHitD", "cPhiVsFRPCHitD", 500, 500);
  hPhiVsFRPCHitD->Draw("COLZ");
  pPhiVsFRPCHitD->Draw("same");


}
