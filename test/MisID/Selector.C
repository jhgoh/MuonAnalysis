#define Selector_cxx
// The class definition in Selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Selector.C")
// root> T->Process("Selector.C","some options")
// root> T->Process("Selector.C+")
//

#include "Selector.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

void Selector::Begin(TTree * /*tree*/)
{
  TString option = GetOption();
  mode_ = option.Data();

}

void Selector::SlaveBegin(TTree * /*tree*/)
{
  TString option = GetOption();
  mode_ = option.Data();

  std::vector<double> ptbins = {4, 5, 6, 8, 10, 15, 20, 30, 200};
  std::vector<double> etabins = {0, 0.9, 1.2, 1.6, 2.4};
  const int nbins = 50;
  double xmin = 0, xmax = 1;
  if      ( mode_ == "ks"   ) { xmin = 0.45; xmax = 0.55; }
  else if ( mode_ == "phi"  ) { xmin = 1.00; xmax = 1.04; }
  else if ( mode_ == "lamb" ) { xmin = 1.10; xmax = 1.13; }
  const double xbinW = 1000*(xmax-xmin);

  std::string xyTitle = Form("Mass (GeV);Candidates per %f MeV", xbinW);

  for ( auto idName : idNames_ ) {
    for ( int leg=1; leg<=2; ++leg ) {
      const string path  = Form("%s/%s_leg%d", mode_.c_str(), idName.c_str(), leg);

      TH1D* hFramePt = new TH1D(Form("%s/pt/hFrame", path.c_str()), "pt;p_{T} (GeV)", ptbins.size()-1, &ptbins[0]);
      hists_[hFramePt->GetName()] = hFramePt;
      for ( int i=1, n=ptbins.size(); i<=n; ++i ) {
        TH1D* hPass = new TH1D(Form("%s/pt/bin%d/hPass", path.c_str(), i),
                               Form("Passing;Mass (GeV);Candidates per %.3f MeV", xbinW), nbins, xmin, xmax);
        TH1D* hFail = new TH1D(Form("%s/pt/bin%d/hFail", path.c_str(), i),
                               Form("Failing;Mass (GeV);Candidates per %.3f MeV", xbinW), nbins, xmin, xmax);
        hists_[hPass->GetName()] = hPass;
        hists_[hFail->GetName()] = hFail;
      }

      TH1D* hFrameEta = new TH1D(Form("%s/abseta/hFrame", path.c_str()), "abseta;|#eta|", etabins.size()-1, &etabins[0]);
      hists_[hFrameEta->GetName()] = hFrameEta;
      for ( int i=1, n=etabins.size(); i<=n; ++i ) {
        TH1D* hPass = new TH1D(Form("%s/abseta/bin%d/hPass", path.c_str(), i),
                               Form("Passing;Mass (GeV);Candidates per %.3f MeV", xbinW), nbins, xmin, xmax);
        TH1D* hFail = new TH1D(Form("%s/abseta/bin%d/hFail", path.c_str(), i),
                               Form("Failing;Mass (GeV);Candidates per %.3f MeV", xbinW), nbins, xmin, xmax);
        hists_[hPass->GetName()] = hPass;
        hists_[hFail->GetName()] = hFail;
      }
    }
  }
}

Bool_t Selector::Process(Long64_t entry)
{
  GetEntry(entry);

  if ( vtx_lxy > 4.0 ) return false;

  for ( auto idName : idNames_ ) {
    for ( int leg=1; leg<=2; ++leg ) {
      const double pt = trk_pt->at(leg-1);
      const double eta = std::abs(trk_eta->at(leg-1));

      if ( pt < 4 or eta > 2.5 ) continue;
      if ( idName == "RPC" and eta > 2.1 ) continue;

      const char* path = Form("%s/%s_leg%d", mode_.c_str(), idName.c_str(), leg);
      auto hFramePt = hists_[Form("%s/pt/hFrame", path)];
      auto hFrameEta = hists_[Form("%s/abseta/hFrame", path)];

      const int binPt = hFramePt->GetXaxis()->FindBin(pt);
      const int binEta = hFrameEta->GetXaxis()->FindBin(eta);

      const bool isPass = mu_idVars_[idName]->at(leg-1) > 0.5 and mu_dR->at(leg-1) < 0.01;

      if ( binPt >= 1 and binPt <= hFramePt->GetNbinsX() ) {
        auto hPass = hists_[Form("%s/pt/bin%d/hPass", path, binPt)];
        auto hFail = hists_[Form("%s/pt/bin%d/hFail", path, binPt)];

        if ( isPass ) hPass->Fill(vtx_mass);
        else hFail->Fill(vtx_mass);
      }
      if ( binEta >= 1 and binEta <= hFrameEta->GetNbinsX() ) {
        auto hPass = hists_[Form("%s/abseta/bin%d/hPass", path, binEta)];
        auto hFail = hists_[Form("%s/abseta/bin%d/hFail", path, binEta)];

        if ( isPass ) hPass->Fill(vtx_mass);
        else hFail->Fill(vtx_mass);
      }
    }
  }

  return kTRUE;
}

void Selector::SlaveTerminate()
{
  for ( auto key=hists_.begin(); key!=hists_.end(); ++key ) {
    auto h = key->second;
    fOutput->Add(h);
  }
}

void Selector::Terminate()
{
  TFile fout(Form("hist_%s.root", mode_.c_str()), "recreate");
  for ( int i=0; i<fOutput->GetEntries(); ++i ) {
    TObject* obj = fOutput->At(i);
    TH1* h = dynamic_cast<TH1*>(obj);
    if ( h ) {
      auto dir = FetchDir(&fout, h->GetName());
      dir->cd();
      h->SetName(BaseName(h->GetName()).c_str());
      h->Write();
    }
  }
}
