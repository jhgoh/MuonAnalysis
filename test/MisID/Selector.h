//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 26 23:38:43 2017 by ROOT version 6.02/13
// from TTree tree/tree
// found on file: 20170226_1/JetHT_2015D/ntuple_1.root
//////////////////////////////////////////////////////////

#ifndef Selector_h
#define Selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1D.h>

#include <vector>
#include <map>
#include <string>

class Selector : public TSelector {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  UChar_t         run;
  UChar_t         lumi;
  ULong64_t       event;
  UChar_t         nPV;
  UChar_t         nSV;
  UChar_t         nGen;
  Float_t         vtx_mass;
  Float_t         vtx_pt;
  Float_t         vtx_lxy;
  Float_t         vtx_vz;
  Float_t         vtx_ndof;
  Float_t         vtx_chi2;
  vector<int>     *trk_pdgId;
  vector<float>   *trk_pt;
  vector<float>   *trk_eta;
  vector<float>   *trk_phi;
  vector<int>     *mu_q;
  vector<float>   *mu_dR;
  vector<float>   *mu_pt;
  Float_t         gen_dR;

  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_lumi;   //!
  TBranch        *b_event;   //!
  TBranch        *b_nPV;   //!
  TBranch        *b_nSV;   //!
  TBranch        *b_nGen;   //!
  TBranch        *b_vtx_mass;   //!
  TBranch        *b_vtx_pt;   //!
  TBranch        *b_vtx_lxy;   //!
  TBranch        *b_vtx_vz;   //!
  TBranch        *b_vtx_ndof;   //!
  TBranch        *b_vtx_chi2;   //!
  TBranch        *b_trk_pdgId;   //!
  TBranch        *b_trk_pt;   //!
  TBranch        *b_trk_eta;   //!
  TBranch        *b_trk_phi;   //!
  TBranch        *b_mu_q;   //!
  TBranch        *b_mu_dR;   //!
  TBranch        *b_mu_pt;   //!
  TBranch        *b_gen_dR;   //!

  std::map<std::string, std::vector<bool>*> mu_idVars_;
  std::map<std::string, TBranch*> b_mu_idVars_;

  Selector(TTree * /*tree*/ =0) : fChain(0) { }
  virtual ~Selector() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  TDirectory* FetchDir(TDirectory* dir, std::string path, const bool ignoreLast=true) const;
  std::string BaseName(const std::string path) const;

  std::string mode_;
  std::map<std::string, TH1D*> hists_;

  const static std::vector<std::string> idNames_;
  double minLxy_, maxLxy_;

  ClassDef(Selector,0);
};

#endif

#ifdef Selector_cxx
const std::vector<std::string> Selector::idNames_ = {
  "Tight", "Medium", "Loose", "Soft", "HighPt",
  "GLB", "TRK", "STA", "RPC",
};

void Selector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  trk_pdgId = 0;
  trk_pt = 0;
  trk_eta = 0;
  trk_phi = 0;
  mu_q = 0;
  mu_dR = 0;
  mu_pt = 0;

  for ( auto idName : idNames_ ) {
    mu_idVars_[idName] = 0;
  }
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
  fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
  fChain->SetBranchAddress("nGen", &nGen, &b_nGen);
  fChain->SetBranchAddress("vtx_mass", &vtx_mass, &b_vtx_mass);
  fChain->SetBranchAddress("vtx_pt", &vtx_pt, &b_vtx_pt);
  fChain->SetBranchAddress("vtx_lxy", &vtx_lxy, &b_vtx_lxy);
  fChain->SetBranchAddress("vtx_vz", &vtx_vz, &b_vtx_vz);
  fChain->SetBranchAddress("vtx_ndof", &vtx_ndof, &b_vtx_ndof);
  fChain->SetBranchAddress("vtx_chi2", &vtx_chi2, &b_vtx_chi2);
  fChain->SetBranchAddress("trk_pdgId", &trk_pdgId, &b_trk_pdgId);
  fChain->SetBranchAddress("trk_pt", &trk_pt, &b_trk_pt);
  fChain->SetBranchAddress("trk_eta", &trk_eta, &b_trk_eta);
  fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
  fChain->SetBranchAddress("mu_q", &mu_q, &b_mu_q);
  fChain->SetBranchAddress("mu_dR", &mu_dR, &b_mu_dR);
  fChain->SetBranchAddress("mu_pt", &mu_pt, &b_mu_pt);
  for ( auto idName : idNames_ ) {
    fChain->SetBranchAddress(Form("mu_is%s", idName.c_str()), &mu_idVars_[idName], &b_mu_idVars_[idName]);
  }
  fChain->SetBranchAddress("gen_dR", &gen_dR, &b_gen_dR);
}

Bool_t Selector::Notify()
{
  return kTRUE;
}

TDirectory* Selector::FetchDir(TDirectory* dir, std::string path, const bool ignoreLast) const
{
  //if ( ignoreLast ) path = path.substr(0, path.rfind('/'));

  std::string::size_type prev_pos = 0, pos = 0;
  while((pos = path.find('/', pos)) != std::string::npos) {
    auto substring = path.substr(prev_pos, pos-prev_pos);
    prev_pos = ++pos;

    TDirectory* subdir = dir->GetDirectory(substring.c_str());
    if ( subdir ) dir = subdir;
    else {
      dir = dir->mkdir(substring.c_str());
    }
  }

  return dir;
}

std::string Selector::BaseName(const std::string path) const
{
  return path.substr(path.rfind('/')+1, std::string::npos);
}

#endif // #ifdef Selector_cxx
