#include <iostream>

void run_selector()
{
  gROOT->SetBatch();
  TProof::Open("");
  
  std::vector<const char*> modes = {"ks", "phi", "lamb"};
  for ( auto mode : modes ) {
    TChain chain(Form("%s/tree", mode));
    chain.Add("TT_powheg/ntuple_*.root");
    //chain.Add("JetHT_2015D/ntuple_*.root");
    chain.SetProof();

    chain.Process("Selector.C+", mode);
  }

}

