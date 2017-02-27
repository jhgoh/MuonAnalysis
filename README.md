## Installation
Install the muon analysis packages
```
source /cvmfs/cms.cern.ch/crab3/crab.sh

cmsrel CMSSW_7_6_6_patch1
cd CMSSW_7_6_6_patch1/src
cmsenv
git-cms-init

git-cms-addpkg MuonAnalysis
cd MuonAnalysis
git clone https://github.com/cms-analysis/MuonAnalysis-TagAndProbe TagAndProbe
git clone https://github.com/jhgoh/MuonAnalysis-MuonIdentification MuonIdentification

cd ..
scram b -j8
```

## Muon misID measurement
This package includes modules and scripts to produce flat ntuples to measure the muon misidentification probability.
The workflow starts from ntuple production, histogramming and fitting.

Testing with small sample (you may have to modify the configuration file to read valid AOD root files.
```
cd $CMSSW_BASE/src/MuonAnalysis/MuonIdentification/test/MisID
cmsRun prod_RD_cfg.py
#cmsRun prod_MC_cfg.py
```

You can submit crab jobs to process full dataset. Currently the configuration file is set to read 2015 data and MC.

```
crab submit crabConfigRD.py
crab submit crabConfigMC.py
```

You can find out TTbar powheg sample and JetHT 2015D ntuple in the eos

```
eos ls /store/user/jhgoh/MuonMisID/20150401_2/JetHT_2015D/JetHT/crab_20160401_143746/160401_123805/0000
eos ls /store/user/jhgoh/MuonMisID/20150401_2/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160401_143809/160401_123825/0000
```

new root files with histograms will be created using run_selector.C
A root file containing invariant mass distributions for each pt, eta bins will be created. This use multiple CPUs using the Proof.

```
root -b -q -l run_selector.C
```

Run the fitter to extract misID probability.

```
python fit.py ks hist_ks.root fit_ks.root
python fit.py phi hist_phi.root fit_phi.root
python fit.py lamb hist_lamb.root fit_lamb.root
```

## RPCMuon efficiency measurement using the Tag and Probe method

Test the ntuple production

```
cd MuonAnalysis/MuonIdentification/test/TagAndProbe
cmsRun tp_prod_MC_cfg.py
root -l tnp.root
```

Submit crab jobs

```
crab submit submitCrabRD.py
crab submit submitCrabMC.py
```

DYJets MG5 and Run2015D singleMuon datasets are under production. Saving them into the eos.

```
/store/user/jhgoh/RPCMuonEff/20150404_1/DYJets_MG/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_20160404_173038/160404_153124/0000:
/store/user/jhgoh/RPCMuonEff/20150404_1/SingleMuon_2015D/SingleMuon/crab_20160404_173352/160404_153427/0000
```
It will take some time to be finished, but you can start the fitting before finishing it.
configuration files to do the tnp fitting is already prepared.

```
cmsRun tp_fit_MC_cfg.py
cmsRun tp_fit_RD_cfg.py
```
