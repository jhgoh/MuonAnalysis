## Installation
Install the muon analysis packages
```
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
source /cvmfs/cms.cern.ch/crab3/crab.sh

crab submit crabConfigRD.py
crab submit crabConfigMC.py
```

You can find out TTbar powheg sample and JetHT 2015D ntuple in the eos,
`/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15DR76-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/AODSIM` and `/JetHT/Run2015D-16Dec2015-v1/AOD`.

```
eos ls /store/user/jhgoh/MuonMisID/20170227_1/JetHT_Run2015D/ntuple_*.root
eos ls /store/user/jhgoh/MuonMisID/20170227_1/TT_powheg/ntuple_*.root
```
(for HYU, `/data/users/jhgoh/MuonMisID/20170227_1/`)

new root files with histograms will be created using run_selector.C
A root file containing invariant mass distributions for each pt, eta bins will be created. This use multiple CPUs using the Proof.

```
root -b -q -l run_selector.C
```

Run the fitter to extract misID probability.

```
python fit.py ks
python fit.py phi
python fit.py lamb
```

This will hist_ks.root and produce fit_ks.root, etc. You can find fit canvases in the output root files.

Plots can be produced with the final step,

```python -i draw.py pion```

or `python -i draw.py kaon`, `python -i draw.py lambda`.

## RPCMuon efficiency measurement using the Tag and Probe method

Test the ntuple production

```
cd MuonAnalysis/MuonIdentification/test/TagAndProbe
cmsRun tp_prod_MC_cfg.py
root -l tnp.root
```

Submit crab jobs

```
source /cvmfs/cms.cern.ch/crab3/crab.sh

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
