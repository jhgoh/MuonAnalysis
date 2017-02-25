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
Testing it

```
cd CMSSW_7_6_3_patch2/src/SKKU/MuonAnalysis/test
cmsRun prod_misID_MC_cfg.py
root -l hist.root
```

You can find out TTbar powheg sample and JetHT 2015D ntuple in the eos

```
eos ls /store/user/jhgoh/MuonMisID/20150401_2/JetHT_2015D/JetHT/crab_20160401_143746/160401_123805/0000
eos ls /store/user/jhgoh/MuonMisID/20150401_2/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160401_143809/160401_123825/0000
```

Submit crab jobs

```
crab submit submitCrabRD.py
crab submit submitCrabMC.py
```

Currently the crab cfg files are set for the JetHT and ttbar powheg samples

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
