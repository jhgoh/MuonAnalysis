---++ Installation
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

---++ Muon misID measurement


테스팅 하기

```
cd CMSSW_7_6_3_patch2/src/SKKU/MuonAnalysis/test
cmsRun prod_misID_MC_cfg.py
root -l hist.root
```

TTbar powheg sample과 JetHT 2015D dataset은 ntuple production완료해서 eos에 저장해 두었음.
```
eos ls /store/user/jhgoh/MuonMisID/20150401_2/JetHT_2015D/JetHT/crab_20160401_143746/160401_123805/0000
eos ls /store/user/jhgoh/MuonMisID/20150401_2/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160401_143809/160401_123825/0000
```

crab job섭밋하는 법

```
crab submit submitCrabRD.py
crab submit submitCrabMC.py
```

현재 crab cfg파일은 JetHT와 TTbar powheg에 대해 세팅 되어 있음. 
---++ RPCMuon efficiency measurement using the Tag and Probe method

ntuple production 테스팅 하기

```
cd MuonAnalysis/MuonIdentification/test/TagAndProbe
cmsRun tp_prod_MC_cfg.py
root -l tnp.root
```

crab job submit하기

```
crab submit submitCrabRD.py
crab submit submitCrabMC.py
```

DYJets MG5와 Run2015D singleMuon dataset에 대해서 production 돌리는중 (task monitoring: SingleMuon, DYJets). 저장은 eos에 하는 중.
```
/store/user/jhgoh/RPCMuonEff/20150404_1/DYJets_MG/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_20160404_173038/160404_153124/0000:
/store/user/jhgoh/RPCMuonEff/20150404_1/SingleMuon_2015D/SingleMuon/crab_20160404_173352/160404_153427/0000
```
완료되기까지는 시간이 걸리지만 일단 일부만 모아서 fitting을 바로 시작할 수 있다.

tnp fitting을 위한 cfg file을 data, mc에 대해 각각 준비해 두었음.
```
cmsRun tp_fit_MC_cfg.py
cmsRun tp_fit_RD_cfg.py
```
