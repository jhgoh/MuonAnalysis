=== RPCMuon efficiency measurement using the Tag and Probe method ===

설치 방법
```
cmsrel CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src
cmsenv
git-cms-init
git clone https://github.com/cms-analysis/MuonAnalysis-TagAndProbe 
git clone https://github.com/hep-skku/SKKU
scram b -j8
```

ntuple production 테스팅 하기

```
cd CMSSW_7_6_3_patch2/src/SKKU/RPCMuonAnalysis/test/TagAndProbe
cmsRun tp_prod_MC_cfg.py
root -l tnp.root
```

crab job submit하기

```
source /cvmfs/cms.cern.ch/crab3/crab.sh
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
