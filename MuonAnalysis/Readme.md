= Muon misID measurement =

설치 방법

cmsrel CMSSW_7_6_3_patch2
cd CMSSW_7_6_3_patch2/src
cmsenv
git-cms-init
git clone https://github.com/cms-analysis/MuonAnalysis 
git clone https://github.com/hep-skku/SKKU
scram b -j8

테스팅 하기

cd CMSSW_7_6_3_patch2/src/SKKU/MuonAnalysis/test
cmsRun prod_misID_MC_cfg.py
root -l hist.root

TTbar powheg sample과 JetHT 2015D dataset은 ntuple production완료해서 eos에 저장해 두었음.
eos ls /store/user/jhgoh/MuonMisID/20150401_2/JetHT_2015D/JetHT/crab_20160401_143746/160401_123805/0000
eos ls /store/user/jhgoh/MuonMisID/20150401_2/TT_powheg/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_20160401_143809/160401_123825/0000

crab job섭밋하는 법

source /cvmfs/cms.cern.ch/crab3/crab.sh
crab submit submitCrabRD.py
crab submit submitCrabMC.py

현재 crab cfg파일은 JetHT와 TTbar powheg에 대해 세팅 되어 있음. 
