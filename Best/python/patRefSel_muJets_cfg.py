#
# This file contains the Top PAG reference selection work-flow for mu + jets analysis.
# as defined in
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopLeptonPlusJetsRefSel_mu#Selection_Version_SelV4_valid_fr
#

import sys

import FWCore.ParameterSet.Config as cms

# setup 'standard' options
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
options.register('runOnMC', False, VarParsing.VarParsing.multiplicity.singleton, VarParsing.VarParsing.varType.bool, "decide if run on MC or data")

process = cms.Process( 'PAT' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


### Data or MC?
runOnMC = options.runOnMC
#runOnMC = False

### Switch on/off selection steps

# Step 0a
useTrigger      = True
# Step 0b
useGoodVertex   = True
# Step 1
useGoodMuon     = True
# Step 2
useMuonVeto     = True
# Step 3
useElectronVeto = True
# Step 4a
use1Jet         = True
# Step 4b
use2Jets        = True
# Step 4c (choice depends on trigger)
use3JetsLoose   = True 
use3JetsTight   = False 
# Step 5
use4Jets        = False
## Step 6
#useBTag         = False

### Trigger matching?
addTriggerMatching = True

### Reference selection

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets import *
# Muons general
#muonsUsePV     = False
#muonEmbedTrack = True
# Muons
#muonCut       = ''
#signalMuonCut = ''
#muonVertexMaxDZ = 0.5
# Electrons
#electronCut = ''
# Jets
#jetCut          = ''
#jetCutPF        = 'pt > 20 && abs(eta) < 2.5 && numberOfDaughters > 1 && neutralHadronEnergyFraction < 0.99 && neutralEmEnergyFraction < 0.99 && (chargedEmEnergyFraction < 0.99 || abs(eta) > 2.4) && (chargedHadronEnergyFraction > 0. || abs(eta) >= 2.4) && (chargedMultiplicity > 0 || abs(eta) >= 2.4)'

#veryLooseJetCut = 'pt > 20.' # transverse momentum (all jets)
#looseJetCut     = 'pt > 35.' # transverse momentum (3rd jet, optional for 'use3JetsLoose = True')
#tightJetCut     = 'pt > 45.' # transverse momentum (leading jets)

# Trigger and trigger object
#triggerSelectionData       = 'HLT_IsoMu20_eta2p1_TriCentralPFJet30_v* OR HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v*'
#triggerSelectionData       = 'HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v* OR HLT_IsoMu17_eta2p1_TriCentralPFJet30_v* '
#triggerSelectionData       = 'HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_v* OR HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet30_30_20_v* OR HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet45_35_25_v* OR HLT_IsoMu17_eta2p1_TriCentralPFJet30_v* OR HLT_IsoMu20_eta2p1_TriCentralPFJet30_v* OR HLT_IsoMu20_eta2p1_TriCentralPFNoPUJet30_v* OR HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet50_40_30_v*'
triggerSelectionData       = 'HLT_IsoMu17_eta2p1_TriCentral* OR HLT_Mu17_eta2p1_TriCentral*'
triggerObjectSelectionData        = 'type("TriggerMuon") && ( path("HLT_IsoMu17_eta2p1_TriCentral*") || path("HLT_Mu17_eta2p1_TriCentral*") )'
#triggerSelectionData       = 'HLT_IsoMu24_eta2p1_v*'
#triggerObjectSelectionData = 'HLT_IsoMu24_eta2p1_v*'
#triggerSelectionMC       = ''
#triggerObjectSelectionMC = ''

### Particle flow

postfix = 'PF'


# subtract charged hadronic pile-up particles (from wrong PVs)
# effects also JECs
usePFnoPU       = True # before any top projection
usePfIsoLessCHS = True # switch to new PF isolation with L1Fastjet CHS

# other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = True # before MET top projection

# cuts used in top projections
# vertices
#pfVertices  = 'goodOfflinePrimaryVertices'
#pfD0Cut     = 0.2
#pfDzCut     = 0.5
# muons
#pfMuonSelectionCut = 'pt > 5.'
useMuonCutBasePF = False # use minimal (veto) muon selection cut on top of 'pfMuonSelectionCut'
#pfMuonIsoConeR03 = False
#pfMuonCombIsoCut = 0.2
# electrons
#pfElectronSelectionCut  = 'pt > 5. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'
useElectronCutBasePF  = False # use minimal (veto) electron selection cut on top of 'pfElectronSelectionCut'
#pfElectronIsoConeR03 = True
#pfElectronCombIsoCut  = 0.2

### JEC levels

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True  # takes effect only on data
useL5Flavor     = False
useL7Parton     = False

typeIMetCorrections = True

### Input

# list of input files
useRelVals = False # if 'False', "inputFiles" is used
inputFiles = [#'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/197/044/E2D96041-00BE-E111-8BC1-001D09F2447F.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/EC2E6A5A-56AF-E111-9081-001D09F2924F.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/E8E3D3F9-55AF-E111-A81F-0025901D624E.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/E284461A-5AAF-E111-881B-001D09F25267.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/DC4C4117-55AF-E111-BE40-001D09F290BF.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/DAF425A5-55AF-E111-AE7B-E0CB4E4408C4.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/D8BD9075-4EAF-E111-AF3F-001D09F28EA3.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/D6777FC4-53AF-E111-8247-0025B32034EA.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/D0E96041-61AF-E111-9351-BCAEC5329705.root'
             #,'/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/195/398/CC6B2BF4-56AF-E111-9D79-001D09F25109.root'
              '/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/752/3C0A66AC-789B-E111-94AB-003048CF99BA.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/774/88364725-889B-E111-B053-001D09F27003.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/806/ACB6E731-C89B-E111-811E-001D09F24664.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/812/145D80B4-D09B-E111-84E4-001D09F23A20.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/818/885D618A-E69B-E111-BBE8-001D09F25041.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/829/3AE91D27-049C-E111-A34A-BCAEC518FF8E.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/830/BA69838D-059C-E111-8DE1-003048D2BC30.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/834/6A191E05-5C9C-E111-8744-5404A63886A8.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/835/12C53245-4F9C-E111-B615-001D09F24FBA.root'
             ,'/store/data/Run2012B/MuHad/AOD/PromptReco-v1/000/193/836/54FABAD3-3E9C-E111-A345-003048D3C90E.root'
             # '/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/770/FCB29466-77C3-E111-8D0C-001D09F2462D.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/772/F65787B4-79C3-E111-971B-0019B9F72F97.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/774/56A19EC6-79C3-E111-884F-001D09F2462D.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/885/E875AEF2-2BC4-E111-8B0E-003048F11C58.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/889/A21ED86E-2DC4-E111-824B-00237DDC5BBC.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/891/94D9A5C2-2CC4-E111-A545-001D09F2906A.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/903/2A0AEA62-EAC4-E111-847B-5404A640A642.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/197/931/5278717C-D1C4-E111-9FCC-003048D3733E.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/198/011/3293599A-6BC5-E111-9FA4-BCAEC5364C62.root'
             #,'/store/data/Run2012C/MuHad/AOD/PromptReco-v1/000/198/022/9018F0EF-DCC5-E111-B098-0025901D6268.root' 
             ] # overwritten, if "useRelVals" is 'True'

# maximum number of events
maxEvents = -1 # reduce for testing

### Conditions

# GlobalTags
globalTagData = 'GR_R_52_V7D::All' # incl. Summer12 JEC and new b-tag SF
globalTagMC   = 'START52_V9C::All' # incl. Summer12 JEC and new b-tag SF

### Output

# output file
#outputFile = 'patRefSel_muJets.root'
outputFile = 'patTuple.root'

# event frequency of Fwk report
fwkReportEvery = 1000

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True


###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###


###
### Basic configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_basics_cff" )
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options.wantSummary = wantSummary
if runOnMC:
  process.GlobalTag.globaltag = globalTagMC
else:
  process.GlobalTag.globaltag = globalTagData


###
### Input configuration
###

if useRelVals:
  from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles
  if runOnMC:
    inputFiles = pickRelValInputFiles( cmsswVersion  = 'CMSSW_5_2_5_cand1'
                                     , relVal        = 'RelValTTbar'
                                     , globalTag     = 'START52_V9'
                                     , maxVersions   = 1
                                     )
  else:
    inputFiles = pickRelValInputFiles( cmsswVersion  = 'CMSSW_5_2_5_cand1'
                                     , relVal        = 'SingleMu'
                                     , dataTier      = 'RECO'
                                     , globalTag     = 'GR_R_52_V7_RelVal_mu2011B'
                                     , maxVersions   = 1
                                     )
process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
process.source.fileNames = inputFiles
process.maxEvents.input  = maxEvents


###
### Output configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
# output file name
process.out.fileName = outputFile
# event content
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out.outputCommands += patEventContent
# clear event selection
process.out.SelectEvents.SelectEvents = []


###
### Cleaning and trigger selection configuration
###

### Trigger selection
if runOnMC:
  triggerSelection = triggerSelectionMC
else:
  if useRelVals:
    triggerSelection = triggerSelectionDataRelVals
  else:
    triggerSelection = triggerSelectionData
from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults
process.step0a = triggerResults.clone(
  triggerConditions = [ triggerSelection ]
)

### Good vertex selection
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
process.step0b = process.goodOfflinePrimaryVertices.clone( filter = True )

### Event cleaning
process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
process.trackingFailureFilter.VertexSource = cms.InputTag( pfVertices )
process.step0c = process.eventCleaning
if runOnMC:
  process.step0c += process.eventCleaningMC
else:
  process.step0c += process.eventCleaningData


###
### PAT/PF2PAT configuration
###

process.load( "PhysicsTools.PatAlgos.patSequences_cff" )

### Check JECs

# JEC set
jecSet = 'AK5PF'
if usePFnoPU:
  jecSet += 'chs'

# JEC levels
if useL1FastJet and useL1Offset:
  sys.exit( 'ERROR: switch off either "L1FastJet" or "L1Offset"' )
jecLevels = []
if useL1FastJet:
  jecLevels.append( 'L1FastJet' )
if useL1Offset:
  jecLevels.append( 'L1Offset' )
if useL2Relative:
  jecLevels.append( 'L2Relative' )
if useL3Absolute:
  jecLevels.append( 'L3Absolute' )
if useL2L3Residual and not runOnMC:
  jecLevels.append( 'L2L3Residual' )
if useL5Flavor:
  jecLevels.append( 'L5Flavor' )
if useL7Parton:
  jecLevels.append( 'L7Parton' )

### Switch configuration

from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
usePF2PAT( process
         , runPF2PAT           = True
         , runOnMC             = runOnMC
         , jetAlgo             = jetAlgo
         , postfix             = postfix
         , jetCorrections      = ( jecSet
                                 , jecLevels
                                 )
         , typeIMetCorrections = typeIMetCorrections
         , pvCollection        = cms.InputTag( pfVertices )
         )

if useMuonCutBasePF:
  pfMuonSelectionCut += ' && %s'%( muonCut )
if useElectronCutBasePF:
  pfElectronSelectionCut += ' && %s'%( electronCut )

getattr( process, 'pfNoPileUp'   + postfix ).enable = usePFnoPU
getattr( process, 'pfNoMuon'     + postfix ).enable = useNoMuon
getattr( process, 'pfNoElectron' + postfix ).enable = useNoElectron
getattr( process, 'pfNoJet'      + postfix ).enable = useNoJet
getattr( process, 'pfNoTau'      + postfix ).enable = useNoTau

if useL1FastJet:
  getattr( process, 'pfPileUpIso' + postfix ).checkClosestZVertex = usePfIsoLessCHS

getattr( process, 'pfMuonsFromVertex' + postfix ).d0Cut = pfD0Cut
getattr( process, 'pfMuonsFromVertex' + postfix ).dzCut = pfDzCut
getattr( process, 'pfSelectedMuons'   + postfix ).cut = pfMuonSelectionCut
getattr( process, 'pfIsolatedMuons'   + postfix ).isolationCut = pfMuonCombIsoCut

if pfMuonIsoConeR03:
  getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix )
                                                                                            )
  getattr( process, 'pfIsolatedMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
  getattr( process, 'pfIsolatedMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                            , cms.InputTag( 'muPFIsoValueGamma03' + postfix )
                                                                                            )
  getattr( process, 'pfMuons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix )
                                                                                    )
  getattr( process, 'pfMuons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
  getattr( process, 'pfMuons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                    , cms.InputTag( 'muPFIsoValueGamma03' + postfix )
                                                                                    )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'muPFIsoValueChargedAll03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'muPFIsoValueGamma03' + postfix )
  getattr( process, 'patMuons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'muPFIsoValueCharged03' + postfix )

getattr( process, 'pfElectronsFromVertex' + postfix ).d0Cut = pfD0Cut
getattr( process, 'pfElectronsFromVertex' + postfix ).dzCut = pfDzCut
getattr( process, 'pfSelectedElectrons'   + postfix ).cut = pfElectronSelectionCut
getattr( process, 'pfIsolatedElectrons'   + postfix ).isolationCut = pfElectronCombIsoCut

if pfElectronIsoConeR03:
  getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
                                                                                                )
  getattr( process, 'pfIsolatedElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
  getattr( process, 'pfIsolatedElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                                , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix )
                                                                                                )
  getattr( process, 'pfElectrons' + postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
                                                                                        )
  getattr( process, 'pfElectrons' + postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
  getattr( process, 'pfElectrons' + postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                        , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix )
                                                                                        )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedAll       = cms.InputTag( 'elPFIsoValueChargedAll03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
  getattr( process, 'patElectrons' + postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )


from PhysicsTools.PatAlgos.tools.coreTools import *

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import *

# remove MC matching, object cleaning, objects etc.
if not runOnMC:
  runOnData( process
           , names = [ 'PFAll' ]
           , postfix = postfix
           )
removeSpecificPATObjects( process
                        , names = [ 'Photons', 'Taus' ]
                        , postfix = postfix
                        ) # includes 'removeCleaning'

# additional event content has to be (re-)added _after_ the call to 'removeCleaning()':
process.out.outputCommands += [ 'keep edmTriggerResults_*_*_*'
                              , 'keep *_hltTriggerSummaryAOD_*_*'
                              # vertices and beam spot
                              , 'keep *_offlineBeamSpot_*_*'
                              , 'keep *_offlinePrimaryVertices*_*_*'
                              , 'keep *_goodOfflinePrimaryVertices*_*_*'
                              ]
if runOnMC:
  process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                , 'keep recoGenParticles_*_*_*'
                                , 'keep *_addPileupInfo_*_*'
                                ]


###
### Additional configuration
###

### Muons

intermediatePatMuons.src = cms.InputTag( 'selectedPatMuons' + postfix )
setattr( process, 'intermediatePatMuons' + postfix, intermediatePatMuons )

goodPatMuons.muonSource   = cms.InputTag( 'intermediatePatMuons' + postfix )
goodPatMuons.vertexSource = cms.InputTag( pfVertices )
setattr( process, 'goodPatMuons' + postfix, goodPatMuons )

step1.src = cms.InputTag( 'goodPatMuons' + postfix )
setattr( process, 'step1' + postfix, step1 )
step2.src = cms.InputTag( 'selectedPatMuons' + postfix )
setattr( process, 'step2' + postfix, step2 )

### Jets

veryLoosePatJets.src = cms.InputTag( 'selectedPatJets' + postfix )
veryLoosePatJets.cut = veryLooseJetCut
setattr( process, 'veryLoosePatJets' + postfix, veryLoosePatJets )
loosePatJets.src = cms.InputTag( 'veryLoosePatJets' + postfix )
loosePatJets.cut = looseJetCut
setattr( process, 'loosePatJets' + postfix, loosePatJets )
tightPatJets.src = cms.InputTag( 'loosePatJets' + postfix )
tightPatJets.cut = tightJetCut
setattr( process, 'tightPatJets' + postfix, tightPatJets )

step4a.src = cms.InputTag( 'tightPatJets' + postfix )
setattr( process, 'step4a' + postfix, step4a )
step4b.src = cms.InputTag( 'tightPatJets' + postfix )
setattr( process, 'step4b' + postfix, step4b )
step4cTight.src = cms.InputTag( 'tightPatJets' + postfix )
setattr( process, 'step4cTight' + postfix, step4cTight )
step4cLoose.src = cms.InputTag( 'loosePatJets' + postfix )
setattr( process, 'step4cLoose' + postfix, step4cLoose )
step5.src = cms.InputTag( 'veryLoosePatJets' + postfix )
setattr( process, 'step5' + postfix, step5  )

### Electrons

step3.src = cms.InputTag( 'selectedPatElectrons' + postfix )
setattr( process, 'step3' + postfix, step3 )

process.out.outputCommands.append( 'keep *_goodPatMuons*_*_*' )
process.out.outputCommands.append( 'keep *_veryLoosePatJets*_*_*' )
process.out.outputCommands.append( 'keep *_loosePatJets*_*_*' )
process.out.outputCommands.append( 'keep *_tightPatJets*_*_*' )


###
### Selection configuration
###

### Muons

getattr( process, 'patMuons' + postfix ).usePV      = muonsUsePV
getattr( process, 'patMuons' + postfix ).embedTrack = muonEmbedTrack

getattr( process, 'selectedPatMuons' + postfix ).cut = muonCut

getattr( process, 'intermediatePatMuons' + postfix ).cut = signalMuonCut

getattr( process, 'goodPatMuons' + postfix ).maxDZ = muonVertexMaxDZ

### Jets

getattr( process, 'selectedPatJets'  + postfix ).cut = jetCut
getattr( process, 'veryLoosePatJets' + postfix ).cut = veryLooseJetCut
getattr( process, 'loosePatJets'     + postfix ).cut = looseJetCut
getattr( process, 'tightPatJets'     + postfix ).cut = tightJetCut

### Electrons

getattr( process, 'patElectrons' + postfix ).electronIDSources = electronIDSources

getattr( process, 'selectedPatElectrons' + postfix ).cut = electronCut


###
### Trigger matching
###

if addTriggerMatching:

  if runOnMC:
    triggerObjectSelection = triggerObjectSelectionMC
  else:
    if useRelVals:
      triggerObjectSelection = triggerObjectSelectionDataRelVals
    else:
      triggerObjectSelection = triggerObjectSelectionData
  ### Trigger matching configuration
  from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
  from TopQuarkAnalysis.Configuration.patRefSel_triggerMatching_cfi import patMuonTriggerMatch
  from PhysicsTools.PatAlgos.tools.trigTools import *
  triggerProducerPF = patTrigger.clone()
  setattr( process, 'patTrigger' + postfix, triggerProducerPF )
  triggerMatchPF = patMuonTriggerMatch.clone( matchedCuts = triggerObjectSelection )
  setattr( process, 'triggerMatch' + postfix, triggerMatchPF )
  switchOnTriggerMatchEmbedding( process
                               , triggerProducer = 'patTrigger' + postfix
                               , triggerMatchers = [ 'triggerMatch' + postfix ]
                               , sequence        = 'patPF2PATSequence' + postfix
                               , postfix         = postfix
                               )
  removeCleaningFromTriggerMatching( process
                                   , sequence = 'patPF2PATSequence' + postfix
                                   )
  getattr( process, 'intermediatePatMuons' + postfix ).src = cms.InputTag( 'selectedPatMuons' + postfix + 'TriggerMatch' )


###
### Scheduling
###

# MVA electron ID

process.load( "EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi" )
process.eidMVASequence = cms.Sequence(
  process.mvaTrigV0
+ process.mvaNonTrigV0
)

# The additional sequence

patAddOnSequence = cms.Sequence(
  getattr( process, 'intermediatePatMuons' + postfix )
* getattr( process, 'goodPatMuons'         + postfix )
* getattr( process, 'veryLoosePatJets'     + postfix )
* getattr( process, 'loosePatJets'         + postfix )
* getattr( process, 'tightPatJets'         + postfix )
)
setattr( process, 'patAddOnSequence' + postfix, patAddOnSequence )

# The paths

process.p = cms.Path()
if useTrigger:
  process.p += process.step0a
process.p += process.goodOfflinePrimaryVertices
if useGoodVertex:
  process.p += process.step0b
process.p += process.step0c
process.p += process.eidMVASequence
process.p += getattr( process, 'patPF2PATSequence' + postfix )
process.p += getattr( process, 'patAddOnSequence' + postfix )
if useGoodMuon:
  process.p += getattr( process, 'step1' + postfix )
if useMuonVeto:
  process.p += getattr( process, 'step2' + postfix )
if useElectronVeto:
  process.p += getattr( process, 'step3' + postfix )
if use1Jet:
  process.p += getattr( process, 'step4a' + postfix )
if use2Jets:
  process.p += getattr( process, 'step4b' + postfix )
if use3JetsTight:
  process.p += getattr( process, 'step4cTight' + postfix )
elif use3JetsLoose:
  process.p += getattr( process, 'step4cLoose' + postfix )
if use4Jets:
  process.p += getattr( process, 'step5' + postfix )
process.out.SelectEvents.SelectEvents.append( 'p' )

##______________________________________________________________________________________________//
### PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process ) # overwrite sequence default "patDefaultSequence", since it is not used in any path
process.patTrigger.addL1Algos = cms.bool(True) # add L1 algorithms' collection
switchOnTrigger( process ) # called once more to update the event content according to the changed parameters!!!

process.hltTrigReport = cms.EDAnalyzer( "HLTrigReport",
    HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT' ),
    reportBy = cms.untracked.string( 'run' ),
    resetBy  = cms.untracked.string( 'run' ),
)

process.MessageLogger.categories.append("HLTrigReport")
process.outPath = cms.EndPath(
    process.hltTrigReport
)
process.p += process.patTrigger

##______________________________________________________________________________________________//
process.out.outputCommands = cms.untracked.vstring(
  'drop *',
  'keep *Vertex*_goodOfflinePrimaryVertices_*_PAT',
  'keep *_goodPatMuons*_*_*',
  'keep *_selectedPatMuonsPF_*_*',
  'keep *Jet*_veryLoosePatJets*_*_*',
  'keep *Jet*_loosePatJets*_*_*',
  'keep *Jet*_tightPatJets*_*_*', 
  'keep *_patMETs*_*_*',
  'keep *TriggerPath*_patTrigger_*_*',
  'keep *_hltTriggerSummaryAOD_*_*',
  'keep *_TriggerResults_*_*',
  'keep *_offlineBeamSpot_*_*',
  'keep *_offlinePrimaryVertices*_*_*',
  )

##______________________________________________________________________________________________//


