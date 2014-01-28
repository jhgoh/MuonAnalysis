import FWCore.ParameterSet.Config as cms

fakeKshort = cms.EDAnalyzer("FakeMuonAnalyzer",
    muon = cms.InputTag("muons"),
    vetoMuon = cms.InputTag("tagMuons"),
    nMaxVetoMuon = cms.uint32(1),
    vetoCone = cms.double(0.5),
    #vertexCand = cms.InputTag("generalV0Candidates", "Kshort"),
    vertexCand = cms.InputTag("kshortVertex"),
    muonIds = cms.PSet(),
    vertexCut = cms.string(""),
    match = cms.string("matchByTrackRef"),
    #maxDR = cms.double(0.02),
    #maxDPt = cms.double(0.1),
    doTree = cms.bool(True),
    doHist = cms.bool(False),
    massMin = cms.double(0.43),
    massMax = cms.double(0.56),
)

fakeKshort.muonIds.all = cms.PSet(
    cut = cms.string(""),
    idSelection = cms.string("All"),
    arbitration = cms.string("NoArbitration"),
)

fakeKshort.muonIds.globalMuon = cms.PSet(
    cut = cms.string("isGlobalMuon"),
    idSelection = cms.string("AllGlobalMuons"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)
fakeKshort.muonIds.globalMuonMedium = cms.PSet(
    cut = cms.string(
        "isGlobalMuon"
        " && globalTrack().normalizedChi2()<10.0"
        " && innerTrack().hitPattern().numberOfValidPixelHits() > 0"
        " && track().hitPattern().trackerLayersWithMeasurement() > 5"
    ),
    idSelection = cms.string("AllGlobalMuons"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)
fakeKshort.muonIds.globalMuonTight = cms.PSet(
    cut = cms.string(
        "isGlobalMuon"
        " && globalTrack().normalizedChi2()<10.0"
        " && innerTrack().hitPattern().numberOfValidPixelHits() > 0"
        " && track().hitPattern().trackerLayersWithMeasurement() > 5"
        " && numberOfMatchedStations() > 1"
        " && globalTrack().hitPattern().numberOfValidMuonHits() > 0"
    ),
    idSelection = cms.string("AllGlobalMuons"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)

fakeKshort.muonIds.globalMuonPT = cms.PSet(
    cut = cms.string("isGlobalMuon"),
    idSelection = cms.string("GlobalMuonPromptTight"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)

fakeKshort.muonIds.globalMuonTest = cms.PSet(
    cut = cms.string("isPFMuon && (isGlobalMuon || isTrackerMuon) "),
    idSelection = cms.string("All"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)

fakeKshort.muonIds.TMOneStationLoose = cms.PSet(
    cut = cms.string("isTrackerMuon"),
    idSelection = cms.string("TMOneStationLoose"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)
fakeKshort.muonIds.TMOneStationTight = cms.PSet(
    cut = cms.string("isTrackerMuon"),
    idSelection = cms.string("TMOneStationTight"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)
fakeKshort.muonIds.TMTwoStationTest = cms.PSet(
    cut = cms.string("isTrackerMuon  && numberOfMatchedStations >= 2"),
    idSelection = cms.string("TrackerMuonArbitrated"),
    arbitration = cms.string("SegmentAndTrackArbitration"),
)

fakeKshort.muonIds.RPCMuLoose = cms.PSet(
    cut = cms.string("isRPCMuon"),
    idSelection = cms.string("RPCMuLoose"),
    arbitration = cms.string("RPCHitAndTrackArbitration"),
)
fakeKshort.muonIds.RPCMuMedium = cms.PSet(
    cut = cms.string("isRPCMuon &&  numberOfMatchedStations('RPCHitAndTrackArbitration') >= 2 "),
    idSelection = cms.string("RPCMuLoose"),
    arbitration = cms.string("RPCHitAndTrackArbitration"),
)
fakeKshort.muonIds.RPCMuTight = cms.PSet(
    cut = cms.string("isRPCMuon &&  numberOfMatchedStations('RPCHitAndTrackArbitration') >= 2 && numberOfMatchedRPCLayers >=3 "),
    idSelection = cms.string("RPCMuLoose"),
    arbitration = cms.string("RPCHitAndTrackArbitration"),
)


fakeLambda = fakeKshort.clone()
#fakeLambda.vertexCand = cms.InputTag("generalV0Candidates", "Lambda")
fakeLambda.vertexCand = cms.InputTag("lambdaVertex")
fakeLambda.massMin = cms.double(1.08)
fakeLambda.massMax = cms.double(1.20)

fakePhi = fakeKshort.clone()
fakePhi.vertexCand = cms.InputTag("phiVertex")
fakePhi.massMin = cms.double(1.00)
fakePhi.massMax = cms.double(1.04)

fakeJpsi = fakeKshort.clone()
fakeJpsi.vertexCand = cms.InputTag("jpsiVertex")
fakeJpsi.massMin = cms.double(2.80)
fakeJpsi.massMax = cms.double(3.40)
