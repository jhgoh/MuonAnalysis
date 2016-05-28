import FWCore.ParameterSet.Config as cms

fitTemplate = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    NumCPU = cms.uint32(1),
    SaveWorkspace = cms.bool(False),

    Variables = cms.PSet(
        mass = cms.vstring("Tag-muon Mass", "76", "125", "GeV/c^{2}"),
        pt = cms.vstring("muon p_{T}", "0", "1000", "GeV/c"),
        eta    = cms.vstring("muon #eta", "-2.5", "2.5", ""),
        abseta = cms.vstring("muon |#eta|", "0", "2.5", ""),
        tag_nVertices = cms.vstring("Number of vertices", "0", "999", ""),
    ),

    Categories = cms.PSet(
        Glb   = cms.vstring("Global", "dummy[pass=1,fail=0]"),
        PF    = cms.vstring("PF Muon", "dummy[pass=1,fail=0]"),
        TM    = cms.vstring("Tracker Muon", "dummy[pass=1,fail=0]"),
        RPC   = cms.vstring("RPC Muon", "dummy[pass=1,fail=0]"),
    ),

    Expressions = cms.PSet(
        LooseDefault = cms.vstring("Loose Default", "PF&&(Glb||TM)", "PF", "Glb", "TM"),
        LooseWithRPC = cms.vstring('Loose with RPC', "PF&&(Glb||TM||RPC)", "PF", "Glb", "TM", "RPC"),
    ),

    Cuts = cms.PSet(
        cutLooseDefault = cms.vstring("cut Loose Default", "LooseDefault", "0.5"),
        cutLooseWithRPC = cms.vstring("cut Loose with RPC", "LooseWithRPC", "0.5"),
    ),

    PDFs = cms.PSet(
        voigtPlusExpo = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
            "Exponential::backgroundPass(mass, lp[0,-5,5])",
            "Exponential::backgroundFail(mass, lf[0,-5,5])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),

        voigtPlusPol = cms.vstring(
            "Voigtian::signal(mass, mean[90,80,100], width[2.495], sigma[3,1,20])",
            "RooChebychev::backgroundPass(mass, {pp0[0, -5, 5], pp1[0, -5, 5]})",
            "RooChebychev::backgroundFail(mass, {pf0[0, -5, 5], pf1[0, -5, 5]})",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]"
        ),
    ),

    binnedFit = cms.bool(True),
    binsForFit = cms.uint32(40),
    saveDistributionsPlot = cms.bool(False),

    InputFileNames = cms.vstring(),
    InputTreeName = cms.string("fitter_tree"),
    InputDirectoryName = cms.string("tpTree"),
    OutputFileName = cms.string('tnp_fit.root'),
    Efficiencies = cms.PSet(
        RPC_voigtExpo = cms.PSet(
            BinToPDFmap = cms.vstring('voigtPlusExpo'),
            BinnedVariables = cms.PSet(
                #pt     = cms.vdouble(  10, 20, 30, 40, 60, 100 ),
                pt     = cms.vdouble( 20, 250 ),
                #abseta = cms.vdouble(  0.0, 0.4, 0.9, 1.2, 1.6, 1.8, 1.9, 2.1, 2.4),
                abseta = cms.vdouble([x*0.1 for x in range(0, 17)])
            ),
            EfficiencyCategoryAndState = cms.vstring('RPC', 'pass'),
            UnbinnedVariables = cms.vstring("mass"),
        ),
#        RPC_voigtPol = cms.PSet(
#            BinToPDFmap = cms.vstring('voigtPlusPol'),
#            BinnedVariables = cms.PSet(
#                pt     = cms.vdouble(  10, 20, 30, 40, 60, 100 ),
#                abseta = cms.vdouble(  0.0, 0.4, 0.9, 1.2, 1.6, 1.8, 1.9, 2.1, 2.4),
#            ),
#            EfficiencyCategoryAndState = cms.vstring('RPC', 'pass'),
#            UnbinnedVariables = cms.vstring("mass"),
#        ),
    ),
)
