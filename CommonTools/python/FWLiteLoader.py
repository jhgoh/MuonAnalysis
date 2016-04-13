#!/usr/bin/env python

import ROOT
import os

#Set up FW Lite for automatic loading of CMS libraries
#and data formats.   As you may have other user-defined setup
#in your rootlogon.C, the CMS setup is executed only if the CMS
#environment is set up.

if 'CMSSW_BASE' in os.environ:
    print "Loading FW Lite setup."
    ROOT.gSystem.Load("libFWCoreFWLite.so")
    FWLiteEnabler.enable()
    ROOT.gSystem.Load("libDataFormatsFWLite.so")
    ROOT.gSystem.Load("libDataFormatsPatCandidates.so")

