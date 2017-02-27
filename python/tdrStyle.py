# tdrStyle
import os
from ROOT import gROOT

gROOT.ProcessLine(".L %s/src/PhysicsTools/TagAndProbe/test/utilities/tdrstyle.C" % os.environ["CMSSW_RELEASE_BASE"])
