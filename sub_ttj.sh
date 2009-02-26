#! /bin/bash
cd /home/cms/sckao/Top/CMSSW_2_2_3/src/PhysicsTools/TtAnalysis/
source /home/cms/sckao/.bashrc
eval `scramv1 runtime -sh`
cmsRun ttjets.py
