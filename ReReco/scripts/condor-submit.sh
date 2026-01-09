#!/bin/bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate deepmuonreco-py312
cmsenv
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(ls -d /cvmfs/cms.cern.ch/$(scram arch)/external/cuda/*/lib64/stubs | tail -n 1)

python condor-submit.py \
    -p ../test/runReReco_cfg.py \
    -i /hdfs/store/mc/Phase2Spring24DIGIRECOMiniAOD/SingleMu_FlatPt-2to100/GEN-SIM-DIGI-RAW-MINIAOD/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2 \
    -o /hdfs/user/joshin/DeepMuonReco/CMSSW_14_0_9/SingleMu_FlatPt-2to100/PU200_Trk1GeV_140X_mcRun4_realistic_v4-v2/RERECO/v02
