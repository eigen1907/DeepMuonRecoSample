#!/bin/bash
set -euo pipefail

cfg="$1"
input_file="$2"
max_events="$3"

: "${DMR_CMSSW_BASE:?Missing DMR_CMSSW_BASE}"
: "${_CONDOR_SCRATCH_DIR:?Not running in a Condor scratch directory}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${DMR_CMSSW_BASE}/src"
eval "$(scram runtime -sh)"
cd "${_CONDOR_SCRATCH_DIR}"

cmsRun "${DMR_CMSSW_BASE}/src/DeepMuonRecoSample/Phase2/test/${cfg}" \
    inputFiles="file:${input_file}" outputFile=output.root maxEvents="${max_events}"
