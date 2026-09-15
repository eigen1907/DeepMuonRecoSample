#!/bin/bash
set -euo pipefail

events="$1"
skip_events="$2"
input_file="$3"
output_file="$4"

: "${DMR_CMSSW_BASE:?Missing DMR_CMSSW_BASE}"
: "${_CONDOR_SCRATCH_DIR:?Not running in a Condor scratch directory}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${DMR_CMSSW_BASE}/src"
eval "$(scram runtime -sh)"

if [[ "${CMSSW_VERSION}" != "CMSSW_14_0_21_patch1" ]]; then
    echo "Expected CMSSW_14_0_21_patch1, got ${CMSSW_VERSION}" >&2
    exit 2
fi

cd "${_CONDOR_SCRATCH_DIR}"
cmsRun "${DMR_CMSSW_BASE}/src/DeepMuonRecoSample/Run3/test/run_ntuple_cfg.py" \
    inputFiles="${input_file}" maxEvents="${events}" \
    skipEvents="${skip_events}" isMC=false outputFile="${output_file}"
